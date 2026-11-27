############################################################################
# polyA_nc_group_compare.R
#
# Content-based triplet pairing + round-spermatid (RS) classification +
# WT vs cKO comparison of polyA RNA nucleo/cytoplasmic distribution.
#
# Input:  folders of single-stain RGB TIFF exports, arbitrary file order:
#           red channel   = polyA FISH
#           green channel = DDX4
#           blue channel  = DAPI
#         Triplets (one field of view) are found by image correlation,
#         NOT by file order.
#
# Usage:  Rscript polyA_nc_group_compare.R <out_dir> <WT_folder> [cKO_folder] ... [prob_dirs]
#         prob_dirs (optional, last argument): one or more folders with ilastik
#         probability files, separated by ";" (e.g. one per group). Files named
#         <group>_<DAPIbase>_Probabilities.h5 or <DAPIbase>_Probabilities.h5,
#         e.g. DDX4-WT_DDX4-0001_Probabilities.h5 (channel 1 = nucleus class).
#         Fields with a matching h5 use ilastik-based nucleus segmentation;
#         fields without one fall back to the adaptive-DAPI method (flagged
#         as seg_method in the output tables).
############################################################################

suppressPackageStartupMessages(library(EBImage))
# options(warn = 1)   # print warnings as they happen (Rscript exits before warnings() is possible)

args    <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else "/mnt/results/polyA_nc_quant/group_comparison"
prob_dirs <- NULL
if (length(args) >= 4) {
  prob_dirs <- strsplit(args[length(args)], ";", fixed = TRUE)[[1]]
  prob_dirs <- trimws(prob_dirs)
  args <- args[-length(args)]
}
folders <- if (length(args) >= 2) args[2:length(args)] else
  c("/mnt/user-uploads/DDX4_polyA/DDX4-WT", "/mnt/user-uploads/DDX4_polyA/DDX4-cKO")

PARAMS <- list(
  sigma_smooth  = 2,
  min_nuc_area  = 685,    # was 120: keep small RS nuclei 在设置大小前自己拿imageJ画一下，别拍脑袋写数字
  max_nuc_area  = 3000,
  seed_brush    = 5,      # fallback seed-erosion brush (used if peak seeds fail)
  # --- DAPI-dim nucleus recovery (adaptive threshold + intensity-peak seeds) ---
  adapt_w       = 51,     # adaptive-threshold window (local mean, px)
  adapt_offset  = 0.02,   # adaptive threshold: dapi_s > local_mean + offset
  peak_r        = 5,      # local-maximum radius for seed detection (px)
  peak_min_dapi = 0.30,   # absolute DAPI floor: rejects background noise peaks
  nuc_erode     = 0,      # 核内信号取样区向内收缩（px）
  annulus_mode  = "ilastik_cyto", # "ddx4" = DDX4 ring; "ilastik_cyto" = class-2 cytoplasm (needs prob h5)
  cyto_class_id = 2,        # h5 channel of the cytoplasm class
  min_ddx4_mean = 0.055,    # lenient per-cell DDX4 gate for ilastik_cyto mode
  max_ring_px   = 50,     # was 30: perinuclear band width for cytoplasm sampling
                          # AND halo display (both modes): 用imageJ量了下，大概是30px左右吧，最大就是50了
  min_elong_ar  = 1.9,
  min_circ      = 1.9,
  min_ddx4_frac = 0.08,   # was 0.15: accept weaker DDX4 rings in dense clusters
  min_ann_pix   = 10,     # was 30: accept smaller cytoplasm samples
  min_ann_signal = 0.009,  # ~2 gray levels: cytoplasm must exceed background
  min_field_cells = 10,    # fields below this germ-cell count fail QC
  # round-spermatid criteria (user-confirmed: small + round + DDX4+)
  rs_area_max   = 1850,
  rs_ar_max     = 1.40,   # was 1.55: stricter roundness (major/minor axis)
  rs_circ_max   = 1.30,   # was 1.4: stricter circularity (1 = perfect circle)
  # field-level nuclear-polyA accumulation metric (segmentation-free)
  bright_polya_thr = 0.60,  # polyA intensity (0-1) counting as "bright"
  bright_min_area  = 180,   # px: ignore speckles
  bright_max_area  = 4000,  # px: keep dense unsplittable clusters
  # --- nuclear polyA speckles: bright polyA puncta INSIDE nuclei ---
  speckle_nuc_fold = 1.2, # speckle threshold = this x mean polyA of the SAME
                          # nucleus (per-nucleus adaptive; no absolute threshold)
  speckle_min_area = 70,    # px: minimum speckle size (nuclear puncta are small)
  # --- ilastik probability-based nucleus segmentation (preferred) ---
  prob_nuc_class = 1,     # h5 channel of the nucleus class (user-confirmed)
  prob_thr      = 0.5,    # nucleus probability threshold
  require_prob  = FALSE,  # TRUE = stop if any field lacks its probability h5
  # content-pairing thresholds
  pair_rg_min   = 0.20,     # polyA <-> DDX4
  pair_br_min   = 0.15,     # DAPI <-> polyA
  sat_thresh    = 254
)


# ------------------------- low-level helpers ------------------------------
read_stain <- function(path) {
  img <- readImage(path)
  if (length(dim(img)) == 2) return(img)
  rng <- apply(img, 3, function(ch) diff(quantile(ch, c(0.01, 0.999))))
  out <- img[, , which.max(rng)]
  colorMode(out) <- Grayscale
  out
}
otsu_mask <- function(ch, sigma) { sm <- gblur(ch, sigma); fillHull(sm > otsu(sm)) }
maxfilt <- function(x, r) {   # local-maximum filter (square window 2r+1)
  nr <- nrow(x); nc <- ncol(x)
  pad <- matrix(-1e6, nr + 2*r, nc + 2*r)
  pad[(r+1):(nr+r), (r+1):(nc+r)] <- x
  out <- x
  for (dx in -r:r) for (dy in -r:r)
    out <- pmax(out, pad[(r+1+dx):(nr+r+dx), (r+1+dy):(nc+r+dy)])
  out
}
# ---- ilastik probability I/O ----
have_rhdf5 <- requireNamespace("rhdf5", quietly = TRUE)
find_prob_file <- function(prob_dirs, grp, dapi_path) {
  dirn <- dirname(dapi_path)
  base <- sub("\\.tiff?$", "", basename(dapi_path), ignore.case = TRUE)
  # "DDX4-WT-1-DAPI" -> "DDX4-WT-1" (drop trailing stain tag)
  stem <- sub("[-_ ]?dapi$", "", base, ignore.case = TRUE)
  cands <- c(file.path(dirn, paste0(base, "_Probabilities.h5")),
             file.path(dirn, paste0(stem, "_Probabilities.h5")),
             file.path(dirn, paste0(stem, "_merged_Probabilities.h5")))
  if (!is.null(prob_dirs))
    for (d in prob_dirs)
      cands <- c(file.path(d, paste0(grp, "_", base, "_Probabilities.h5")),
                 file.path(d, paste0(grp, "_", stem, "_Probabilities.h5")),
                 file.path(d, paste0(base, "_Probabilities.h5")),
                 file.path(d, paste0(stem, "_Probabilities.h5")),
                 file.path(d, paste0(stem, "_merged_Probabilities.h5")), cands)
  hit <- cands[file.exists(cands)]
  if (length(hit)) hit[1] else NA_character_
}
read_prob_nuc <- function(path, ref, cls) {
  a <- rhdf5::h5read(path, "exported_data")
  if (length(dim(a)) == 3) {           # locate the class axis, then extract it
    d <- dim(a); caxis <- which(d == min(d))[1]
    a <- if (caxis == 1) a[cls, , ] else if (caxis == 2) a[, cls, ] else a[, , cls]
  }
  a <- as.matrix(a)
  if (!all(dim(a) == dim(ref))) {      # EBImage orientation: rows = x
    if (all(dim(t(a)) == dim(ref))) a <- t(a) else
      stop("probability map size ", paste(dim(a), collapse = "x"),
           " does not match image ", paste(dim(ref), collapse = "x"))
  }
  Image(a)
}
clean_mask <- function(m, min_area, max_area) {
  lab <- bwlabel(m)
  if (max(lab, 0) == 0) return(lab)
  sz <- table(lab); ids <- as.integer(names(sz))
  keep <- ids[ids > 0 & sz >= min_area & sz <= max_area]
  lab[!(lab %in% keep)] <- 0L
  bwlabel(lab > 0)
}
nuc_features <- function(lab) {
  f  <- computeFeatures.shape(lab)
  fm <- computeFeatures.moment(lab)
  area  <- as.numeric(f[, "s.area"])
  perim <- if ("s.perimeter" %in% colnames(f)) as.numeric(f[, "s.perimeter"]) else NA
  major <- pmax(as.numeric(fm[, "m.majoraxis"]), 1e-6)
  ecc   <- pmin(as.numeric(fm[, "m.eccentricity"]), 0.999999)
  data.frame(id = seq_len(max(lab, 0)), area = area, perim = perim,
             ar = major / pmax(major * sqrt(1 - ecc^2), 1e-6),
             circ = perim^2 / (4 * pi * pmax(area, 1)),
             x = as.numeric(fm[, "m.cx"]), y = as.numeric(fm[, "m.cy"]))
}
ds <- function(x, f = 4) x[seq(1, dim(x)[1], f), seq(1, dim(x)[2], f)]
corr_ab <- function(a, b) {
  if (sd(a) < 1e-9 || sd(b) < 1e-9) return(0)
  a <- a - mean(a); b <- b - mean(b)
  as.numeric(mean(a * b) / (sd(a) * sd(b) + 1e-9))
}

# ------------------------- content-based pairing --------------------------
# classify each file by dominant channel, then greedily match R<->G<->B
pair_folder <- function(folder, P) {
  if (!dir.exists(folder))
    stop("Group folder does not exist (check the path for stray quote characters):\n  ",
         folder, call. = FALSE)
  ff <- list.files(folder, full.names = TRUE)
  ff <- ff[grepl("\\.(tif|tiff|png|jpg)$", ff, ignore.case = TRUE)]
  if (!length(ff))
    stop("No image files (.tif/.tiff/.png/.jpg) found directly inside\n  ", folder,
         "\nIf the images are in a subfolder, pass the subfolder path instead.",
         call. = FALSE)
  imgs   <- list()
  stains <- character(length(ff))
  for (i in seq_along(ff)) {
    img <- readImage(ff[i])
    if (length(dim(img)) == 3) {
      ch <- apply(img, 3, mean); k <- which.max(ch)
      imgs[[i]] <- ds(img[, , k])
    } else {           # single-channel grayscale: cannot tell the stain by color
      ch <- c(mean(img), 0, 0); k <- 1L
      imgs[[i]] <- ds(img)
    }
    stains[i] <- if (max(ch) < 5 / 255) "?" else c("R", "G", "B")[k]
  }
  info <- data.frame(file = ff, stain = stains, stringsAsFactors = FALSE)
  rs <- info[info$stain == "R", ]; gs <- info[info$stain == "G", ]
  bs <- info[info$stain == "B", ]; other <- info[info$stain == "?", ]
  im_r <- imgs[info$stain == "R"]; im_g <- imgs[info$stain == "G"]; im_b <- imgs[info$stain == "B"]
  if (nrow(other)) for (f in other$file) message("UNCLASSIFIED (skipped): ", basename(f))
  if (!nrow(rs) || !nrow(gs) || !nrow(bs))
    stop(sprintf(paste0(
      "Pairing needs polyA(R), DDX4(G) and DAPI(B) single-stain files in every folder, ",
      "but in\n  %s\nfound R=%d, G=%d, B=%d. If the files are multi-page merged TIFFs or ",
      "single-channel grayscale images the stain cannot be told apart by color - ",
      "point the script at the original single-stain RGB export folders instead."),
      folder, nrow(rs), nrow(gs), nrow(bs)), call. = FALSE)

  used_g <- logical(nrow(gs))
  pairs_rg <- NULL
  cand <- expand.grid(ri = seq_len(nrow(rs)), gi = seq_len(nrow(gs)))
  cand$rho <- as.numeric(mapply(function(ri, gi)
    corr_ab(im_r[[ri]], im_g[[gi]]), cand$ri, cand$gi))
  cand <- cand[order(-cand$rho), ]
  used_r <- logical(nrow(rs))
  for (k in seq_len(nrow(cand))) {
    if (cand$rho[k] < P$pair_rg_min) break
    ri <- cand$ri[k]; gi <- cand$gi[k]
    if (used_r[ri] || used_g[gi]) next
    used_r[ri] <- TRUE; used_g[gi] <- TRUE
    pairs_rg <- rbind(pairs_rg, data.frame(r_file = rs$file[ri], g_file = gs$file[gi],
                                           rho_rg = cand$rho[k], ri_hint = ri))
  }
  unpaired <- c(rs$file[!used_r], gs$file[!used_g])
  for (f in unpaired) message("UNPAIRED (R/G, skipped): ", basename(f))

  trips <- NULL
  if (!is.null(pairs_rg)) {
    used_b <- logical(nrow(bs))
    cand2 <- expand.grid(pi = seq_len(nrow(pairs_rg)), bi = seq_len(nrow(bs)))
    cand2$rho <- as.numeric(mapply(function(pi, bi)
      corr_ab(im_r[[pairs_rg$ri_hint[pi]]], im_b[[bi]]), cand2$pi, cand2$bi))
    cand2 <- cand2[order(-cand2$rho), ]
    used_p <- logical(nrow(pairs_rg))
    for (k in seq_len(nrow(cand2))) {
      if (cand2$rho[k] < P$pair_br_min) break
      pi <- cand2$pi[k]; bi <- cand2$bi[k]
      if (used_p[pi] || used_b[bi]) next
      used_p[pi] <- TRUE; used_b[bi] <- TRUE
      trips <- rbind(trips, data.frame(
        dapi = bs$file[bi], polya = pairs_rg$r_file[pi],
        ddx4 = pairs_rg$g_file[pi], rho_rg = pairs_rg$rho_rg[pi],
        rho_br = cand2$rho[k]))
    }
    for (pi in which(!used_p)) message("UNPAIRED (B, skipped): ", basename(pairs_rg$r_file[pi]))
    for (bi in which(!used_b)) message("UNPAIRED (B, skipped): ", basename(bs$file[bi]))
  }
  trips
}
read_stain_data <- function(path) ds(read_stain(path))

# ------------------------- per-field quantification -----------------------
process_set <- function(tag, files, out_dir, P) {
  dapi  <- read_stain(files$dapi);  polya <- read_stain(files$polya)
  ddx4  <- read_stain(files$ddx4)
  hw <- dim(dapi); h <- hw[1]; w <- hw[2]

  use_prob <- !is.null(files$prob) && !is.na(files$prob) && file.exists(files$prob)
  if (use_prob && !have_rhdf5)
    stop("rhdf5 not installed; run BiocManager::install(\"rhdf5\") to read probability files")
  seg_method <- if (use_prob) "ilastik_prob" else "dapi_adaptive"

  if (use_prob) {
    # ilastik nucleus-class probability: continuous surface splits dense clusters
    nuc_p  <- read_prob_nuc(files$prob, dapi, P$prob_nuc_class)
    prob_s <- gblur(nuc_p, 1)
    nuc_mask <- fillHull(prob_s > P$prob_thr)
    nuc_lab <- clean_mask(nuc_mask, P$min_nuc_area, 10^6)
    d <- distmap(nuc_lab)
    # seeds = local maxima of the (continuous) nucleus probability surface
    set.seed(1)
    pn <- prob_s + array(runif(length(prob_s), 0, 1e-4), dim(prob_s))
    seeds <- bwlabel(pn >= maxfilt(pn, P$peak_r) & nuc_lab > 0 & prob_s > P$prob_thr)
    if (max(seeds, 0) == 0)
      seeds <- bwlabel(erode(nuc_lab > 0, makeBrush(P$seed_brush, "disc")))
    nuc_lab <- propagate(d, seeds, mask = nuc_lab > 0)
    nuc_lab <- clean_mask(nuc_lab > 0, P$min_nuc_area, P$max_nuc_area)
  } else {
    # adaptive threshold catches DAPI-dim / polyA-bright nuclei missed by global Otsu
    dapi_s <- gblur(dapi, P$sigma_smooth)
    nuc_mask <- fillHull(thresh(dapi_s, w = P$adapt_w, offset = P$adapt_offset) & dapi_s > P$peak_min_dapi)
    # pre-clean: drop only tiny noise; oversized clumps are kept so peak seeds can split them
    nuc_lab <- clean_mask(nuc_mask, P$min_nuc_area, 10^6)
    d <- distmap(nuc_lab)
    # seeds = local maxima of smoothed DAPI inside the mask (splits touching nuclei)
    seeds <- bwlabel(dapi_s >= maxfilt(dapi_s, P$peak_r) & nuc_lab > 0 & dapi_s > P$peak_min_dapi)
    if (max(seeds, 0) == 0)
      seeds <- bwlabel(erode(nuc_lab > 0, makeBrush(P$seed_brush, "disc")))
    nuc_lab <- propagate(d, seeds, mask = nuc_lab > 0)
    nuc_lab <- clean_mask(nuc_lab > 0, P$min_nuc_area, P$max_nuc_area)
  }
  border_ids <- unique(c(nuc_lab[1, ], nuc_lab[h, ], nuc_lab[, 1], nuc_lab[, w]))
  border_ids <- setdiff(border_ids, 0)
  if (length(border_ids)) nuc_lab <- rmObjects(nuc_lab, border_ids)

  nf <- nuc_features(nuc_lab)
  if (nrow(nf) == 0) { message("No nuclei for ", tag); return(NULL) }

  nf$elongated <- nf$ar >= P$min_elong_ar | nf$circ >= P$min_circ
  keep_ids <- nf$id[!nf$elongated]; drop_ids <- nf$id[nf$elongated]
  keep_mask <- nuc_lab * (nuc_lab %in% keep_ids) > 0

  ddx4_mask <- otsu_mask(ddx4, P$sigma_smooth)
  ddx4_mask <- erode(dilate(ddx4_mask, makeBrush(5, "disc")), makeBrush(5, "disc"))

  # cytoplasm definition: DDX4 ring (default) or ilastik cytoplasm class
  use_cyto2 <- identical(P$annulus_mode, "ilastik_cyto") && use_prob
  if (use_cyto2) {
    cyto2_p  <- read_prob_nuc(files$prob, dapi, P$cyto_class_id)
    cyto2_m  <- fillHull(gblur(cyto2_p, 1) > P$prob_thr)
    bg_region <- (nuc_lab == 0) & !cyto2_m
  } else bg_region <- (nuc_lab == 0) & !ddx4_mask
  bg_polya  <- median(polya[bg_region])

  nuc_er <- erode(keep_mask, makeBrush(2 * P$nuc_erode + 1, "box")) > 0
  seeds  <- nuc_lab * (nuc_lab %in% keep_ids)
  cap    <- dilate(keep_mask, makeBrush(2 * P$max_ring_px + 1, "disc")) > 0
  if (use_cyto2) {
    # learned cytoplasm restricted to the perinuclear band, Voronoi-assigned
    vor  <- propagate(Image(1, hw), seeds = seeds, mask = (cyto2_m & cap) | keep_mask)
    vor[!(cyto2_m & cap) | keep_mask] <- 0L
  } else {
    ann_mask <- ddx4_mask & cap & !keep_mask
    vor  <- propagate(Image(1, hw), seeds = seeds, mask = ann_mask | keep_mask)
    vor[!ann_mask] <- 0L
  }

  er_lab <- nuc_lab * nuc_er; er_lab[!(er_lab %in% keep_ids)] <- 0L
  ids <- keep_ids
  res <- data.frame(set = tag, id = ids)
  res$seg_method <- seg_method
  res$area_px <- nf$area[match(ids, nf$id)]
  res$x <- nf$x[match(ids, nf$id)]; res$y <- nf$y[match(ids, nf$id)]

  tmean <- function(img, lab) { v <- img[lab > 0]; k <- lab[lab > 0]
    m <- tapply(v, k, mean); as.numeric(m[as.character(ids)]) }
  tnpix <- function(lab) { t <- table(lab[lab > 0]); as.numeric(t[as.character(ids)]) }

  res$nuc_polya_raw <- tmean(polya, er_lab)
  res$ann_polya_raw <- tmean(polya, vor)
  res$nuc_npix <- tnpix(er_lab); res$nuc_npix[is.na(res$nuc_npix)] <- 0
  res$ann_npix <- tnpix(vor);    res$ann_npix[is.na(res$ann_npix)] <- 0
  res$ann_ddx4_frac <- tnpix(vor * ddx4_mask) / pmax(res$ann_npix, 1)
  res$ann_ddx4_frac[is.na(res$ann_ddx4_frac)] <- 0
  res$ann_ddx4_mean <- tmean(ddx4, vor); res$ann_ddx4_mean[is.na(res$ann_ddx4_mean)] <- 0
  res$ring_width_px <- res$ann_npix / pmax(nf$perim[match(ids, nf$id)], 1)

  res$nuc_polya_bc <- res$nuc_polya_raw - bg_polya
  res$ann_polya_bc <- res$ann_polya_raw - bg_polya
  res$nc_ratio <- res$nuc_polya_bc / pmax(res$ann_polya_bc, 1e-6)
  # germ-cell gate: DDX4 fraction of ring (ddx4 mode) or mean DDX4 intensity of
  # the assigned cytoplasm (ilastik_cyto mode, lenient: keeps DDX4-weak germ cells)
  germ_gate <- if (use_cyto2) res$ann_ddx4_mean >= P$min_ddx4_mean
               else res$ann_ddx4_frac >= P$min_ddx4_frac
  res$germ_cell <- res$ann_npix >= P$min_ann_pix &
                   germ_gate &
                   !is.na(res$nuc_polya_bc) & !is.na(res$ann_polya_bc) &
                   is.finite(res$nc_ratio)
  res$is_round_spermatid <- res$germ_cell &
                   res$area_px < P$rs_area_max &
                   nf$ar[match(ids, nf$id)] < P$rs_ar_max &
                   nf$circ[match(ids, nf$id)] < P$rs_circ_max
  # per-cell shape metrics (px-based) so RS thresholds can be tuned from the CSV
  res$ar   <- nf$ar[match(ids, nf$id)]
  res$circ <- nf$circ[match(ids, nf$id)]

  # --- nuclear polyA speckles: bright polyA puncta inside kept nuclei ---
  # a "speckle" = connected region of polyA >= (speckle_nuc_fold x mean polyA
  # of that same nucleus) within the nucleus mask. The per-nucleus relative
  # threshold compensates for field-to-field brightness differences.
  # NOTE: uniformly bright nuclei have NO pixel above 1.2 x their own mean, so
  # they score 0 speckles by design - their accumulation shows in
  # frac_nuc_polya_bright / nuc_polya_bc instead.
  sp_thr_map <- P$speckle_nuc_fold * res$nuc_polya_raw[match(as.vector(seeds), ids)]
  sp_thr_map[is.na(sp_thr_map)] <- Inf
  sp_mask <- matrix(as.vector(polya) >= sp_thr_map,
                    dim(polya)[1], dim(polya)[2]) & keep_mask
  message(sprintf("  speckle threshold: %.2f x mean nuclear polyA (per nucleus)",
                  P$speckle_nuc_fold))
  sp_lab  <- bwlabel(sp_mask > 0)
  if (max(sp_lab, 0) > 0) {
    sz <- table(sp_lab); s_ids <- as.integer(names(sz))
    s_keep <- s_ids[s_ids > 0 & sz >= P$speckle_min_area]
    sp_lab[!(sp_lab %in% s_keep)] <- 0L
  }
  sp_lab <- bwlabel(sp_lab > 0)
  n_speckles <- max(sp_lab, 0)
  # assign each speckle to its (majority-overlap) nucleus, count per cell
  res$n_speckles <- 0L
  res$speckle_npix <- 0    # total speckle pixels per cell (sum of speckle areas)
  res$speckle_intensity <- NA_real_  # mean polyA inside speckles of that nucleus
  if (n_speckles > 0) {
    idx <- which(as.vector(sp_lab) > 0)
    spv <- as.vector(sp_lab)[idx]; nuv <- as.vector(nuc_lab)[idx]
    tb  <- table(spv, nuv)
    nuc_of_sp <- as.integer(colnames(tb)[apply(tb, 1, which.max)])
    cnt <- tapply(nuc_of_sp, nuc_of_sp, length)
    hit <- as.integer(names(cnt))
    res$n_speckles[match(hit, ids)] <- as.integer(cnt)
    szs <- table(spv)                       # pixel count of each speckle
    npix <- tapply(as.numeric(szs), nuc_of_sp[as.integer(names(szs))], sum)
    hit2 <- as.integer(names(npix))
    res$speckle_npix[match(hit2, ids)] <- as.numeric(npix)
    # mean polyA intensity inside each nucleus's speckle pixels (raw scale);
    # expand per-speckle nucleus assignment to per-pixel before grouping
    nu_of_px <- rep(nuc_of_sp, as.integer(table(spv)))
    spi <- tapply(as.vector(polya)[idx], nu_of_px, mean)
    hit3 <- as.integer(names(spi))
    res$speckle_intensity[match(hit3, ids)] <- as.numeric(spi)
  }
  # per-cell speckle load: fraction of eroded-nucleus pixels above threshold
  res$nuc_speckle_frac <- tmean(Image(sp_mask), er_lab)
  res$nuc_speckle_frac[is.na(res$nuc_speckle_frac)] <- 0

  # QC overlay: cyan=kept nuclei, red=elongated excluded, orange=round spermatids
  col3 <- rgbImage(dapi, polya, ddx4)
  if (length(drop_ids))
    col3 <- paintObjects(bwlabel(nuc_lab * (nuc_lab %in% drop_ids) > 0),
                         col3, col = "red", thick = 1)
  rs_ids <- res$id[res$is_round_spermatid]
  if (length(rs_ids))
    col3 <- paintObjects(bwlabel(nuc_lab * (nuc_lab %in% rs_ids) > 0),
                         col3, col = "orange", thick = 2)
  col3 <- paintObjects(bwlabel(keep_mask), col3, col = "cyan", thick = 1)
  if (n_speckles > 0)
    col3 <- paintObjects(sp_lab, col3, col = "#75A025", thick = 1)  # green = nuclear polyA speckles
  ann_disp <- if (use_cyto2) cyto2_m & cap & !keep_mask else ann_mask
  for (i in 1:3) {
    plane <- col3[, , i]; am <- as.vector(ann_disp)
    plane[am] <- 0.75 * plane[am] + 0.25 * c(1, 1, 0)[i]
    col3[, , i] <- plane
  }
  writeImage(col3, file.path(out_dir, paste0(tag, "_QC_overlay.png")))

  # NC map: only ROUND SPERMATIDS are tinted by N/C (blue = N/C<=1, magenta =
  # N/C>1). Larger non-RS germ cells (spermatocyte-stage) are left unpainted
  # so they cannot be mistaken for RS on the map.
  ncmap <- as.array(rgbImage(dapi, polya, ddx4))
  ncd <- res[res$germ_cell & res$is_round_spermatid, ]
  for (i in seq_len(nrow(ncd))) {
    sel <- as.vector(nuc_lab == ncd$id[i])
    tint <- if (isTRUE(ncd$nc_ratio[i] > 1)) c(0, 0.45, 1) else c(1, 0.35, 0.95)
    for (ch in 1:3) {
      plane <- ncmap[, , ch]
      plane[sel] <- 0.35 * plane[sel] + 0.65 * tint[ch]
      ncmap[, , ch] <- plane
    }
  }
  writeImage(rgbImage(ncmap[, , 1], ncmap[, , 2], ncmap[, , 3]),
             file.path(out_dir, paste0(tag, "_NC_map.png")))

  list(per_cell = res, summary = data.frame(
    set = tag,
    seg_method = seg_method,
    n_nuclei = nrow(nf),
    n_elongated_excluded = length(drop_ids),
    n_germ_cells = sum(res$germ_cell),
    n_round_spermatids = sum(res$is_round_spermatid),
    median_NC_RS = if (sum(res$is_round_spermatid)) median(res$nc_ratio[res$is_round_spermatid]) else NA,
    mean_NC_RS   = if (sum(res$is_round_spermatid)) mean(res$nc_ratio[res$is_round_spermatid]) else NA,
    median_NC_allround = if (sum(res$germ_cell)) median(res$nc_ratio[res$germ_cell]) else NA,
    # --- field-level "nuclear polyA accumulation" metric (segmentation-free) ---
    # captures the focal phenotype: tight clusters of small nuclei with strong
    # nuclear polyA that cannot be split by single-plane DAPI.
    n_bright_nuc = {
      m <- fillHull(gblur(polya, 2) > P$bright_polya_thr & gblur(dapi, 2) > P$peak_min_dapi)
      lab <- bwlabel(m)
      if (max(lab, 0) > 0) {
        sz <- table(lab); ids <- as.integer(names(sz))
        keep <- ids[ids > 0 & sz >= P$bright_min_area & sz <= P$bright_max_area]
        lab[!(lab %in% keep)] <- 0L
      }
      max(bwlabel(lab > 0), 0)
    },
    frac_nuc_polya_bright = mean(polya[keep_mask] > P$bright_polya_thr),
    # --- nuclear speckle metrics ---
    n_speckles_total = n_speckles,
    n_cells_with_speckles = sum(res$n_speckles > 0),
    frac_cells_with_speckles = mean(res$n_speckles[res$germ_cell] > 0),
    median_speckles_per_positive_cell =
      if (any(res$n_speckles > 0)) median(res$n_speckles[res$n_speckles > 0]) else 0L,
    mean_speckle_area_px =
      if (n_speckles > 0) mean(table(sp_lab[sp_lab > 0])) else 0,
    median_speckle_intensity = if (any(!is.na(res$speckle_intensity)))
      median(res$speckle_intensity, na.rm = TRUE) else NA,
    bg_polya = bg_polya))
}

# ------------------------------- main --------------------------------------
SCRIPT_VERSION <- "v14-speckfix"   # bump when output columns/definitions change:
                                 # a checkpoint written by another version is discarded
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# archive the exact parameters of this run (overwritten each run)
write.csv(data.frame(param = names(PARAMS), value = unlist(PARAMS)),
          file.path(out_dir, "params_used.csv"), row.names = FALSE)
ck_file <- paste0(out_dir, "/shared-workspace/polyA_nc_group_checkpoint.rds")
dir.create(paste0(out_dir, "/shared-workspace"), showWarnings = FALSE, recursive = TRUE)
if (file.exists(ck_file)) {
  ck <- readRDS(ck_file)
  if (isTRUE(ck$version == SCRIPT_VERSION) && isTRUE(identical(ck$params, PARAMS))) {
    all_pc <- ck$all_pc; all_su <- ck$all_su; pairlog <- ck$pairlog
    message("Resuming from checkpoint: ", length(all_su), " fields already done.")
  } else {
    all_pc <- list(); all_su <- list(); pairlog <- NULL
    message("Checkpoint was written by a different script version - starting fresh.")
  }
} else { all_pc <- list(); all_su <- list(); pairlog <- NULL }

for (fld in folders) {
  grp <- basename(fld)
  done_here <- sum(grepl(paste0("^", grp, "_"), names(all_su)))
  if (done_here) message("Group ", grp, ": ", done_here, " field(s) already done, processing the rest.")
  message("==== Group: ", grp, " ====")
  trips <- pair_folder(fld, PARAMS)
  if (is.null(trips) || nrow(trips) == 0) {
    message("No complete triplets in ", grp); next
  }
  message("Paired ", nrow(trips), " triplets in ", grp)
  trips$group <- grp
  pairlog <- rbind(pairlog, transform(trips, set = sub("\\.tif$", "", basename(dapi))))
  for (k in seq_len(nrow(trips))) {
    tag <- paste0(grp, "_", sub("\\.tif$", "", basename(trips$dapi[k])))
    if (tag %in% names(all_su)) next   # field already in checkpoint
    prob_f <- find_prob_file(prob_dirs, grp, trips$dapi[k])
    if (is.na(prob_f)) {
      if (PARAMS$require_prob)
        stop("No probability h5 for ", tag, " (require_prob = TRUE)")
      message(tag, ": no probability h5 -> adaptive-DAPI segmentation")
    } else message(tag, ": ilastik probability segmentation (", basename(prob_f), ")")
    files_k <- list(dapi = trips$dapi[k], polya = trips$polya[k],
                    ddx4 = trips$ddx4[k], prob = prob_f)
    out <- tryCatch(process_set(tag, files_k, out_dir, PARAMS),
                    error = function(e) { message("ERROR in ", tag, ": ",
                      conditionMessage(e)); NULL })
    if (!is.null(out)) {
      all_pc[[tag]] <- out$per_cell; all_su[[tag]] <- out$summary
      saveRDS(list(version = SCRIPT_VERSION, params = PARAMS, all_pc = all_pc,
                   all_su = all_su, pairlog = pairlog), ck_file)
    }
  }
}

if (!length(all_su)) stop("No triplets processed.")

per_cell <- do.call(rbind, all_pc)
summ     <- do.call(rbind, all_su)
# group = the folder name the field came from (set = <group>_<field>).
# Match against the actual group-folder names instead of splitting on "_",
# which breaks when the group name itself contains underscores (URH_WT_merge).
set2group <- function(sets) {
  pref <- paste0(basename(folders), "_")
  out <- rep(NA_character_, length(sets))
  for (i in seq_along(sets)) {
    hit <- pref[startsWith(sets[i], pref)]
    if (length(hit)) out[i] <- sub("_$", "", hit[which.max(nchar(hit))])
  }
  out[is.na(out)] <- sub("_.*$", "", sets[is.na(out)])  # legacy naming fallback
  out
}
summ$group <- set2group(summ$set)
per_cell$group <- set2group(per_cell$set)
if (!"seg_method" %in% names(summ)) summ$seg_method <- NA_character_
if (!"seg_method" %in% names(per_cell)) per_cell$seg_method <- NA_character_

# re-apply validity filters on the (possibly checkpointed) per-cell table:
# cytoplasm signal must rise above background, otherwise the N/C ratio is
# undefined (division by ~0) and the cell is excluded
valid <- per_cell$germ_cell & per_cell$ann_polya_bc >= PARAMS$min_ann_signal
n_dropped <- sum(per_cell$germ_cell, na.rm = TRUE) - sum(valid, na.rm = TRUE)
message("Excluded ", n_dropped, " germ cells with cytoplasm signal at background level.")
per_cell$germ_cell <- valid
per_cell$is_round_spermatid <- per_cell$is_round_spermatid & valid

# recompute per-field summaries from the filtered table
mk <- function(d) data.frame(
  set = d$set[1], group = d$group[1],
  n_germ_cells = sum(d$germ_cell),
  n_round_spermatids = sum(d$is_round_spermatid),
  median_NC_RS = if (sum(d$is_round_spermatid)) median(d$nc_ratio[d$is_round_spermatid]) else NA,
  mean_NC_RS   = if (sum(d$is_round_spermatid)) mean(d$nc_ratio[d$is_round_spermatid]) else NA,
  median_NC_allround = if (sum(d$germ_cell)) median(d$nc_ratio[d$germ_cell]) else NA)
summ_rs <- do.call(rbind, lapply(split(per_cell, per_cell$set), mk))
summ$median_NC_RS <- summ_rs[match(summ$set, summ_rs$set), "median_NC_RS"]
summ$mean_NC_RS   <- summ_rs[match(summ$set, summ_rs$set), "mean_NC_RS"]
summ$n_round_spermatids <- summ_rs[match(summ$set, summ_rs$set), "n_round_spermatids"]
summ$n_germ_cells <- summ_rs[match(summ$set, summ_rs$set), "n_germ_cells"]
summ$median_NC_allround <- summ_rs[match(summ$set, summ_rs$set), "median_NC_allround"]
summ$field_pass <- summ$n_germ_cells >= PARAMS$min_field_cells
for (s in summ$set[!summ$field_pass])
  message("FIELD FAILS QC (low germ-cell count, excluded from stats): ", s)
if (length(unique(stats::na.omit(summ$seg_method))) > 1)
  message("WARNING: fields segmented with MIXED methods (",
          paste(unique(summ$seg_method), collapse = " / "),
          "); group stats mix both. Supply probability h5 for all fields for a clean run.")
write.csv(per_cell, file.path(out_dir, "per_cell_quantification.csv"), row.names = FALSE)
write.csv(summ,     file.path(out_dir, "summary_per_image.csv"), row.names = FALSE)
write.csv(pairlog,  file.path(out_dir, "pairing_log.csv"), row.names = FALSE)

# --------------------- WT vs cKO comparison (RS) --------------------------
rs_med <- tapply(summ$median_NC_RS[summ$field_pass], summ$group[summ$field_pass],
                 function(x) x[!is.na(x)])
wt_name  <- grep("WT", names(rs_med), value = TRUE)[1]
ko_name  <- grep("KO", names(rs_med), value = TRUE)[1]
stats <- NULL
if (!is.na(wt_name) && !is.na(ko_name)) {
  a <- rs_med[[ko_name]]   # cKO first: diff = cKO - WT
  b <- rs_med[[wt_name]]
  tt <- t.test(a, b); wt <- wilcox.test(a, b)
  stats <- data.frame(
    comparison = paste(ko_name, "vs", wt_name),
    n1 = length(a), n2 = length(b),
    median_cKO = median(a), median_WT = median(b),
    mean_diff_cKO_minus_WT = mean(a) - mean(b),
    ci_lo = tt$conf.int[1], ci_hi = tt$conf.int[2],
    p_ttest = tt$p.value, p_wilcoxon = wt$p.value)
  write.csv(stats, file.path(out_dir, "group_comparison_stats.csv"), row.names = FALSE)
  print(stats)

  # --- field-level nuclear polyA accumulation (segmentation-free) ---
  sfp <- summ[summ$field_pass, ]
  bn_a <- sfp$n_bright_nuc[sfp$group == ko_name]
  bn_b <- sfp$n_bright_nuc[sfp$group == wt_name]
  if (length(bn_a) > 2 && length(bn_b) > 2) {
    tt2 <- t.test(bn_a, bn_b); wt2 <- suppressWarnings(wilcox.test(bn_a, bn_b))
    stats2 <- data.frame(
      comparison = paste(ko_name, "vs", wt_name, "(nuclear polyA bright objects per field)"),
      n1 = length(bn_a), n2 = length(bn_b),
      median_cKO = median(bn_a), median_WT = median(bn_b),
      mean_diff_cKO_minus_WT = mean(bn_a) - mean(bn_b),
      ci_lo = tt2$conf.int[1], ci_hi = tt2$conf.int[2],
      p_ttest = tt2$p.value, p_wilcoxon = wt2$p.value)
    write.csv(stats2, file.path(out_dir, "group_comparison_bright_nuclei.csv"), row.names = FALSE)
    print(stats2)
  }

  # --- stats restricted to ilastik-segmented fields (when both groups have >=3) ---
  sfp_il <- summ[summ$field_pass & summ$seg_method == "ilastik_prob" &
                 !is.na(summ$median_NC_RS), ]
  if (length(unique(sfp_il$group)) == 2 && all(table(sfp_il$group) >= 3)) {
    a3 <- sfp_il$median_NC_RS[sfp_il$group == ko_name]
    b3 <- sfp_il$median_NC_RS[sfp_il$group == wt_name]
    tt3 <- t.test(a3, b3); wt3 <- suppressWarnings(wilcox.test(a3, b3))
    stats3 <- data.frame(
      comparison = paste(ko_name, "vs", wt_name, "(ilastik-segmented fields only)"),
      n1 = length(a3), n2 = length(b3),
      median_cKO = median(a3), median_WT = median(b3),
      mean_diff_cKO_minus_WT = mean(a3) - mean(b3),
      ci_lo = tt3$conf.int[1], ci_hi = tt3$conf.int[2],
      p_ttest = tt3$p.value, p_wilcoxon = wt3$p.value)
    write.csv(stats3, file.path(out_dir, "group_comparison_stats_ilastik_only.csv"),
              row.names = FALSE)
    print(stats3)
  }
}

# --------------------- figures (SVG) --------------------------------------
pmain <- if (!is.null(stats)) stats$p_ttest else NA
svg(file.path(out_dir, "NC_WT_vs_cKO_boxplot.svg"), width = 4.5, height = 5)
par(mar = c(5, 4.5, 2, 1), family = "Liberation Sans")
sfp <- summ[summ$field_pass & !is.na(summ$median_NC_RS), ]
bp <- boxplot(median_NC_RS ~ group, sfp, col = c("#0279EE55", "#FD9BED55"),
              border = "#000000", outpch = NA, ylab = "median polyA N/C ratio (round spermatids)",
              xlab = "", main = "")
btc <- tapply(sfp$median_NC_RS, sfp$group, identity)
for (i in seq_along(btc)) {
  j <- jitter(rep(i, length(btc[[i]])), amount = 0.06)
  points(j, btc[[i]], pch = 16, col = "#00000090")
}
if (!is.na(pmain)) {
  lab <- if (pmain < 0.001) "***" else if (pmain < 0.01) "**" else if (pmain < 0.05) "*" else "ns"
  text(1.5, max(sfp$median_NC_RS, na.rm = TRUE) * 1.06, lab, cex = 1.4)
}
dev.off()

svg(file.path(out_dir, "NC_RS_distribution_overlay.svg"), width = 5, height = 5)
par(mar = c(5, 4.5, 2, 1), family = "Liberation Sans")
rs_cells <- per_cell[per_cell$is_round_spermatid & is.finite(per_cell$nc_ratio), ]
xr <- quantile(rs_cells$nc_ratio, c(0.005, 0.995))
rs_trim <- rs_cells[rs_cells$nc_ratio >= xr[1] & rs_cells$nc_ratio <= xr[2], ]
cols <- c("#0279EE", "#FD9BED"); names(cols) <- sort(unique(rs_trim$group))
x0 <- rs_trim$nc_ratio[rs_trim$group == names(cols)[1]]
x1 <- rs_trim$nc_ratio[rs_trim$group == names(cols)[2]]
x0 <- x0[is.finite(x0)]; x1 <- x1[is.finite(x1)]
if (length(x0) < 10 || length(x1) < 10) {
  message("NOTE: skipping NC_RS_distribution_overlay.svg (too few RS cells: ",
          length(x0), " vs ", length(x1), ")")
} else {
d0 <- density(x0); d1 <- density(x1)
plot(d0, col = cols[1], lwd = 2, main = "",
     xlab = "polyA N/C ratio (round spermatids)", ylab = "Density",
     xlim = range(c(d0$x, d1$x)), ylim = range(c(d0$y, d1$y)))
polygon(d0, col = paste0(cols[1], "40"), border = NA)
lines(d1, col = cols[2], lwd = 2); polygon(d1, col = paste0(cols[2], "40"), border = NA)
abline(v = 1, lty = 3, col = "#00000080")
legend("topright", bty = "n", fill = cols, legend = names(cols))
}
dev.off()

message("Done. Results in ", out_dir)
print(summ)
