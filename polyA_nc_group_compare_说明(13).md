# polyA_nc_group_compare.R 参数与输出说明

对 WT 与 URH49-cKO 小鼠睾丸 IF/FISH 三色图像（polyA FISH / DDX4 / DAPI）做批量单细胞核质（N/C）定量与组间比较。支持可选的 ilastik 概率图辅助分割。

## 运行方式

```
Rscript.exe polyA_nc_group_compare.R <输出目录> <WT文件夹> [cKO文件夹 ...] [概率h5文件夹]
```

- 参数一律用**位置参数**（不要写 `name=value`），路径含中文/空格时加引号
- 第 4 个参数（可选）：概率 h5 文件夹，**多个文件夹用分号 `;` 分隔**（如每组一个：`"H:\h5\WT;H:\h5\cKO"`）；查找顺序：`<组名>_<DAPI主干>_Probabilities.h5`、`<DAPI主干>_Probabilities.h5`、`<DAPI主干>_merged_Probabilities.h5`（主干 = DAPI 文件名去掉结尾 `-DAPI`/`_DAPI`，如 `DDX4-WT-1-DAPI` → `DDX4-WT-1`），也可直接把 h5 放在 DAPI 同目录（则无需第 4 参数）
- **依赖**：EBImage、rhdf5（`BiocManager::install(c("EBImage","rhdf5"))`）
- **改参数后必须删除旧 checkpoint**（输出目录外的 `polyA_nc_group_checkpoint.rds`），否则脚本沿用旧结果直接跳过

## 输入约定

- 每个数据文件夹 = 一只小鼠/一组；内含单色导出的 8 位 RGB TIFF：红=polyA，绿=DDX4，蓝=DAPI
- 三色文件按**图像内容相关性**自动配对（不依赖文件名顺序）
- 可选概率文件：命名含 DAPI 主干（`DDX4-WT-1-DAPI` → `DDX4-WT-1`）：`DDX4-WT-1_merged_Probabilities.h5` 或 `DDX4-WT-1_Probabilities.h5`，放在概率文件夹（可加组名前缀）或 DAPI 同目录均可；通道 1 = 细胞核、通道 2 = 胞质

## 参数（脚本开头 PARAMS 列表，第 32 行起）

### 分割

| 参数 | 默认 | 说明 |
|---|---|---|
| `sigma_smooth` | 2 | 高斯平滑 σ（像素） |
| `min_nuc_area` / `max_nuc_area` | 685 / 3000 | 核面积范围（px），超出丢弃 |
| `seed_brush` | 5 | 种子腐蚀刷尺寸（回退用） |
| `adapt_w` / `adapt_offset` | 51 / 0.02 | DAPI 自适应阈值：窗口 / 偏移（无概率文件时用） |
| `peak_r` | 5 | 分核种子的局部极大值半径（px） |
| `peak_min_dapi` | 0.30 | DAPI 绝对下限，抑制背景噪声峰 |
| `prob_nuc_class` / `prob_thr` | 1 / 0.5 | ilastik 核类别通道号 / 概率阈值 |
| `require_prob` | FALSE | TRUE 时缺概率文件的场直接报错（保证全组统一分割法） |
| `nuc_erode` | 0 | 核内信号取样区向内收缩（px） |

### 胞质与 N/C 计算

| 参数 | 默认 | 说明 |
|---|---|---|
| `annulus_mode` | "ilastik_cyto" | `"ddx4"`=DDX4⁺ 环形胞质；`"ilastik_cyto"`=ilastik 类别2 胞质 ∩ 核外 `max_ring_px` 带（需概率文件） |
| `cyto_class_id` | 2 | ilastik 胞质类别通道号 |
| `max_ring_px` | 50 | 胞质取样带距核边界最大宽度（px），也是 QC 图光晕宽度 |
| `min_ann_pix` | 10 | 胞质取样最少像素，不足剔除该细胞 |
| `min_ddx4_frac` | 0.08 | ddx4 模式：环内 DDX4⁺ 像素占比下限（生殖细胞门控） |
| `min_ddx4_mean` | 0.055 | ilastik_cyto 模式：分配胞质的 DDX4 平均强度下限（宽松门控） |
| `min_ann_signal` | 0.009 | 背景校正后胞质 polyA 下限，防止比值爆炸 |

### 圆形精子细胞（RS）判定（小 + 圆 + DDX4⁺）

| 参数 | 默认 | 说明 |
|---|---|---|
| `rs_area_max` | 1850 | 核面积上限（px） |
| `rs_ar_max` | 1.40 | 长短轴比上限 |
| `rs_circ_max` | 1.30 | 圆度上限（周长²/4π·面积） |

 elongated（精子头等长形核）排除：`min_elong_ar = 1.9` 或 `min_circ = 1.9`。

### 场级"核内 polyA 积聚"指标（不依赖单核分割，捕捉局灶表型）

| 参数 | 默认 | 说明 |
|---|---|---|
| `bright_polya_thr` | 0.60 | polyA"亮"强度阈值（0–1） |
| `bright_min_area` / `bright_max_area` | 180 / 4000 | 亮核物体面积范围（px） |

### 核内 polyA speckle（每个核内的亮斑点，逐核+逐场）

| 参数 | 默认 | 说明 |
|---|---|---|
| `speckle_nuc_fold` | 1.15 | 每核自适应阈值：polyA ≥ 该核自身平均强度 × 1.15 的像素记为斑点。天然消除场间/细胞间亮度差异。注意：整核均匀亮的核没有像素能超过自身均值 ×1.15，按定义 0 个斑点——这类核的积聚用 `frac_nuc_polya_bright`/`nuc_polya_bc` 衡量 |
| `speckle_min_area` | 70 | 单个 speckle 最小面积（px），过滤噪声点 |

定义：核掩膜内 polyA ≥ 阈值的连通区；**无面积上限**（整核均亮的细胞记为 1 个大 speckle）；speckle 按多数重叠分配到所属核。

### 其他

| 参数 | 默认 | 说明 |
|---|---|---|
| `min_field_cells` | 10 | 生殖细胞少于该数的场不通过 QC，不计入组间统计 |
| `pair_rg_min` / `pair_br_min` | 0.20 / 0.15 | 三色配对的相关系数下限 |
| `sat_thresh` | 254 | 饱和像素判定阈值 |

## N/C 比值的定义

N/C = (核内 polyA 均值 − 场背景) / (胞质 polyA 均值 − 场背景)

- 核内取样 = 核掩膜内缩 `nuc_erode` 像素（当前 0 = 全核）；胞质取样 = DDX4⁺ 距核 ≤`max_ring_px` 的环（或 ilastik 类别2 胞质 ∩ 核外带），经 Voronoi 分配到最近的核
- 场背景 = 非核、非 DDX4⁺ 区域的 polyA 中位数
- **均值口径**（浓度比），非总强度

## 输出文件（写入输出目录）

| 文件 | 内容 |
|---|---|
| `per_cell_quantification.csv` | 每细胞（核）一行，列含义见下表 |

#### per_cell_quantification.csv 列说明

**标识与分割**

| 列 | 含义 |
|---|---|
| `set` | 组名_场名（如 `URH-WT-merge_DDX4-WT-4-DAPI`），该核来自哪个组哪一场 |
| `id` | 该场内本核的编号 |
| `seg_method` | 本场分割方式：`ilastik_prob`（概率图）或 `dapi_adaptive`（缺 h5 时回退） |
| `group` | 组名（按组文件夹名前缀匹配） |
| `x`, `y` | 核质心像素坐标 |

**形态学（像素级）**

| 列 | 含义 |
|---|---|
| `area_px` | 核面积 = 核掩膜内像素数（RS 判据之一） |
| `ar` | 长短轴比，1=正圆（RS 要求 < `rs_ar_max`） |
| `circ` | 圆度 = 周长²/(4π·面积)，1=正圆（RS 要求 < `rs_circ_max`） |

**强度测量（0–1 标度）**

| 列 | 含义 |
|---|---|
| `nuc_polya_raw` | 核内 polyA 平均强度（全核采样，未校正） |
| `ann_polya_raw` | 胞质取样区内 polyA 平均强度（未校正） |
| `nuc_npix` / `ann_npix` | 核 / 胞质环的采样像素数 |
| `nuc_polya_bc` | 核内 polyA 背景校正值 = `nuc_polya_raw − bg_polya`（N/C 分子） |
| `ann_polya_bc` | 胞质 polyA 背景校正值 = `ann_polya_raw − bg_polya`（N/C 分母） |
| `nc_ratio` | **N/C = `nuc_polya_bc / ann_polya_bc`**，主要读数（浓度比） |
| `ann_ddx4_frac` | 胞质环内 DDX4⁺ 像素占比（`ddx4` 模式生殖细胞门控用） |
| `ann_ddx4_mean` | 胞质环内 DDX4 平均强度（`ilastik_cyto` 模式门控用，≥ `min_ddx4_mean`） |
| `ring_width_px` | 胞质环平均宽度 = `ann_npix / 核周长`（QC：过窄 = 胞质取样不足） |

**分类标签**

| 列 | 含义 |
|---|---|
| `germ_cell` | TRUE = 生殖细胞（DDX4 门控通过 + 背景校正后信号有效 + 胞质像素 ≥ `min_ann_pix`） |
| `is_round_spermatid` | TRUE = 圆精子细胞（germ_cell 且面积 < `rs_area_max` 且 `ar` < `rs_ar_max` 且 `circ` < `rs_circ_max`）；**组间统计只用这些行** |

**核内 speckle（阈值 = 该核 polyA 均值 × `speckle_nuc_fold`）**

| 列 | 含义 |
|---|---|
| `n_speckles` | 该核内斑点个数（整核均匀亮的核按定义 = 0） |
| `nuc_speckle_frac` | 核内超过斑点阈值的像素占比（"斑点化程度"） |
| `speckle_npix` | 该核所有斑点的像素总数（斑点面积负载） |
| `speckle_intensity` | 该核斑点像素的 polyA 平均强度（斑点亮度，原始标度；无斑点 = NA） |
| `summary_per_image.csv` | 每场一行：seg_method、核数、排除的长形核数、生殖细胞数、RS 数、RS 中位/均值 N/C、全场生殖细胞中位 N/C、n_bright_nuc（场级亮核物体数）、frac_nuc_polya_bright（核内亮区占比）、**n_speckles_total（场内核 speckle 总数）**、**n_cells_with_speckles / frac_cells_with_speckles（有 speckle 的细胞数/占比）**、**median_speckles_per_positive_cell（阳性细胞的中位 speckle 数）**、**mean_speckle_area_px（speckle 平均面积）**、**median_speckle_intensity（本场有 speckle 细胞的 speckle 强度中位数）**、bg_polya、group、field_pass |
| `pairing_log.csv` | 三色配对结果与相关系数 |
| `group_comparison_stats.csv` | 场级 RS 中位 N/C 的 WT vs cKO：中位数、均值差、95%CI、t 检验与 Wilcoxon p 值 |
| `group_comparison_bright_nuclei.csv` | 场级亮核物体数的组间比较（局灶表型主读数） |
| `group_comparison_stats_ilastik_only.csv` | 仅用 ilastik 分割场的组间统计（每组 ≥3 场时输出） |
| `NC_WT_vs_cKO_boxplot.svg` | 场级 RS 中位 N/C 箱线图（可编辑 SVG） |
| `NC_RS_distribution_overlay.svg` | 两组 RS N/C 分布叠加密度图 |
| `<组>_<场名>_QC_overlay.png` | 每场 QC 图：青=保留核、红=长形核排除、橙=RS、绿=核内 speckle、黄晕=胞质环 |
| `<组>_<场名>_NC_map.png` | 每场 N/C 着色图：蓝=核质比>1，品红=≤1 |

## 统计单位与注意事项

- 统计单位是**场（视野）**而非细胞；每组仅 1 只小鼠时，场为技术重复，结论属**探索性**
- 两种分割/胞质方法混用时日志会打 WARNING；正式比较请让所有场用同一方法（补齐概率文件或设 `require_prob = TRUE`）
- 组间比较前请确认两组参数完全一致
