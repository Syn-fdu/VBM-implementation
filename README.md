# Version10.1 folder layout without `extensions/`

这版按照新的要求重构：

1. zip 解压后得到一个外层文件夹：`version10.1/`
2. 项目内部没有 `extensions/` 文件夹。
3. `main.R` 直接调用根目录下四个可以单独运行的 extension 脚本：
   - `extension_01_hidden_strength_vbm_msm.R`
   - `extension_02_vbm_ps_misspecification.R`
   - `extension_03_good_overlap_instability.R`
   - `extension_04_rgm_vbm_msm.R`
4. 四个 enhancedplot 脚本统一放在 `scripts/` 下，并由对应 R 脚本自动启动。
5. 输出统一写入：`output/version10_extension_results/`

## 推荐运行方式

```r
setwd("C:/Users/11150/Desktop/生物统计/version10.1")
source("main.R")
```

执行 `main.R` 时，默认会顺带运行四个 extension。开关仍在：

```r
functions/config.R
```

其中默认：

```r
run_version10_extensions = TRUE
run_extension_01_hidden_strength_vbm_msm = TRUE
run_extension_02_vbm_ps_misspecification = TRUE
run_extension_03_good_overlap_instability = TRUE
run_extension_04_rgm_vbm_msm = TRUE
```

## 单独运行某个 extension

```r
source("extension_01_hidden_strength_vbm_msm.R")
source("extension_02_vbm_ps_misspecification.R")
source("extension_03_good_overlap_instability.R")
source("extension_04_rgm_vbm_msm.R")
```

## 主要目录结构

```text
version10.1/
  main.R
  extension_01_hidden_strength_vbm_msm.R
  extension_02_vbm_ps_misspecification.R
  extension_03_good_overlap_instability.R
  extension_04_rgm_vbm_msm.R
  functions/
  rgm_model/
  scripts/
  data/
  output/
```
