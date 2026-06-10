# ============================================================
# Version10.2 root-level extension runner
# ============================================================
# Standalone run:
#   source("run_all_version10_extensions.R")
#
# main.R does not depend on this file; it directly calls the four root-level
# extension scripts. This file is only a convenience launcher.
# ============================================================

source("functions/config.R")

if (isTRUE(config$run_version10_extensions)) {
  if (isTRUE(config$run_extension_01_hidden_strength_vbm_msm)) {
    source("extension_01_hidden_strength_vbm_msm.R")
  }
  if (isTRUE(config$run_extension_02_vbm_ps_misspecification)) {
    source("extension_02_vbm_ps_misspecification.R")
  }
  if (isTRUE(config$run_extension_03_good_overlap_instability)) {
    source("extension_03_good_overlap_instability.R")
  }
  if (isTRUE(config$run_extension_04_rgm_vbm_msm)) {
    source("extension_04_rgm_vbm_msm.R")
  }
}
