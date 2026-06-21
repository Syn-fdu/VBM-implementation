# RGM / MSM / VBM comparison

|---|---|---|---|---|---|
| VBM | R2 | Variance of weights | Hidden perturbation is measured through extra weight variance. | Interpretable covariate benchmark scale. | Can be unstable when observed weight variance is near zero. |
| MSM | Gamma | Bounded odds-ratio multiplier | Ideal control weights differ from observed weights by a bounded multiplier. | Classic sharp marginal sensitivity model. | Gamma calibration can be conservative and ratio-driven. |
| RGM-sharp | T | Total variation/L1 mass shift | The normalized ideal control-weight distribution lies within TV distance T. | Finite-sample transport interpretation: move probability mass from low to high outcomes. | T calibration requires distributional benchmark interpretation. |
| RGM-conservative | T | Total variation/L1 range bound | The same TV radius is combined with an outcome range bound. | Simple closed-form conservative bound. | May be loose when the sample outcome range is wide. |
