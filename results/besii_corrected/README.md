# Corrected calculation with published BES inputs

These results combine the `numerical-updates` implementation with the published inclusive BES input. This is an **intermediate, superseded result set**. The current erratum uses [paper_method_corrected](../paper_method_corrected/README.md), with robust LOWESS, Gaussian MLE summaries and separately propagated uncertainty components. These files are retained for independent numerical checks.

- `results.csv`: 209 measured energies, ISR and VP factors, dressed and Born cross sections in nb, and propagated input standard deviations.
- `vp_dense.csv`: VP evaluated on the refined solution and exported at five subdivisions per measured interval, including the lower region for the 10.4–11.2 GeV plot.
- `*_covariance.csv`: covariance matrices, in the order of `results.csv`.
- `comparison.csv`: paired comparison with the historical algorithm and original charm input, using the same 9999 random samples and mean definition.
- `fixed_inputs.csv`: a separate fixed-input calculation for diagnostics. Its low-energy background uses the stored ensemble mean, not an asserted experimental nominal input.
- `validation.toml`, `independent_methods.toml`, and `tests.log`: convergence, independent integration, and regression evidence.

Central values are the means of 9999 paired random samples, excluding the first fixed bottom sample. The background first sample is not assumed to be nominal. Final cross sections apply the correction factors to the unsmoothed observed inputs. The numerical update retains signed linear estimates. Thus the difference from the historical positive multiplicative estimator is not entirely a numerical integration error.

The uncertainty model is inherited from the stored inputs; smoothing/model uncertainties and missing parameter uncertainties are not included. No complete statistical/systematic split is claimed. Resonance parameters have not been refitted.

Run from the repository root:

```sh
julia --project=. scripts/run_besii_corrected.jl
julia --project=. scripts/check_besii_corrected_methods.jl
julia --project=. test/runtests.jl
```

The pQCD segment starts at 5.05 GeV, after the last BES input at 5.00 GeV. Production VP uses refinement 20, background refinement 2, a 200 GeV cutoff and 0.125 GeV tail spacing, with relative iteration tolerance $10^{-11}$. The combined fixed-input refinement sensitivity is approximately $5.20\times10^{-7}$ in the VP factor; this is an empirical resolution estimate, not a rigorous error bound.
