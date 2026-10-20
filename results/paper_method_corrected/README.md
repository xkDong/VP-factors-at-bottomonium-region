# Corrected results with the requested article methods

The main erratum uses this result set. It retains published BES data and the corrected VP constants and integration kernel, uses robust local-linear LOWESS (span 0.045, start 10.6659 GeV, three bisquare residual-reweighting passes), and stops the ISR/VP iteration at a maximum relative factor change below $10^{-4}$ (0.01%). ISR uses a signed Jacobi iteration of the same triangular convolution system. No input or estimate is clipped or discarded.

The seed is 20261020 and all 10000 samples are random Gaussian draws. Central values and widths are unbinned Gaussian maximum-likelihood estimates. The covariance normalization uses the same sample count as the fitted variances. This fit does not establish normality: see the per-point skewness and excess kurtosis files. Near-constant factors (width below $10^{-10}$) are excluded from aggregate shape diagnostics.

## Uncertainty components

- **stat**: bottomonium and BES statistical errors only; all systematic nuisances and the low-energy background fixed.
- **sys**: bottomonium and BES systematic errors only; all statistical nuisances and the low-energy background fixed.
- **low_background**: the unresolved combined low-energy error alone; it is not relabelled statistical or systematic.
- **total**: all components varied simultaneously. Its covariance is calculated directly, not forced to equal the sum of component covariances.

Bottomonium point-to-point errors are independent; common systematic nuisances are shared within each experiment. In accord with the article's declared background approximation, each BES component is fully correlated across energies, using independent statistical and systematic nuisances. The low-energy combined error uses one common nuisance. These are adopted correlation assumptions, not recovered experimental covariance matrices. Additional smoothing/model uncertainties are not included.

`results.csv` gives factors and cross sections with total fitted widths; `error_components.csv` provides all four uncertainty components. The `*_covariance.csv` and `*_gaussian_fit.csv` files provide covariance and fit diagnostics. The low-energy input is copied from the archived `inputR.mat` export (663 points); its fourth column is zero, and the third-column uncertainty is deliberately left unclassified.

Run from the repository root:

```sh
julia --project=. scripts/run_paper_method_corrected.jl
julia --project=. scripts/check_paper_resolution.jl
julia --project=. test/runtests.jl
```

Numerical validation is in `validation.toml` and `resolution.toml`. The unchanged corrected convolution/PV functionals also have independent integration checks in `results/besii_corrected/independent_methods.toml`. Old result sets remain available for provenance and are not the current erratum figures. No resonance fit was rerun.

## Reader tables

- [ISR](isr.csv) and [VP](vp.csv): factors and total fitted widths at measured energies.
- [Cross sections](cross_sections.csv): observed, dressed and Born cross sections and total errors, in nb.
- [Dense VP curve](vp_dense.csv): the plotted central curve and uncertainty band.
- [All uncertainty components](error_components.csv): conditional statistical, systematic, low-background and simultaneous total widths.

Generate the reader tables and synchronize the manuscript figure inputs with `python3 scripts/export_erratum.py`.
