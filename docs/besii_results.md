# BES input replacement: results and limitations

> Historical input-only BES replacement, superseded by the current erratum results in `results/paper_method_corrected/`. See the repository README for current figures and data.

All $209$ bottomonium experimental energies were recomputed with the historical algorithm and the published BES charm-region inputs. The calculation used two ISR ensembles of $10000$ samples and paired old/new VP ensembles of $10000$ samples each. The comparison baseline is `paper-reproduction` commit `5cabb61`.

| Quantity | Change relative to the historical input |
| --- | --- |
| ISR correction and its uncertainty | Unchanged; historical replay checked |
| Dressed cross section and its uncertainty | Unchanged; identical raw samples and ISR factors |
| VP mean | Decreases by $0.001829\%$–$0.001893\%$ |
| VP standard deviation | Increases by $5.318\%$–$6.688\%$ |
| Born cross-section mean | Increases by $0.001772\%$–$0.001922\%$ |
| Born standard deviation | Relative change between $-0.004734\%$ and $+0.011688\%$ |

These percentages compare means with means and standard deviations with standard deviations. The largest VP mean change occurs at $10.5628\,\mathrm{GeV}$; the largest Born mean change occurs at $10.6308\,\mathrm{GeV}$. The maximum Monte Carlo standard error of the *paired VP shift* is $1.81\times10^{-6}$ in absolute VP units. It is distinct from the propagated experimental error, approximately $0.00075$–$0.00083$ after replacement.

## Definitions

The final cross sections are formed using the original unsmoothed experimental samples:

$$\sigma_{t}^{\rm dressed}=\sigma_{t}^{\rm observed}C_{{\rm ISR},t},\qquad
\sigma_{t}^{\rm Born}=\frac{\sigma_{t}^{\rm dressed}}{F_{{\rm VP},t}}.$$

The ensemble mean of the ratio is reported; it is not replaced by the ratio of means. Thus the Born shift need not be the exact negative of the VP-mean shift. Both correction factors and the experimental inputs are paired within each sample. Covariance CSV files have matching row/column energy ordering; the first row contains column energies, and all $209$ rows follow `results.csv`.

The historical normalization $\sigma_{\mu\mu}=86.85/s$ in nb is retained for reported $R_b$. ISR means follow the historical statistical ensemble convention, while ISR errors come from the total ensemble. VP, dressed and Born means and errors use all historical total samples. This includes a fixed first bottom sample, whereas the first low-energy and charm-background samples are random. No first-sample output is mislabeled as a nominal central prediction.

## Validation

- $35/35$ Julia tests passed, including published-input checks, the splice at $5.05\,\mathrm{GeV}$ and direct-sum agreement.
- All $76$ BES 2002 points match the primary paper's central values and statistical errors exactly. PDG-derived systematic errors differ from the rounded paper values by at most $0.00192$ in R. The BES 2000 endpoint is $R=3.47\pm0.32\pm0.29$ after rounding.
- Historical VP samples reproduce the saved reference within $4.89\times10^{-15}$ absolute.
- The new VP matrix calculation agrees with the literal trapezoid-sum implementation within $6.67\times10^{-16}$ for samples $1,2,101,10000$.
- VP iteration $3$ versus $6$ changes the fixed-input diagnostic by at most $4.81\times10^{-8}$ absolute.
- ISR iteration sensitivity is still present: $30$ versus $60$ states changes the correction by up to $0.001671$ absolute. Keeping that setting is part of this branch's historical scope, not a convergence claim.

This validates the **input-replacement comparison in the historical implementation**, not the correctness of all historical formulas or a new absolute-accuracy target. In particular, the hadronic imaginary-term coefficient and pion conversion constant have not been corrected here. Numerical/model uncertainties, smoothing dependence and the true experimental covariance remain outside the quoted errors. The fully correlated statistical and systematic components are inherited from the original calculation, not asserted as the BES experimental covariance.

## Files

- `results.csv`: current correction factors, cross sections and errors.
- `comparison.csv`: historical/revised summaries and relative changes.
- `*_covariance.csv`: paired total-sample covariance matrices.
- `paired_vp_shift.csv`: shift and its Monte Carlo standard error.
- `fixed_input_iterations.csv`: iteration diagnostic, explicitly not a nominal prediction.
- `validation.toml`, `source_checks.json`: numerical and primary-source checks.
- `figures/`: input comparison, overview with error bands, and relative changes.

The run is self-contained. The thesis, desktop archive, external disk and MATLAB are not runtime dependencies. The baseline inputs under `reference/legacy_charm/` are used only to produce and verify the comparison.
