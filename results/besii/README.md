# Published BES input: reader data

> Historical input-only BES replacement, superseded by the current erratum results in `results/paper_method_corrected/`. See the repository README for current figures and data.

These files contain the results of the `besii-paper-consistency` branch. They use published BES inclusive R measurements in the charm region and retain the historical numerical algorithm. The reference calculation is `paper-reproduction` at commit `5cabb61`.

## Ready-to-use tables

Every result table contains the same 209 bottomonium energies, in ascending order. CSV files have one header row, decimal numeric values, and no index column.

| File | Contents |
| --- | --- |
| [results.csv](results.csv) | Complete current results, including factors, cross sections and errors |
| [isr.csv](isr.csv) | ISR correction and total propagated error |
| [vp.csv](vp.csv) | VP factor and total propagated error |
| [cross_sections.csv](cross_sections.csv) | Dressed/Born cross sections and bottom-quark R values, with errors |
| [comparison.csv](comparison.csv) | Old/new ensemble means, errors and percentage changes |
| [paired_vp_shift.csv](paired_vp_shift.csv) | Paired VP mean shift and its Monte Carlo standard error |
| [Published BES input](../../data/besii/input.csv) | 77 input energies, inclusive R, statistical/systematic errors and publication year |

The three shorter tables are exact column selections of `results.csv`, with no rounding or recomputation. The full table is the canonical source.

## Column definitions

| Columns | Meaning and units |
| --- | --- |
| `energy_GeV` | Center-of-mass energy in GeV |
| `ISR_correction` | Multiplicative correction taking observed to dressed cross section |
| `VP_factor` | Vacuum-polarization factor; divide dressed cross section by it to obtain Born |
| `sigma_dressed_nb`, `sigma_Born_nb` | Cross sections in nb |
| `R_dressed`, `R_Born` | Bottom-quark cross-section ratios, dimensionless |
| `*_error`, `*_error_nb` | Propagated sample standard deviations, in the corresponding quantity's units |

The historical normalization is

$$\sigma_{\mu\mu}(s)=\frac{86.85}{s}\,\mathrm{nb},\qquad s=E_{\mathrm{cm}}^2\text{ in GeV}^2.$$

Cross sections were computed for each toy using the original unsmoothed experimental sample:

$$\sigma_t^{\mathrm{dressed}}=\sigma_t^{\mathrm{observed}}C_{\mathrm{ISR},t},\qquad
\sigma_t^{\mathrm{Born}}=\sigma_t^{\mathrm{dressed}}/F_{\mathrm{VP},t}.$$

The cross-section columns are means of these paired products/ratios. Multiplying the separately reported means will not in general reproduce them exactly. ISR means use the historical statistical ensemble; all other means and all reported standard deviations use total ensembles. Each ensemble has 10,000 samples. Reported errors are not Monte Carlo standard errors and do not include unquantified numerical/model systematics.

## Covariance matrices

[ISR](isr_covariance.csv) · [VP](vp_covariance.csv) · [Dressed cross section](sigma_dressed_covariance.csv) · [Born cross section](sigma_born_covariance.csv)

Each matrix is 209 by 209 after skipping the header row. Header entries are column energies; rows follow the same energy order. There is **no leading row-label column**. Entries have squared units (dimensionless for correction factors, square nb for cross sections). Diagonal entries equal squared total errors in `results.csv`. These are separate within-observable covariance matrices, not a joint covariance across all observables.

## Figures and validation

- [Overview with propagated error bands](figures/overview.png) ([PDF](figures/overview.pdf)).
- [Relative changes](figures/relative_changes.png) ([PDF](figures/relative_changes.pdf)).
- [Charm input comparison](figures/input_comparison.png) ([PDF](figures/input_comparison.pdf)).
- [Numerical validation](validation.toml), [source checks](source_checks.json), [output checks](output_checks.json).
- [Scope and known limitations](../../docs/besii_results.md).

`fixed_input_iterations.csv` is a diagnostic with one fixed background draw, not a nominal central-value table. Older directories outside `results/besii` retain historical results. This branch is an input-only revision; it does not incorporate the numerical fixes on `numerical-updates`.

After rerunning the Julia calculation, refresh the shorter tables with `python3 scripts/export_besii_tables.py` from the repository root.
