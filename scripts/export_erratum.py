"""Export the current erratum tables and synchronize its self-contained plot inputs."""
from pathlib import Path
import csv
from decimal import Decimal, localcontext
import shutil

root = Path(__file__).resolve().parents[1]
out = root / "results/paper_method_corrected"
support = root / "docs/erratum_besii/support"
support.mkdir(parents=True, exist_ok=True)

def read(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))

def write(name, rows, fields):
    with (out / name).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        w.writerows({k: row[k] for k in fields} for row in rows)

rows = read(out / "results.csv")
observed = read(root / "data/total_input.csv")[1:]
assert len(rows) == len(observed) == 209
for row, obs in zip(rows, observed):
    assert abs(float(row["energy_GeV"]) - float(obs["energy_GeV"])) < 1e-12
    row["sigma_observed_nb"] = obs["sigma_nb"]
    # Decimal arithmetic keeps this derived column identical across Python versions/platforms.
    with localcontext() as context:
        context.prec = 50
        variance = sum(Decimal(v)**2 for k, v in obs.items() if k not in ("energy_GeV", "sigma_nb"))
        row["sigma_observed_error_nb"] = format(variance.sqrt(), ".16g")
write("isr.csv", rows, ["energy_GeV", "ISR_correction", "ISR_total_error"])
write("vp.csv", rows, ["energy_GeV", "VP_factor", "VP_total_error"])
write("cross_sections.csv", rows, ["energy_GeV", "sigma_observed_nb", "sigma_observed_error_nb", "sigma_dressed_nb", "sigma_dressed_error_nb", "sigma_Born_nb", "sigma_Born_error_nb"])
for name in ["results.csv", "vp_dense.csv", "error_components.csv", "validation.toml", "resolution.toml"]:
    shutil.copyfile(out / name, support / ("paper_method_" + name))
for name in ["isr", "vp", "sigma_dressed", "sigma_born"]:
    shutil.copyfile(out / ("total_" + name + "_covariance.csv"), support / ("paper_method_" + name + "_covariance.csv"))
shutil.copyfile(root / "data/total_input.csv", support / "observed_input.csv")
print("Exported 209 rows and synchronized erratum plot inputs.")
