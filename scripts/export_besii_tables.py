"""Export exact column selections from the canonical BES result table."""
from pathlib import Path
import csv

out = Path(__file__).resolve().parents[1] / "results/besii"
with (out / "results.csv").open() as stream:
    rows = list(csv.DictReader(stream))
assert len(rows) == 209
columns = {
    "isr.csv": ["energy_GeV", "ISR_correction", "ISR_total_error"],
    "vp.csv": ["energy_GeV", "VP_factor", "VP_total_error"],
    "cross_sections.csv": ["energy_GeV", "sigma_dressed_nb", "sigma_dressed_error_nb", "sigma_Born_nb", "sigma_Born_error_nb", "R_dressed", "R_dressed_error", "R_Born", "R_Born_error"],
}
for name, fields in columns.items():
    with (out / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows({key: row[key] for key in fields} for row in rows)
    with (out / name).open() as stream:
        exported = list(csv.DictReader(stream))
    assert len(exported) == len(rows)
    assert all(row[key] == original[key] for row, original in zip(exported, rows) for key in fields)
print("Exported and verified three tables (209 rows each).")
