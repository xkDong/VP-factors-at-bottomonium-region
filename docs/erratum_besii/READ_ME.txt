Current erratum manuscript

main.tex and main.pdf are the current manuscript. figures/ contains the three displayed figures; support/ contains the current plotting data. The Julia results are in results/paper_method_corrected at the repository root.

From the repository root run python3 scripts/export_erratum.py, then python3 docs/erratum_besii/revision/generate_revision.py. To compile, run pdflatex twice on main.tex from this directory. The TeX build uses only local assets.

Historical local drafts and audit notes are not part of this distribution.
