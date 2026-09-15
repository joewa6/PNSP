# Manuscript figures and tables

This directory contains the processed numerical data and plotting scripts for
the manuscript.

## Reproduce the outputs

From the repository root:

```bash
python paper/scripts/make_main_figures.py
python paper/scripts/make_supplementary_figures.py
python paper/scripts/make_tables.py
python paper/scripts/verify_release.py
```

Outputs are written to `paper/generated/`.

| Script | Outputs |
|---|---|
| `make_main_figures.py` | Figures 3–8 |
| `make_supplementary_figures.py` | Figures S2–S13 |
| `make_tables.py` | Tables S1–S2 |
| `verify_release.py` | Data and output validation |

## Data organization

`processed_data/` contains the numerical inputs used by the plotting scripts.
`manifest.csv` records their SHA-256 checksums and dimensions.

The principal screen contains 203 contaminant systems, three replicas per
system, and ten contaminant copies per replica. Production trajectories were
15 ns long and saved every 10 ps.

For each copy, occupancy is the fraction of frames in which the minimum
contaminant–PET heavy-atom distance is less than 4.0 Å. System occupancy is the
arithmetic mean of the 30 copy occupancies. A contact event is a maximal
uninterrupted run of bound frames. Mean system residence time is calculated
from the event-positive copy-level mean event durations pooled across the three
replicas.
