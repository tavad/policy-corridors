# Policy corridors as a diagnostic: the R code and data

The R code and data reproducing every figure and table in:

> Tavadyan, A. and Tavadyan, A. (2026). *Policy Corridors as a Diagnostic Tool:
> Cross-Country Evidence from Six Decades of Adherence.*

Archived on Zenodo (DOI: ⟨ZENODO_DOI⟩).

![Adherence by indicator and development group, 2008–2023.](plots/F1_adherence.png)

## Contents

| Path | What it is |
|---|---|
| `R/panel.R` | The canonical 217-economy sample drawn from the detailed panel (drops the eight World Bank aggregates). Sourced by the scripts below; not run on its own |
| `R/historical_test.R` | Long-panel adherence, 1960–2023: by era, growth take-offs, and the event-time series (**Figures 2 and 4, Table 3, §4.2 and §4.5**). Writes `data/historical_20260816/` |
| `R/necessity_test.R` | The corridor test (below / inside / above each bound) and the necessity check (**Figure 3, Table 5, §4.4 and §4.6**). Writes `data/necessity_20260816/` |
| `R/joint_compliance.R` | The joint-score forward outcomes, the pre-trends check and Benjamini–Hochberg multiplicity (**§4.6 and Appendix A, joint-score group**). Writes `data/joint_results_table.csv` |
| `R/make_tables.R` | Renders **Tables 2, 3, 4, 5 and B1** from the panels and the intermediates above. Writes `outputs/` |
| `R/make_figures.R` | Renders **Figures 1–4** from `outputs/T3_adherence.csv` and the intermediates. Writes `plots/` |
| `data/magic_squares_data_clean.csv` | The detailed panel, 2008–2023, 217 economies (five indicators). Frozen January 2026 |
| `data/wdi_long_1960_2023.csv` | The long panel, 1960–2023, three back-reaching indicators (World Bank WDI, pulled 4/16 August 2026) |
| `plots/` | Figures 1–4 as published (PNG) |
| `outputs/` | Where `make_tables.R` writes Tables 2–5 and B1 (pandoc pipe format); empty until you run it |

The generated table files are named by build order, not by manuscript number:
`T3_adherence` → **Table 2**, `T4_adherence_by_era` → **Table 3**, `T6_shock_years` → **Table 4**,
`T5_corridor_test` → **Table 5**, `T2_data_coverage` → **Table B1**. Table 1 is the framework's
normative bounds and is stated in the manuscript, not computed.

## Reproducing the figures and tables

Base R only. Install the `here` package (used to anchor paths to this directory via the
empty `.here` marker) and, for the figures, `ragg` (falls back to base `png()` if absent).
No API keys, no other dependencies. Run from the repository root, in this order:

```sh
Rscript R/historical_test.R      # -> data/historical_20260816/
Rscript R/necessity_test.R       # -> data/necessity_20260816/
Rscript R/joint_compliance.R     # -> data/joint_results_table.csv
Rscript R/make_tables.R          # -> outputs/  (Tables 2, 3, 4, 5, B1)
Rscript R/make_figures.R         # -> plots/    (Figures 1-4)
```

`historical_test.R` and `necessity_test.R` each write a folder of intermediate CSVs
(`data/historical_20260816/`, `data/necessity_20260816/`), created on run and git-ignored;
`joint_compliance.R` writes `data/joint_results_table.csv`. `make_tables.R` and `make_figures.R`
read the intermediates. Every plotted and tabulated value is read from a file at run
time, so a figure or table cannot drift from the numbers behind it.

## The data

Two panels, both from the World Bank World Development Indicators.

The **detailed panel** (`magic_squares_data_clean.csv`, 2008–2023) is frozen as extracted in
January 2026. A fresh WDI pull does **not** reproduce it: the World Bank revises the series and
has renamed identifiers, so a later release moves observations. The manuscript reports that
adherence is stable under this revision (the largest single cell moves 2.2 points, the two group
means 0.1 and 0.9), but the frozen file is what the published numbers were built from and is the
file to use.

The **long panel** (`wdi_long_1960_2023.csv`) was pulled from WDI on 16 August 2026 for the three
indicators that reach back to 1960: inflation (`FP.CPI.TOTL.ZG`), broad money as a share of GDP
(`FM.LBL.BMNY.GD.ZS`) and the trade balance (`NE.RSB.GNFS.ZS`), together with GDP growth
(`NY.GDP.MKTP.KD.ZG`) and GDP per capita in constant 2015 US dollars (`NY.GDP.PCAP.KD`) for the
outcome and the development classification. Unemployment and the real lending rate do not reach
back far enough, which is why the long panel scores three corridors rather than five. A re-pull
will differ for the same reason as the detailed panel; the frozen file reproduces the paper.
