# VAR-MERTON — Replication package

Replication code and data for the paper:

> **Generalized Impulse Responses of Portfolio Default Probabilities:  
> A Modular Framework with an Application to Geopolitical Risk** — G. Flament, C. Hurlin, Q. Lajaunie, Y. Pull.

Working paper: [SSRN 7106918](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7106918)
(doi:[10.2139/ssrn.7106918](https://doi.org/10.2139/ssrn.7106918)).


The framework propagates an identified macro-financial innovation through a
Bayesian VAR, maps it into a latent systematic credit factor via a satellite
equation, and converts it into closed-form generalized impulse responses of the
mean, the quantiles (PD-at-Risk) and the expected shortfall of portfolio
default probabilities (Merton-Vasicek map). The application studies U.S.
geopolitical-risk (GPR) shocks.

## Requirements

- R (version pinned in `renv.lock`); package versions are restored by `renv`.
- No internet access is required at run time: all input data are frozen in
  `data/raw/` (see *Data* below).
- A multi-core machine is recommended; posterior computations use up to 60%
  of the available logical cores (option `robustness.core_fraction`). The
  reported results were produced on an Apple M4 Max (128 GB RAM).

## How to run

Start R from the repository root:

```r
install.packages("renv")
renv::restore()

source("run_all.R")
```

`run_all.R` is the single entry point; there is no other orchestrator. It
executes the full pipeline (data construction, estimation, all exercises,
paper deliverables) and prints a per-step recap with timings. Indicative
runtime: several hours on a recent laptop; the heaviest steps are the
information-set robustness grids and the per-window model averaging.

After a run, a lightweight sanity check verifies that the expected
deliverables exist with sane shapes — including that every figure name the
LaTeX source includes is present in `output/paper_replication/figures/`. It
re-estimates nothing and exits non-zero on failure:

```r
Rscript tests/smoke_test.R
```

A second, self-contained test checks the closed-form expected-shortfall
response against numerical integration:

```r
Rscript tests/test_pd_es_closedform.R
```

The last pipeline step stops with an error if any paper deliverable is missing,
rather than leaving a stale copy in place. To inspect a deliberately partial
run, set `options(paper.strict_deliverables = FALSE)` before sourcing
`code/shared/make_paper_outputs.R`.

## Paper deliverables

The final pipeline step (`code/shared/make_paper_outputs.R`) assembles every
figure and table included in the paper, under the exact file names referenced
by the LaTeX source, into:

```text
output/paper_replication/
  figures/main/            main-text figures
  figures/robustness/      appendix figures (one folder per VAR information set)
  tables/                  LaTeX tables
  paper_numbers.csv        every scalar quoted in the text
```

To compile the paper, copy `output/paper_replication/figures/` next to the
`.tex` source.

## Pipeline map

`run_all.R` runs the full pipeline:

| Step (order of `run_all.R`) | Produces | Paper objects |
|---|---|---|
| `code/01_process_gpr.R`, `code/02_build_dataset.R` | quarterly dataset (`data/processed/`) | sample of Section 3.1 |
| `code/shared/var_information_set_diagnostics.R` | BVAR stability diagnostics | stability table (appendix) |
| `code/dralacbn/02_*.R` | baseline specification + information-set robustness | dynamic responses, satellite tables, appendix figures |
| `code/dralacbn/03_*.R`, `code/dralacbn/04_*.R` | direct-channel / orthogonality test (relaxes `eta`$\perp$`u`) | control-function table, PD with/without direct channel, `lambda` posterior |
| `code/dralacbn/07_*.R`, `code/dralacbn/09_*.R` | per-window BMA short-default-sample responses (**heavy step**) | Section 3.4, shorter-default-histories figures (2005, 2015 windows) |
| `code/dralacbn/08_dralacbn_perfect_foresight.R` | perfect-foresight bias | Section 3.5 |
| `code/dralacbn/10_dralacbn_state_dependence.R` | historical episodes, state dependence | Section 3.3 |
| `code/review/04_simulation_benchmark.R` | closed form vs forward Monte Carlo, at a fixed parameter vector | Section 2.5 benchmark figure |
| `code/shared/variable_definitions.R` | variable glossary | documentation |
| `code/shared/collate_paper_outputs.R` | collated outputs (`output/paper/`) | working copies |
| `code/shared/make_paper_outputs.R` | `output/paper_replication/` | all paper figures/tables, paper-named |

```r
source("run_all.R")
```

## Sample conventions and reproducibility

The estimation sample ends at `DATA_END_DATE` (set in `config.R`, currently
2024-12-31). The truncation is applied centrally in `code/02_build_dataset.R` and
wherever default-rate series are loaded, so every block shares the same
cutoff. **If you change `DATA_END_DATE`, delete `output/shared/var_kernels/`**
so that the cached BVAR posterior kernels are re-estimated on the new sample.

All randomness is seeded (`SEED*` constants in `config.R`) with the sampling
kind pinned in `code/00_setup.R`, so a full re-run reproduces the reported
numbers exactly under the pinned package versions.

## Data

All inputs are frozen files under `data/raw/`, inventoried with their sources
in `data/raw/MANIFEST.md`. No live download takes place during a replication
run (`ALLOW_DOWNLOADS = FALSE` in `config.R`). Main series:

- **Geopolitical Risk Index** (daily; monthly companion file), Caldara and
  Iacoviello (2022), <https://www.matteoiacoviello.com/gpr.htm>;
- **FRED series** (delinquency rate `DRALACBN`, national accounts, prices,
  employment, interest rates), Federal Reserve Bank of St. Louis;
- **EPU index** (Baker, Bloom and Davis), **VIX** (FRED: `VIXCLS`), and an
  S&P 500 snapshot, used only in the robustness information sets.

The raw files are included for replication convenience; they remain the
property of their providers and are subject to the providers' terms of use.
In particular, the frozen **S&P 500 snapshot** (`data/raw/sp500_GSPC_snapshot.csv`,
used by the `uncertainty` robustness information set) is redistributed here
only to make a replication run self-contained; it remains subject to S&P Dow
Jones Indices' terms. It can otherwise be regenerated at run time by setting
`ALLOW_DOWNLOADS = TRUE` in `config.R` (the code re-fetches `^GSPC`).

## Repository layout

```text
config.R                  all paths and parameters (single source of truth)
run_all.R                 entry point: full pipeline
dependencies.R            package list (informational; use renv.lock)
R/functions.R             legacy helpers kept for reference (not on the
                          pipeline path)
code/
  00_setup.R              shared initialization (sourced by every script)
  01_*, 02_*              data construction
  README.md               notes on the code layout
  shared/                 estimation engine, GIRF machinery, writers,
                          make_paper_outputs.R (paper deliverables)
  dralacbn/               baseline application (delinquency proxy):
                          02 (application + robustness), 03+04 (direct-channel
                          orthogonality test), 07+09 (per-window BMA
                          short-default-sample figures), 08 (perfect foresight),
                          10 (state dependence)
  review/                 exercises written in response to referee points; 04
                          (closed form vs forward simulation) is on the paper
                          pipeline and produces the Section 2.5 figure. The
                          others (z-factor normality, conditional normality,
                          ECL illustration, calibration uncertainty) are
                          diagnostics, run on demand, not part of run_all.R
  plots_paper_descriptive.R  descriptive plots, run on demand
data/raw/                 frozen input data (see MANIFEST.md)
data/processed/           built datasets (regenerated by the pipeline)
output/                   generated results (regenerated by the pipeline)
tests/
  smoke_test.R            post-run sanity checks on the deliverables
  test_pd_es_closedform.R closed-form expected shortfall vs numerical
                          integration (self-contained)
```

## License

The source code is released under the MIT License (see `LICENSE`). The data
files under `data/raw/` are **not** covered by that license: they remain the
property of their respective providers (FRED, Caldara-Iacoviello, Baker-Bloom-
Davis, CBOE, S&P Dow Jones Indices) and are redistributed in frozen form for
replication convenience only, subject to those providers' terms of use. See
`data/raw/MANIFEST.md` for per-series provenance.

## Citation

If you use this code, please cite the paper:

> Flament, G., C. Hurlin, Q. Lajaunie and Y. Pull, "Generalized Impulse
> Responses of Portfolio Default Probabilities: A Modular Framework with an
> Application to Geopolitical Risk." Available at SSRN:
> <https://ssrn.com/abstract=7106918>, doi:10.2139/ssrn.7106918.

The Geopolitical Risk Index should be cited as Caldara, D. and M. Iacoviello
(2022), "Measuring Geopolitical Risk," *American Economic Review*, 112(4),
1194-1225.
