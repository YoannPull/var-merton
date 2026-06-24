# Code structure

Repository-level entry point (run from the repository root):

```r
source("run_all.R")
```

Pipeline scripts (in execution order):

```text
code/01_process_gpr.R
code/02_build_dataset.R
code/shared/var_information_set_diagnostics.R
code/dralacbn/02_dralacbn_information_set_robustness.R
code/dralacbn/03_dralacbn_direct_channel.R
code/dralacbn/04_dralacbn_direct_channel_paper_outputs.R
code/dralacbn/07_dralacbn_short_sample_bma.R
code/dralacbn/09_dralacbn_short_sample_single_figures.R
code/dralacbn/08_dralacbn_perfect_foresight.R
code/dralacbn/10_dralacbn_state_dependence.R
code/shared/variable_definitions.R
code/shared/collate_paper_outputs.R
code/shared/make_paper_outputs.R
```

Shared functions and engines live in `code/shared/`.
