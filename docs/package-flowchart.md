# tdhia Package Flow Chart

This is a simplified local review map for the `tdhia` R package.

Primary visual files:

- `docs/package-flowchart.pdf`: printable one-page Graphviz flow chart.
- `docs/package-flowchart.svg`: browser-friendly version of the same chart.
- `docs/package-flowchart.dot`: editable Graphviz source.
- `docs/package-flowchart.html`: small local viewer for the SVG.

## Main Workflow

```text
Raw IDAT files + study metadata
        |
        v
load_idata_to_probes()
        |
        v
probe_beta
        |
        v
filter_probes()
        |
        v
filtered probe_beta
        |
        v
convert_probes_to_cpgs()
        |
        v
cpg_beta
        |
        v
convert_cpgs_to_icrs()
        |
        v
icr_beta
        |
        v
Analysis options
```

## Analysis Options

- `analyze_case_control()`: t-tests per CpG or ICR with FDR adjustment.
- `analyze_association()`: parallel GLMs with optional imputation.
- `cpg_dml_test()`: LIMMA CpG-level differential methylation analysis.
- `icr_dmr_test()`: ICR-level p-value aggregation and annotation.
- `skat_icr_test()`: ICR CpG-set tests.
- `pls_icr_test()`: PLS CpG-set tests with resampling.

## Review Checks

- Confirm sample ID alignment before modeling.
- Keep object contracts stable:
  - `probe_beta_df`: probe ID x sample ID.
  - `probe_pval_df`: probe ID x sample ID.
  - `cpg_beta_df`: CpG ID x sample ID.
  - `icr_beta_df`: ICR ID x sample ID.
- Review threshold defaults before changing behavior:
  - `max_sig_pval`
  - `max_probe_fail_rate`
  - `max_patient_fail_rate`
  - `max_icr_fail_rate`
  - `n_p_adj`
- Note: no `tests/` directory was present in this checkout during review.

## Regenerate

```bash
dot -Tpdf docs/package-flowchart.dot -o docs/package-flowchart.pdf
dot -Tsvg docs/package-flowchart.dot -o docs/package-flowchart.svg
```
