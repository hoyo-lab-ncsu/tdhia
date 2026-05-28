# TruDiagnostic Human Imprintome Array (`tdhia`)

![banner image](inst/banner.jpg)

## Description

`tdhia` is an R package for analyzing custom Infinium methylation microarrays that characterize the **human imprintome**. The imprintome includes imprinted control regions (ICRs), where collections of nearby CpG sites show parent-of-origin-dependent monoallelic methylation patterns. These regions may influence nearby gene expression and disease susceptibility across the lifespan.

## Pipeline Summary

The pipeline takes raw TruDiagnostic Imprintome array IDAT files and study metadata, then produces probe-, CpG-, and ICR-level methylation summaries. The main output used for many analyses is a table or matrix of methylation beta values, including mean percent methylation for each ICR region.

A local one-page workflow diagram is available at:

- `docs/package-flowchart.pdf`
- `docs/package-flowchart.svg`
- `docs/package-flowchart.html`

### IDAT File Naming

IDAT files are expected to be named by patient/sample ID and to exist in red/green channel pairs.

For a theoretical patient ID of `28137192_R7398`, the expected files are:

1. `28137192_R7398_Red.idat`
2. `28137192_R7398_Grn.idat`

The code is generally tolerant of suffix capitalization, but consistent `Red` and `Grn` naming is recommended.

## User Instructions

### Initial Setup And Installation

1. Download and install the latest version of [R](https://www.r-project.org/) and [RStudio](https://posit.co/download/rstudio-desktop/).
2. On Windows, also install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).
3. Download this repository by either:
   - Clicking the green **Code** button on GitHub and downloading a ZIP file. This is simple, but not recommended because updates require re-downloading the repository.
   - Cloning the repository with Git or GitHub Desktop. This is recommended because it allows you to sync local code with remote updates.
4. Clone or download the companion repository with example scripts: `tdhia_scripts`.

### RStudio

1. Open RStudio.
2. Set the working directory to the `tdhia` repository:
   - **Session** > **Set Working Directory** > **Choose Directory**
3. Install `devtools`:

```r
install.packages("devtools")
```

4. Install `BiocManager`:

```r
install.packages("BiocManager")
```

5. Set the Bioconductor version:

```r
BiocManager::install(version = "3.18")
```

6. Load the package in development mode:

```r
devtools::load_all()
```

7. If dependency installation fails, try:

```r
BiocManager::valid()
update.packages()
```

Then run `devtools::load_all()` again.

8. Analyze the simulated dataset in the vignette:
   - `vignettes/intro-to-tdhia.Rmd`
   - Review the text and code chunks to understand the inputs and outputs for each pipeline function.

### Basic Analysis Script With Real Data

1. Open an example analysis script from the `tdhia_scripts` repository, such as:
   - `basic_scripts/imprintome_metal_exposure.R`
2. Keep the working directory set to the `tdhia` package repository.
3. Confirm that you have access to the relevant IDAT files and study metadata.
4. Update the script paths for your machine:
   - `idat_dir_paths`: full path to the folder containing IDAT files.
   - `output_dir_path`: full path to the folder where outputs should be saved.
   - `study_data`: full path to the study metadata file.
5. Run the script.

IDAT processing can take several minutes. Some scripts cache results, for example in a file such as:

```text
[output_dir_path]/heavymetals_placenta_probe_beta.rda
```

If you need to reprocess the raw IDAT files, delete or overwrite the cached `.rda` file. Otherwise, the analysis script may load the cached results.

## Technical Overview Of Current Pipeline

### Summary

1. **Importation and processing**

   Import IDAT files and process them with the standard Sesame pipeline using the TruDiagnostic Imprintome array manifest.

   Function: `load_idata_to_probes()`

   Main outputs:
   - `probe_beta_df`: beta values, `probe_id x sample`
   - `probe_pval_df`: Sesame detection p-values, `probe_id x sample`

2. **Filter probes**

   Select probes that are both relevant and high quality.

   Function: `filter_probes()`

   Main filtering steps:
   - discard probes that do not map to a distinct genomic CpG site
   - set beta values to `NA` when the Sesame detection p-value is greater than the threshold, commonly `0.2`
   - discard probes that fail the p-value threshold in too many samples
   - optionally discard samples/patients that fail across too many probes

   QC happens in layers:

   1. **Mapping filter**: remove probes that are not mapped to a unique genomic location or are not CpG probes.
   2. **Measurement-level signal filter**: identify individual probe/sample beta values whose detection p-value is above `max_sig_pval`.
   3. **Probe-level fail-rate filter**: remove probes that fail signal detection in too many samples.
   4. **Sample-level fail-rate filter**: remove samples/patients that fail signal detection across too many probes.
   5. **Optional probe design filters**: restrict probes by design score or ICR confidence level.

   Be careful when comparing direct `filter_probes()` calls with the wrapper `tdhia_pipeline()`. Their defaults are not identical:

   | Setting | `filter_probes()` default | `tdhia_pipeline()` default |
   | --- | --- | --- |
   | `max_sig_pval` | `0.2` | `0.2` |
   | `set_failed_betas_na` | `TRUE` | `FALSE` |
   | `max_probe_fail_rate` | `0.25` | `0.2` |
   | `max_patient_fail_rate` | `0.30` | `0.25` |
   | `min_design_score` | `NA` | `NA` |

   In practice, this means the one-step wrapper is somewhat stricter on probe and patient fail rates, but it does not set individual failed beta values to `NA` before the probe-level fail-rate filtering unless that argument is changed. If QC behavior matters for a study, set these arguments explicitly in the analysis script rather than relying on defaults.

3. **Convert to CpG measurements**

   Convert the probe beta matrix to a CpG beta matrix by averaging probes that map to the same CpG site. `NA` values are ignored when possible; if all contributing values are `NA`, the CpG value is returned as `NA`.

   Function: `convert_probes_to_cpgs()`

   Main output:
   - `cpg_beta_df`: beta values, `CpG_id x sample`

   CpGs that do not map to an ICR are usually discarded. After probe filtering, coverage is commonly around 7,000 to 9,000 CpG sites out of roughly 10,000 CpG sites.

   A CpG-level QC report can be generated after this step:

   ```r
   qc_cpg <- tdhia::tdhia_qc_report(
     data_beta = data_beta,
     output_dir_path = "tdhia_results",
     level = "cpg"
   )
   ```

   The report saves sample-level and feature-level QC metrics, plus a PDF with missingness histograms, beta-density overlays, a sample correlation heatmap, and PCA plots. When the relevant objects are available, it also adds detection p-value QC, retention counts, design score summaries, ICR coverage, replicate-probe agreement, and simple sample outlier flags. The same helper can summarize `level = "probe"` or `level = "icr"` when those matrices are available.

4. **Option: convert to ICR measurements**

   Convert the CpG beta matrix to an ICR beta matrix by averaging CpGs that map to the same ICR.

   Function: `convert_cpgs_to_icrs()`

   Main output:
   - `icr_beta_df`: beta values, `ICR_id x sample`

   This is useful for summary and regional analyses. CpG-level analysis may be preferred when preserving site-level signal is important, because averaging across an ICR can mask localized CpG-level effects.

5. **Option: threshold beta values for hemi-methylation state**

   CpG or ICR beta values can be converted from continuous values to binary hemi-methylation states using a threshold around 0.5, such as `0.5 +/- 0.15`.

6. **Statistical analysis**

   Use the CpG or ICR beta matrix and the study metadata to test associations between methylation and study variables.

   The beta matrix can be used either as a response or as a predictor, depending on the research question.

   Examples:
   - `cpg_ids ~ metal_exposure + other_predictors + confounders`
   - `birth_weight ~ cpg_ids + other_predictors + confounders`

   Model fitting is parallelized over CpG or ICR sites. Continuous responses use Gaussian models, while binary responses can use binomial models.

7. **Missing values and imputation**

   Some modeling functions cannot handle missing values directly. The pipeline includes options for removing rows with missing data or imputing missing values with `mice`. Review this choice carefully for each analysis.

8. **P-value correction**

   P-values are corrected for multiple testing. Some analyses currently adjust CpG-level results using the number of covered ICRs, while a stricter CpG-level correction may use the number of CpG sites tested. Review the comparison count before interpreting significance.

9. **Output**

   Typical outputs include significant CpG or ICR results, FDR-corrected p-values, ICR mappings, annotations, and diagnostic plots.

Note: this package often uses "matrix" in the broad R sense, where the object may be a `data.frame` with feature rows and sample columns.

## Study Metadata File Requirements

Study metadata may be supplied as CSV, XLSX, or another tabular format depending on the analysis script. It should include identifiers that can be matched to IDAT basenames or beta matrix column names.

Common required columns:

1. `source`: source/study label with no spaces, such as `SHIP` or `SMKE`.
2. `patient_study_id`: study-specific patient/sample identifier.
3. `Patient_ID`: original patient/sample ID used to locate the IDAT file pair.

Additional covariates, mediators, confounders, and batch variables should be included as needed for statistical analysis.

For beadchip correction in `cpg_dml_test()`, the metadata must also include:

- `Bead`
- `Col`
- `Row`

## Example Summary Output

```text
Example outputs

Probe filtering:
Probe manifest: manifest file has a total of 22819 probes.
Probe filter: discarding 46% probes (12381/22819) because they do not
              map uniquely to the genome or a CpG site.
              10438 probes now remain.
Probe filter: 20% of probe beta measurements (397097/1993658) failed
              the signal max p-value threshold of 0.20, setting them to NA.
Probe filter: 22% of probes (2320/10438) had a signal p-value fail rate
              above threshold of 25%, setting all beta values to NA for
              those probes.
Probe filter: discarded 22% probes (2320/10438) because all measurements
              are now NA.
              8118 probes now remain.

CpG filtering:
CpG filter: 8118 probes in dataset mapped to 7174 unique CpG sites.
CpG filter: discarded 0% of CpG sites (0/7174) because they do not map
            uniquely to the genome. 7174 CpG sites remain.
CpG filter: discarded 7% of CpG sites (485/6689) because they do not map
            to an ICR. 6689 CpG sites remain.

ICR filtering:
ICR Filter: discarded 1% of ICRs (11/1088) because their signal fail rate
            was greater than 20%.
            1029 ICRs still remain.
```
