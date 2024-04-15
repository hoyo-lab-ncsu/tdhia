# TruDiagnostic human imprintome array (tdhia)

![banner image](inst/banner.jpg)
<br><br>   
   
**Description**: this software is designed to analyze custom 450k infinium methylation microarrays that characterize the **human imprintome**. The imprintome is the set of imprinted control regions within the genome where where a collection of proximal CpG sites exhibit fixed monoallelic status in methylation state in a parent-of-origin–dependent fashion. These regions are thought to strongly influence expression of adjacently located genes and may play a role in disease susceptability throughout the human lifespan.

**Pipeline Summary**: takes Imprintome array raw data files (.IDAT extension) and metadata associated with the study, produces a table that shows the mean % methylation for each ICR region. 
<br><br> 
_Note: IDAT files are assumed to be named by patient_ID, and they should exist in pairs, one file for red channel flourescence, and one file for green channel flourescenet. So for the theoretical patient_ID of **"28137192_R7398"**, the IDAT files would be named:_
1. **28137192_R7398**_red.IDAT
2. **28137192_R7398**_grn.IDAT  
_> suffix might be lower or upper case in filenames, code ignores cased for suffix._
<br><br>
    
## Initial Setup and Installation
  
1. Download and install the latest version of [R](https://www.r-project.org/), [R Studio](https://posit.co/download/rstudio-desktop/). If on windows, also install [R Tools](https://cran.r-project.org/bin/windows/Rtools/).
2. Download this repository to your computer by either:
   1. Clicking the green **[Code]** button at the top of this page and download repository as a zip file (not recommended because you have to redownload for any updates to the code).
   2. Cloning this repository using a github client (**recommended option**, cloning allows you to sync the code on your computer with any updates added here). For windows and Mac, you can install [Github Desktop](https://desktop.github.com/). Once installed clock the green **[Code]** button on this page, and select "Open in Github Desktop", a dialogue box will appear to "Clone" this repository. It will download the repository to a default location (usually in ~/Documents/Github/tdhia/). You can download updates to this repository by clicking "Fetch Origin" and "Sync" within github desktop, and all the files will be automatically updated.
3. Clone or download the github that has **example scripts** for this package, using the same steps as previously: [tdhia_scripts](https://github.com/bacorli2/tdhia_scripts).
4. **Set Working Directory to TDHIA**: Open RStudio, click [ **Session** ] >> [ **Set Working Directory** ] >> [ **Browse** ] to the location of the tdhia repository (step 2, in this example case located in ~/Documents/Github/tdhia/).
5. **Install Devtools**: run in the command line: >install.packages("devtools")
6. **Install BiocManager**:

Example scripts (need to convert to vignettes): https://github.com/bacorli2/tdhia_scripts
    
## Technical Overview of Current Pipeline


### Summary:
1. **Importation and Processing**: Import IDAT files, process with standard sesame pipeline using the True Diagnostic Imprintome array manifest. Calculate matrix of beta values (probe_beta_matrix) and sesame detection p-values (probe_beta_pval). Dimensions for both: probe_id x samples. _Function: load_idata_to_probes.R_.  
2. **Filter probes**: select for relevant probes that also yield a high quality signal. Function: filter_probes.R.  
   1. Probes are discarded if they do not map to a distinct cpg site, (need specific genomic location and suffix with “cg*” for probe name).  
   2. Set any beta values to NA if the **sesame signal detection p-value is greater than 0.2**.  
   3. Discard probes that fail the 0.2 p-value threshold more than 20% of samples (an entire row of the probe_beta_matrix is discarded if the p-value > 0.2 for more than 20% of the samples). These are discarded because the probe is deemed unreliable.   
3. **Convert to CpG measurements**: the probe_beta_matrix is converted to a cpg_beta_matrix with calculating the mean beta value when multiple probes map to the same cpg site (NA values disregarded). So, if multiple probes map to a single cpg site, we calculate the mean beta value for each patient, if there are any NA values, we disregard them. If there are only NA values, then we return NA for this step. _Dimensions for cpg_beta_matrix_: cpg site id x patients. Any probes that map to a cpg site that does not map to an ICR is also discarded. _We typically get coverage for 7000 - 9000 cpg sites after probe filtering, out the total ~10,000 cpg sites. Function: convert_probes_to_cpgs.R_.  
4. **Option: Convert to ICR measurements**: use the same simple averaging, we can convert the cpg_beta_matrix to an icr_beta_matrix (for statistical analysis). I don’t typically use this because this would assume that all the mappings of cpg sites to icr sites on a biological level are correct (or else true signal will quickly get masked with any mistakes). I usually perform any analysis at the cpg level and then look up which ICR sites they belong to for summarizing results. _We typically get coverage for 1030 of the 1088 total ICR sites. Function: convert_cpgs_to_icrs.R_.  
5. **Option: Threshold the beta values for hemi-methylation state**. We can convert the _cpg_beta_matrix_ or _icr_beta_matrix_ from a continuous variable to binary by a naïve threshold within some range (such as 0.5 +- 0.15).  
6. **Statistical Analysis**: Use the _cpg_beta_matrix_ and the _study metadata_ (additional data for each sample/patient) to test for associations of individual cpg sites and study metadata variables.  
   1. **Define models for analysis**, the _cpg_beta_matrix_ can either be used as a response variable or predictor depending on the situation. We produce a series of models where the fitting is parallelized over the different cpg sites (since they each need a separate  model). If the response variable is **continuous** (i.e. beta measurements for each CPG or ICR site), we use a **gaussian model**, if the response variable is **binary** (thresholded for hemi-methylation state), we use a **binomial mode**.  
      1. **Example**: testing for association between methylation state of cpg sites in offspring versus the heavy metal exposure of the parents. **Models**: cpg_ids ~ metal_exposure + other_predictors + confounders.  
      3. **Example**: testing for association with methylation state of cpg sites in offspring versus their birth weight. **Models**: birth_weight ~ cpg_ids + other_predictors + confounders.  
   2. **Filling in NA values (imputation)**: The GLM package cannot handle any NA values in the dataset, so we must fill them in (it would initially appear that glm() can handle NA values, but after digging deep into the discussion threads it looks like it cannot do so, the function discards entire rows/columns if there are NA values even if you tell it differently). We do this with imputation with the MICE package. Basically, we take our statistical model above, the cpg_beta_matrix, and have MICE produce multiple imputations and use them to produce reasonable values to fill in any missing datapoints. I don’t know if this is a proper approach or not.  
   3. **Fit models and calculate statistics**: use glm package to fit a series of models, one for each cpg or icr site (parallelized for multicore support).   
   4. **P-value Correction**: the p-values for predictors are corrected for FDR. We currently correct the cpg level analysis with n~= 1000 (number of ICR sites covered by cpg sites that passed QC in step 2c). A more rigorous correction would be n~= 9000 (number of cpg sites analyzed, a proper FDR correction). Maybe the n=1000 correction is too liberal/ forgiving.  
   5. **Output**: list of cpg sites with statistically significant FDR corrected p-values, and what ICR sites they map to.  

Note: The words matrix specially refers to dataframes in R (since that is an R specific term).

    
## Study Metadata File Requirements
  
The study metadata is an xlsx file that has a couple of require columns for the pipeline. At this point, the rest of the columns are ignored, but in the future they will be preserved in the output file for statistical analysis. Currently, you will have to merge the output data with the study metadata to include all the mediators and cofactors for the statistical analysis.
  
_Required columns_:  
1. **source**: a string with no spaces that describes the source of the study (examples; "SHIP" or "SMKE").  
2. **patient_study_id**: this is a simple integer ID for the patient that is specific to the study/source.  
3. **Patient_ID**: the original medical record number of the patient (this is used to open the pair of IDAT files mentioned previously, it is the base name of the IDAT files).  
<br><br>
    
## Example summary output
  
```
Example outputs
Probe filtering:
Probe manifest: manifest file has a total of 22819 probes.
Probe filter: discarding 46% probes ( 12381/ 22819) in dataset  b/c they don't
                map uniquely to the genome or a CpG site.
                10438 probes now remain.
Probe filter: 20% of probe beta measurements ( 397097/ 1993658) failed
                 the signal max p-value threshold of 0.20, setting them to NA.
Probe filter: 22% of probes ( 2320/ 10438) had a signal p-value
                fail rate above threshold of 25%, setting all beta values to NA
                for those probes.
Probe filter: discarded 22% probes ( 2320/ 10438) because all measurements are now NA.
              8118 probes now remain. 

Cpg Filtering:
CpG filter: 8118 of probes in dataset mapped to 7174 unique CpG sites id data.
CpG filter: discarded 0% of CpG sites (0/ 7174) because they
      do not map uniquely to the genome. 7174 CpG sites remain.
CpG filter: discarded 7% of CpG sites (485/ 6689) because they
      do not map to an ICR. 6689 CpG sites remain.

ICR Filtering
ICR Filter: discarded 1% of ICRs ( 11/ 1088) b/c their signal fail rate was > 20%.
              1029 ICRs still remain.
```              
