# TruDiagnostic human imprintome array (tdhia)

![banner image](inst/banner.jpg)
<br><br>   
   
**Description**: this software is designed to analyze custom 450k infinium methylation microarrays that characterize the **human imprintome**. The imprintome is the set of imprinted control regions within the genome where where a collection of proximal CpG sites exhibit fixed monoallelic status in methylation state in a parent-of-origin–dependent fashion. These regions are thought to strongly influence expression of adjacently located genes and may play a role in disease susceptability throughout the human lifespan.

**Pipeline Summary**: takes Imprintime array raw data files (.IDAT extension) and metadata associated with the study, produces a table that shows the mean % methylation for each ICR region. 
<br><br> 
_Note: IDAT files are assumed to be named by patient_ID, and they should exist in pairs, one file for red channel flourescence, and one file for green channel flourescenet. So for the theoretical patient_ID of **"28137192_R7398"**, the IDAT files would be named (suffix might be lower or upper case, code ignored cased for suffix):_
1. **28137192_R7398**_red.IDAT
2. **28137192_R7398**_grn.IDAT
<br><br>
    
## User Instructions
  
This pipeline is designed to be easy to use with basic R experience.
  
1. To use the pipeline, first open the fold for the repository
2. Open **"run_this_first.R"** file within RStudio, and execute this file as source (there is a **Source** button next to run button on upper right side of editor, don't copy and paste code into command line).
3. This code set's the working directory to the base folder of the code repository, downloads and installs some packages if required, and creates a csv file in the base folder, named **"file_paths.csv"**.
4. The user can open this csv file, and copy & paste in the required file paths in the second column for each of the variables (example csv file filled in found [here](media/example_file_paths.csv)), with one variable defined per row:  
   A. **idat_dir_paths**: full file system path to the folder containing all of the IDAT files.  
   B. **study_meta_path**: full file system path to the study metadata file (contains info about which IDAT files to process, ID numbers for patients etc., see Study Metadata File Requirements below).    
   C. **output_dir_path**: full file system path to the folder where the output files will be written.  
5. The user then opens **"process_imprintome.R"** and runs the file from source. This will:  
   A. **Process the IDAT files** and calculate the beta values from each probe.  
   B. **Average probe beta value(s)** that map to the same CpG site.  
   C. **Average CpG beta value(s)** that map to the same ICR.  
   D. Pay attention if there are any **warning messages** that test for problems within the dataset.  
<br><br>
    
## Technical Overview of Current Pipeline

1. To be added.
<br><br>
    
## Study Metadata File Requirements
  
The study metadata is an xlsx file that has a couple of require columns for the pipeline. At this point, the rest of the columns are ignored, but in the future they will be preserved in the output file for statistical analysis. Currently, you will have to merge the output data with the study metadata to include all the mediators and cofactors for the statistical analysis.
  
_Required columns_:  
1. **source**: a string with no spaces that describes the source of the study (examples; "SHIP" or "SMKE").  
2. **patient_study_id**: this is a simple integer ID for the patient that is specific to the study/source.  
3. **Patient_ID**: the original medical record number of the patient (this is used to open the pair of IDAT files mentioned previously, it is the base name of the IDAT files).  
<br><br>
    
## Output Files Description
  
1. **probe_beta_matrix.csv**: table containing beta values for each probe for each patient.  
   * _Rows_: each probe, labeled with probe_id from methylation array.  
   * _Columns_: each patient, with global_patient_id created from the study metadata file (merges the "source" and "patient_study_id" columns into a single string.  
3. **cpg_beta_matrix.csv**: table containing beta values for each cpg site for each patient.  
   * _Rows_: each cpg site, labeled with cpg_id from cpg-ICR metadata file.  
   * _Columns_: each patient, with global_patient_id.  
5. **icr_beta_matrix.csv**: table containing beta values for each ICR site for each patient.  
   * _Rows_: each icr, labeled with icr_id from cpg-ICR metadata file.  
   * _Columns_: each patient, with global_patient_id.  
7. **summary_icr_beta_matrix.csv**: outputs mean and standard deviation of beta value for each ICR across all patients in study.  
<br><br><br><br>
    
## !! Security Notice !!
1. This repository is a _developmental version and private_.
2. It is designed to **not include any HIPAA/ protected study metadata**, it only contains the general functions, and metadata associated with imprintome arrays, and scripts to process particular studies (but not the study data).
3. However this can't be _automatically gauranteed_- if a developed adds a commit with protected data at any time, it will be accessible even if the file gets deleted via version control history.
4. If you want to release this code/ make it public (to include with a publication), you need permission from Dr. Hoyo and the repository has to be carefully reviewed to make sure it is **completely devoid** of protected data.
5. Only the current commit should be released (copied to a new public github repository), and not include past version controlled code changes (included in the .git folder, do not copy that to avoid exposing all of the code changes made during development of repository).

