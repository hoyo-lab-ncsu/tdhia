


here::i_am("README.md")
devtools::load_all()

# Load each patient file, filter for overlap between manifest and wgbs

# Define path where raw WGBS tables are located
data_dir <- paste0(dirname(here::here()), "/data/tdhia_brain/TruDx_17brain_black_and_white_case_control_dereje")


bed_list <- list.files(path = data_dir, pattern = ".bed$")


# Load metadata for imprintome array and the location of cpg sites
mft <- tdhia::manifest_v1A2_design_scores
cpg_atlas <- tdhia::mapping_cpg_icr_ids
cpg_atlas$chr_cpg_start = paste0(cpg_atlas$CpG_chr, "::", cpg_atlas$CpG_start)
unq_cpg_atlas <- dplyr::distinct(dplyr::select(cpg_atlas, -c(CpG_Probe,X)))

# The WGBS data has unique sites with every file
# Therefore the stored data from each patient should be based on the sites in the
#  imprintome manifest, and not based on the sites from each patient.


df_ratio_list <- df_count_list <- list()
for (n in seq_along(bed_list)) {
  cat(sprintf("Processing file %.0f/%.0f\n", n, length(bed_list)))
  # Read in BED file for WGBS data
  df_wgbs <- as.data.frame(data.table::fread(paste0(data_dir, "/", bed_list[n])))
  df_wgbs$chr_cpg_start = paste0("chr", df_wgbs[["#chrom"]], "::", df_wgbs$start)


  # Filter rows of WGBS that match cpg sites in atlas
  df_wgbs$in_imprintome <- df_wgbs$chr_cpg_start  %in% unq_cpg_atlas$chr_cpg_start
  sub_df_wgbs <- df_wgbs[df_wgbs$in_imprintome,]


  # Add cpg_id to unq_cpg_atlast dataframe so that cpg side have same order and
  #  IDs across patients
  merged_df_wgbs <-
    dplyr::left_join(x = dplyr::select(unq_cpg_atlas, c("chr_cpg_start", "CpG_id")),
              y = dplyr::select(sub_df_wgbs, "chr_cpg_start", "start", "end", "ratio","totalC"), by = "chr_cpg_start",
              keep = FALSE, multiple = "first")


  # Extract the ratio data into a data.frame column (with rownames as cpg_id)
  df_ratio_list[[n]] <- data.frame("ratio"= merged_df_wgbs$ratio, row.names = merged_df_wgbs$CpG_id)
  df_count_list[[n]] <- data.frame("totalC"= merged_df_wgbs$totalC, row.names = merged_df_wgbs$CpG_id)

}
df_ratio <- do.call(cbind, df_ratio_list)
df_count <- do.call(cbind, df_count_list)

base_names <- sapply(strsplit(bed_list, "_"), function(x) x[1])

# alzheimer_study_wbgs_shared_cpg_ratio <- df_ratio
# colnames(alzheimer_study_wbgs_shared_cpg_ratio) <- base_names
# 
# alzheimer_study_wbgs_shared_cpg_count <- df_count
# colnames(alzheimer_study_wbgs_shared_cpg_count) <- base_names


colnames(df_ratio) <- base_names
colnames(df_count) <- base_names


alzheimer_study_wbgs_shared_cpg = list()
alzheimer_study_wbgs_shared_cpg$ratio <- df_ratio
alzheimer_study_wbgs_shared_cpg$count <- df_count
  

usethis::use_data(alzheimer_study_wbgs_shared_cpg, overwrite = TRUE)
# devtools::load_all()

