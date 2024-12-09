


here::i_am("README.md")


# Load each patient file, filter for overlap between manifest and wgbs

# Define path where raw WGBS tables are located
data_dir <- paste0(dirname(here::here()), "/data/tdhia_brain/TruDx_17brain_black_and_white_case_control_dereje")


bed_list <- list.files(path = data_dir, pattern = ".bed$")


# Load metadata for imprintome array and the location of cpg sites
mft <- tdhia::manifest_v1A2_design_scores
cpg_atlas <- tdhia::mapping_cpg_icr_ids
cpg_atlas$chr_cpg_start = paste0(cpg_atlas$CpG_chr, "::", cpg_atlas$CpG_start)

df_ratio_list = list()
for (n in seq_along(bed_list)) {
cat(sprintf("Processing file %.0f/%.0f\n", n, length(bed_list)))
# Read in BED file for WGBS data
df_wgbs <- as.data.frame(data.table::fread(paste0(data_dir, "/", bed_list[n])))
df_wgbs$chr_cpg_start = paste0("chr", df_wgbs[["#chrom"]], "::", df_wgbs$start)


# Filter rows of WGBS that match cpg sites in atlas
df_wgbs$in_imprintome <- df_wgbs$chr_cpg_start  %in% cpg_atlas$chr_cpg_start
sub_df_wgbs <- df_wgbs[df_wgbs$in_imprintome,]


# dup_flag <- duplicated(select(sub_df_wgbs, "chr_cpg_start", "start", "end", "ratio"))
# Add cpg_id to thhe wgbs results table
merged_df_wgbs <-
  left_join(x = dplyr::select(sub_df_wgbs, "chr_cpg_start", "start", "end", "ratio"),
          y = dplyr::select(cpg_atlas, c("chr_cpg_start", "CpG_id")), by = "chr_cpg_start",
           keep = FALSE, multiple = "first")


# Extract the ratio data into a data.frame column (with rownames as cpg_id)
df_ratio_list[[n]] <- data.frame("ratio"= merged_df_wgbs$ratio, row.names = merged_df_wgbs$CpG_id)

}





