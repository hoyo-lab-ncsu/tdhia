


here::i_am("README.md")


# Load each patient file, filter for overlap between manifest and wgbs

# Define path where raw WGBS tables are located
data_dir <- paste0(dirname(here::here()), "/data/tdhia_brain/TruDx_17brain_black_and_white_case_control_dereje")


bed_list <- list.files(path = data_dir, pattern = ".bed$")


# Load metadata for imprintome array and the location of cpg sites
mft <- tdhia::manifest_v1A2_design_scores
cpg_atlas <- tdhia::mapping_cpg_icr_ids
cpg_atlas$chr_cpg_start = paste0(cpg_atlas$CpG_chr, "::", cpg_atlas$CpG_start)


n = 1

# t1 <- Sys.time()
# df_wgbs <- as.data.frame(read.table(paste0(data_dir, "/", bed_list[n])))
# t1 - Sys.time()


df_wgbs <- as.data.frame(data.table::fread(paste0(data_dir, "/", bed_list[n])))
df_wgbs$chr_cpg_start = paste0("chr", df_wgbs[["#chrom"]], "::", df_wgbs$start)


# Filter rows of WGBS that match cpg sites in atlas
df_wgbs$in_imprintome <- df_wgbs$chr_cpg_start  %in% cpg_atlas$chr_cpg_start
sub_df_wgbs <- df_wgbs[df_wgbs$in_imprintome,]


test <- merge(x = sub_df_wgbs, y = dplyr::select(cpg_atlas, c("chr_cpg_start", "CpG_id")), by = "chr_cpg_start")

chr_cpg_start

chr_cpg_start, 

start, end, ratio

