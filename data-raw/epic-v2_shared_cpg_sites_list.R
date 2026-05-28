

# Download Link
# Download this file and unzip
#https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/methylationepic/MethylationEPIC%20v2%20Files.zip



# BiocManager::install("IlluminaHumanMethylationEPICv2manifest")

# Read in csv file form download, moved to data-raw folder
df <-  as.data.frame(data.table::fread("data-raw/EPIC-8v2-0_A1.csv"))

# load manifest
mft = tdhia::manifest_v1A2_design_scores

# find shared sites
epicv2_shared_cpg_sites <- df%>% filter(Name %in% mft$Name)

# Save data
usethis::use_data(epicv2_shared_cpg_sites)