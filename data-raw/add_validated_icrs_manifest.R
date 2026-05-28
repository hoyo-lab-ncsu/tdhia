


# Set here path
here::i_am("README.md")


# read.csv(file = )


df_valid <- as.data.frame(data.table::fread(here::here("data-raw","validated_ICRs_11152024.bed")))
colnames(df_valid) <- c("chr", "start", "finish", "icr")


mft <- tdhia::manifest_v1A2_design_scores


# Same chromosome and withint the start and end range
is_validated <- sapply(1:nrow(mft), function(x) {any((
  mft[x,]$CHR==df_valid$chr) &
    (mft[x,]$MAPINFO > df_valid$start &
       mft[x,]$MAPINFO < df_valid$finish ))})


