


# Set here path
here::i_am("README.md")

devtools::load_all()

# Two dataframes, compares based on column
match_report <- function(x,y, fname) {
  cat("Frequency of entries of X within Y:\n")
  x_in_y <- sapply(X = x[[fname]],FUN = function(k) length(which(y[[fname]]==k)))
  print(table(x_in_y))

  cat("Frequency of entries of X within Y:\n")
  y_in_x <- sapply(X = y[[fname]],FUN = function(k) length(which(x[[fname]]==k)))
  print(table(y_in_x))

  x$Design.Score = rep(NA, nrow(x))
  for (n in 1:nrow(x)) {
    temp <- y$Design.Score[which(x[[fname]][n] == y[[fname]])]
    if (length(temp!=0))  x$Design.Score[n] <- temp
  }

  return(x)
}


# Load manifest and discard probes that don't map to genome
#______________________________________________________________________________
mft <-tdhia::manifest_v1A2
# Keep probes that map to a non zero genome position
# mft = mft[!is.na(mft$MAPINFO)& mft$MAPINFO>0,]
mft$cpg_id <- stringr::str_replace(string = mft$Probe_ID, "_.{4}$", "")
# Keep probes that map to a cg site
# mft = mft[grepl(pattern = "^cg[0-9]",x = mft$Probe_ID),]
# Add ICR mapping to mft
icr_mapping <- tdhia::mapping_cpg_icr_ids
mft$icr_id = sapply(mft$cpg_id, function(x) icr_mapping$ICR_id[which(x==icr_mapping$CpG_id)][1])
# Removes probes that don't match to ICRs
# mft <- mft[!is.na(mft$icr_id),]


# Import True Diagnostic Design Scores, Version 1
#_______________________________________________________________________________
df_pass <- read.csv(here::here("data-raw", "TruDx_AllDesigns_Pass_Score_Threshold.csv"))
df_pass$cpg_id <- stringr::str_replace(string = df_pass$Assay_Design_Id, "_.{3}$", "")
df_fail <- read.csv(here::here("data-raw", "TruDx_AllDesigns_Fail_Score_Threshold.csv"))
df_fail$cpg_id <- stringr::str_replace(string = df_fail$Assay_Design_Id, "_.{3}$", "")
df_pass_fail <- rbind(df_pass, df_fail)

# Count how many probes match to unique and shared cpg sites
# table(table(mft$cpg_id))
new_mft1 <- match_report(mft, df_pass, fname = "cpg_id")
new_mft2 <- match_report(mft, df_pass_fail, fname = "cpg_id")
# new_mft <- match_report(mft, design_scores$Pass_Score_Threshold,
#                         fname = "Top_Sequence")



# Import True Diagnostic Design Scores, Version 2
#_______________________________________________________________________________
design_scores_v2 <- read.csv(here::here("data-raw", "TruDx_Imprintome_Rnd3_GRCh38_FinalDesigns_v2_plus_DesignScore.csv"))
design_scores_v2$cpg_id <- stringr::str_replace(string = design_scores_v2$Assay_Design_Id, "_.{3}$", "")
design_scores_v2$Design.Score <- as.numeric(design_scores_v2$Design.Score)


new_mft <- left_join(x = mft, y=select(design_scores_v2, c(cpg_id, Design.Score)),
          by = "cpg_id", keep = FALSE, multiple = "first",relationship = "many-to-one",
          na_matches = "never",unmatched = "drop")

# Confirm that new mft with design.scores added has the same data as old mft
# for all other columns
all((mft==new_mft[,1:25]) | is.na(mft))

# Add new manifest to package data
manifest_v1A2_design_scores <- new_mft
# usethis::use_data(manifest_v1A2_design_scores)






