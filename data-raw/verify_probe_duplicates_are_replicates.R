
#' @description This script verifies that all probes that target the same cpg site are exact 
#' replicates (same A and B allele sequences). This is useful for how the data is 
#' merged between replicates in the processing pipeline.
#' @param mft this is the internal manifest for the package, it is the original 
#' manifest with the cpg_id column added.
#' 
library(dplyr)

# Internal manifest
mft <- tdhia::manifest_v1A2_design_scores

# Find the cpg_ids with multiple probe_ids that map to them
cpg_table <- table(mft$cpg_id)
shared_cpg_ids <- names(cpg_table)[cpg_table==2]
# 1342 shared cpg_ids

# Check how many of shared probes are proper cg probes
sum(grepl(pattern = "^cg[0-9]", x = shared_cpg_ids))
# 1329/1342 are proper cg probes


# Dataframe of probes with shared cpg_ids
shared_mft <- filter(mft, cpg_id %in% shared_cpg_ids)
# 1392 probe_ids target a shared cpg_id

# Filter out missing sequence for B, these are control probes
shared_mft <- filter(shared_mft, !is.na(AlleleB_ProbeSeq))

# Grab one shared cpg+id for testing
# shared_example <- filter(shared_mft, cpg_id == shared_mft$cpg_id[1])
# shared_example
# mft$cpg_id


# Define function that tests if a dataframe of probe_ids are replicates or not 
# based on whether the probe A allele and B allele sequences match.
seq_test = function(a,b) {
  # Test that first A sequence matches all A OR B sequences
  a_match <- a[1] == a | a[1] == b
  # Test that first B sequence matches all B OR A sequences
  b_match <- b[1] == b | b[1] == a
  all_matches <- a_match & b_match
  return(all_matches)
}

# Use dplyr to group probes by shared cpg_id, perform test
replicate_test <- shared_mft %>% group_by(cpg_id) %>% 
  summarize(mean_U = mean(U), is_replicate = all(seq_test(AlleleA_ProbeSeq, AlleleB_ProbeSeq)))
# See how many probes failes test (they all passed)
table(replicate_test$is_replicate)
# 695 cpg_ids have replicate probes (all entries in replicate_test)
