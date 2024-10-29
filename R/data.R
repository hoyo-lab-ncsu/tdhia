


#' True Diagnostic Imprintome Manifest File
#'
#' A dataframe that maps illumina probe ids used in the imprintome array to CpG
#' sites in the human genome (build GRCh37).
#'
#' @format A data frame with 22819 rows and 23 columns
#' \describe{
#'   \item{Probe_ID}{Unique probe identifier with four added alphabetic and
#'   numeric characters to denote top or bottom strand (T/B), converted or opposite
#'   strand (C/O), Infinium probe type (1/2), and the number of synthesis for
#'   representation of the probe on the array (1,2,3,…,n). Denoted as IlmnID in the
#'   GenomeStudio Manifest and Probe_ID in the SeSAMe Manifest.}
#'   \item{Name}{The locus target identifier (cg, ch, rp, mu, rs) followed by an
#'   eight-digit code that relates to the probe sequence. If an eight-digit code
#'   has not yet been generated, standard genomic coordinates follow the locus
#'   target identifier.}
#'   \item{U}{For Infinium I bead types, this is the Address ID for the probe
#'   specific for the A allele, which is the unmethylated allele. For Infinium II
#'    bead types, the Address ID for the probe used for both A and B alleles (ie,
#'     AddressB_ID and AlleleB_ProbeSeq columns are empty). Denoted as
#'     “AddressA_ID” in the GenomeStudio Manifest and “U” in the SeSAMe manifest}
#'   \item{AlleleA_ProbeSeq}{The sequence of the probe identified in AddressA_ID
#'    column.}
#'   \item{M}{For Infinium I bead types, this is the address ID for the probe
#'   that is specific for the B allele, which is the methylated allele. Denoted
#'    as “AddressB_ID” in the GenomeStudio Manifest and “M” in the SeSAMe
#'    manifest.}
#'   \item{AlleleB_ProbeSeq}{For the Infinium I bead type, the sequence of the
#'   probe identified in AddressB_ID column.}
#'   \item{Next_Base}{For Infinium I probes, the nucleotide immediately following
#'    the CpG. Blank for Infinium II.}
#'   \item{COLOR_CHANNEL}{For Infinium I probes, the color channel of the
#'   Next_Base signal.}
#'   \item{col}{For Infinium I probes, the color channel of the “Next_Base”
#'   signal. The red and green are abbreviated to R and G, respectively.}
#'   \item{Probe_Type}{Either cg, ch, mu, rp, or rs to denote CpG, CpH,
#'   multi-unique, repetitive element, or SNP probes. Control probes denoted in
#'    Controls section in the GenomeStudio Manifest and prefixed with “ctl” in
#'    the SeSAMe Manifest.}
#'   \item{Strand_FR}{The forward (F) or reverse (R) designation of the design
#'   strand. Strand_FR is dependent on the genome build used to prepare the array
#'    and manifest.}
#'   \item{Strand_TB}{Either top (T) or bottom (B) specifying whether the probe
#'   is positioned upstream (in smaller coordinates) or downstream (in greater
#'   coordinates) of the target base. Strand_TB is not dependent on the genome
#'   build used to prepare the array and manifest.}
#'   \item{Strand_CO}{Either converted (C) or opposite (O) depending on whether
#'   the probe queries the original bisulfite converted DNA strand or the
#'   opposite strand that results from amplification of the originally converted
#'   DNA stand. Strand_CO is not dependent on the genome build used to prepare
#'   the array and manifest.}
#'   \item{Infinium_Design}{Numeric Version of Infinium_Design_Type}
#'   \item{Infinium_Design_Type}{Infinium I (2 probes/locus) or Infinium II
#'   (1 probe/locus).}
#'   \item{CHR}{Chromosome containing the CpG (Build 37).}
#'   \item{MAPINFO}{Chromosomal coordinates of the CpG (Build 37).}
#'   \item{Species}{Species that samples were drawn from}
#'   \item{Genome_Build}{Genome Build referenced for this manifest.}
#'   \item{Source_Seq}{The original, genomic sequence used for probe design
#'   before bisulfite conversion.}
#'   \item{Forward_Sequence}{Plus (+) strand sequences (5’-3’) flanking the
#'   target base.}
#'   \item{Top_Sequence}{Illumina’s standardized TOP strand nomenclature applied
#'   to an interrogated dinucleotide site. e.g. CpG, CpH}
#'   \item{Rep_Num}{Used to distinguish multiple assays that target the same
#'    genomic site.}
#' }
#' @source True Diagnostic
"manifest_v1A2"



#' Mapping between CpG site IDs and Imprint Control Region IDs
#'
#' A table that maps CpG sites to specific imprint control regions.
#'
#' @format A data frame with 22819 rows and 23 columns.
#' \describe{
#'   \item{X}{   }
#'   \item{CpG_chr}{   }
#'   \item{CpG_start}{  }
#'   \item{CpG_stop}{  }
#'   \item{ICR_chr}{  }
#'   \item{ICR_start}{  }
#'   \item{ICR_stop}{  }
#'   \item{CpG_Probe}{  }
#'   \item{ICR_id}{  }
#'   \item{CpG_id}{  }
#' }
#' @source https://doi.org/10.1080/15592294.2022.2091815
"mapping_cpg_icr_ids"


#' Imprintome probe design scores
#'
#' A named list of passed and failed design scores for imprintome probes.
#'
#' @format A named list.
#' \describe{
#'   \item{Design_Score_Summary}{   }
#'   \item{Pass_Score_Threshold}{   }
#'   \item{Fail_Score_Threshold}{   }
#'   \item{Failed_Designs}{   }
#' }
#' @source none
"design_scores"



#' Sigset Model
#'
#' An example imprintome chip dataset from the SeSame package.
#'
#' @format none
#' \describe{
#'   \item{sigset_template}{   }
#'   \item{sigset_mask}{   }
#'   \item{sigset_values}{   }
#' }
#' @source none
"sigset_model"





#' Manifest V1A2, probe design scores included
#'
#' manifest table for custom imprintome array with illumina design scores added.
#'
#' @format none
#' \describe{
#'   \item{sigset_template}{   }
#'   \item{sigset_mask}{   }
#'   \item{sigset_values}{   }
#' }
#' @source none
"manifest_v1A2_design_scores"



