


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


#' Study data for alzheimer brain samples
#'
#'
#'
#' @format none
#' \describe{
#'   \item{patient_id}{base name of truediagnostic files   }
#'   \item{group}{either control group or those with alzheimers}
#'   \item{race}{race of patient (black or white)}
#' }
#' @source none
"df_alzheimer_brain_study_data"



#' Alzheimer’s Disease Brain Methylation Data
#'
#' This dataset contains the probe-level beta values for 16 Alzheimer’s Disease brain samples.
#'
#' @name ad_brain_n16_probe_beta_matrix
#' @docType data
#' @keywords datasets
#' @usage NULL
#' @format A matrix with rows as probes and columns as samples. The columns represent the samples, and each column contains raw beta values for that sample.
#' \describe{
#'   \item{Column names}{Sample IDs corresponding to the IDAT files.}
#'   \item{Metadata}{
#'     The following table describes the disease status and race of each sample:
#'     \tabular{lll}{
#'       \strong{Sample Name} \tab \strong{Disease} \tab \strong{Race} \cr
#'       207344530004_R01C01 \tab Control \tab black \cr
#'       207344530004_R03C01 \tab Control \tab black \cr
#'       207344530004_R05C01 \tab Alzheimer_Disease \tab black \cr
#'       207344530004_R07C01 \tab Control \tab black \cr
#'       207344530004_R09C01 \tab Alzheimer_Disease \tab black \cr
#'       207344530004_R11C01 \tab Control \tab black \cr
#'       207344530004_R02C01 \tab Alzheimer_Disease \tab black \cr
#'       207344530004_R04C01 \tab Alzheimer_Disease \tab black \cr
#'       207344530004_R06C01 \tab Control \tab white \cr
#'       207344530004_R08C01 \tab Alzheimer_Disease \tab white \cr
#'       207344530004_R10C01 \tab Control \tab white \cr
#'       207344530004_R12C01 \tab Control \tab white \cr
#'       207344530004_R01C02 \tab Alzheimer_Disease \tab white \cr
#'       207344530004_R03C02 \tab Control \tab white \cr
#'       207344530004_R05C02 \tab Alzheimer_Disease \tab white \cr
#'       207344530004_R07C02 \tab Alzheimer_Disease \tab white
#'     }
#'   }
#' }
#' @source 16 Alzheimer's brain samples with matching WGBS that used in array design paper
"ad_brain_n16_probe_beta_matrix"




#' WGBS ratio data from brain tissue from alzheimer study
#'
#' This dataset contains ratio (methylation) values fro whole genome bisulfite sequencing.
#'
#' @name alzheimer_study_wbgs_shared_cpg
#' @docType data
#' @keywords datasets
#' @usage NULL
#' @format A matrix with rows as probes and columns as samples. The columns
#' represents the samples, and each column contains the ratio values for that sample.
#' \describe{
#'   \item{Column names}{Sample IDs corresponding to the IDAT files.}
#'   \item{Metadata}{
#'     The following table describes the disease status and race of each sample:
#'     \tabular{lll}{
#'       \strong{Sample Name} \tab \strong{Disease} \tab \strong{Race} \cr
#'       207344530004_R01C01 \tab Control \tab black \cr
#'       207344530004_R03C01 \tab Control \tab black \cr
#'       207344530004_R05C01 \tab Alzheimer_Disease \tab black \cr
#'       207344530004_R07C01 \tab Control \tab black \cr
#'       207344530004_R09C01 \tab Alzheimer_Disease \tab black \cr
#'       207344530004_R11C01 \tab Control \tab black \cr
#'       207344530004_R02C01 \tab Alzheimer_Disease \tab black \cr
#'       207344530004_R04C01 \tab Alzheimer_Disease \tab black \cr
#'       207344530004_R06C01 \tab Control \tab white \cr
#'       207344530004_R08C01 \tab Alzheimer_Disease \tab white \cr
#'       207344530004_R10C01 \tab Control \tab white \cr
#'       207344530004_R12C01 \tab Control \tab white \cr
#'       207344530004_R01C02 \tab Alzheimer_Disease \tab white \cr
#'       207344530004_R03C02 \tab Control \tab white \cr
#'       207344530004_R05C02 \tab Alzheimer_Disease \tab white \cr
#'       207344530004_R07C02 \tab Alzheimer_Disease \tab white
#'     }
#'   }
#' }
#' @source 16 Alzheimer's brain samples with matching WGBS that used in array design paper
"alzheimer_study_wbgs_shared_cpg"



#' WGBS ratio data from brain tissue from alzheimer study
#'
#' This dataframe contains sample_id information for wgbs files.
#' @name alzheimer_study_wbgs_shared_cpg_metadata
#' @docType data
#' @keywords datasets
#' @usage NULL
#' @format A dataframe with sample ID for wgbs data
#' \describe{
#'   \item{ID}{ ICR ID (see 10.1080/15592294.2022.2091815).}
#'   \item{FileName}{Filename of WGBS outputfile.}
#'   \item{Treatment}{Factor that partitions patients based on race (black/white)
#'    and disease state (control/ AD).}
#'   \item{sampleID}{ Study ID of participant.}
#'
#'   \item{X}{ Different study/sample ID (not sure what this is).}
#'   \item{Sample}{ Another sample ID that match sampleID with a character code
#'   that denotes the study the patient is from.}
#'   \item{Chip}{ ID of imprintome chip.}
#'   \item{Beadchip}{ ID of imprintome bead chip.}
#'
#'   \item{Row_Col}{ Location of beadchip on the imprintome array.}
#'   \item{Notes}{ Misc Notes.}
#'   \item{Patient.ID}{ Patient ID for the imprintome raw data files (Each ent.}
#'   \item{wgbs_base_names}{Sample IDs corresponding to the IDAT files.}
#'
#'   \item{AD.Case.Control}{ Factor that divides patients between control and disease group.}
#'   \item{Ethnicity}{ Ethnicity of patient, either non-hispanic black (NHB) or non-hispanic white (NHW).}
#'   \item{Age}{ Age of participant in years.}
#'   \item{Sex}{ Biological sex of participant.}
#'
#'   \item{Brain.Region}{ Location of brain for tissue sample.}
#'   \item{Time}{ Timepoint (of sample?).}
#'   \item{AD.Neuropathologic.changes}{ Severity of alzheimer's disease (None, Low, Intermediate, High.}
#'   \item{race}{ Same as ethnicity.}
#'
#'   \item{is_control}{ Boolean variable that indicates control patient.}
#'   \item{disease_state}{ Boolean variable, where true indicate presence of alzheimer.}
#'   \item{study_group}{ Same as treatment.}
#' }
#' @source 16 Alzheimer's brain samples with matching WGBS that used in array design paper
"alzheimer_study_wbgs_shared_cpg_metadata"


#' Zinc finger data for imprintome ICRs
#'
#' This dataframe lists any zinc finger sites close to each ICR
#'
#' @name imprintome_icr_zinc_finger
#' @docType data
#' @keywords datasets
#' @format A data-frame that lists up to five of known zinc finger locations
#' that is close to the genomic location of the ICR.
#' \describe{
#'   \item{icr}{ ICR ID (see 10.1080/15592294.2022.2091815).}
#'   \item{location}{ Genomic location of ICR (range of genomic coordinates).}
#'   \item{zf1}{ ID of zinc finger site, #1}
#'   \item{zf2}{ ID of zinc finger site, #2}
#'   \item{zf3}{ ID of zinc finger site, #3}
#'   \item{X}{   ID of zinc finger site, #4}
#'   \item{X.1}{ ID of zinc finger site, #5}
#' }
#' @source 16 Alzheimer's brain samples with matching WGBS that used in array design paper
"imprintome_icr_zinc_finger"

#' WGBS ratio data from brain tissue from alzheimer study
#'
#' This dataframe lists nearest transcripts for all ICRs
#' @name imprintome_icr_nearest_transcripts
#' @docType data
#' @keywords datasets
#' @format A dataframe with sample ID for wgbs data
#' \describe{
#'   \item{ID}{ ICR ID}
#'   \item{Genomic.Coordinates}{ Genomic location of ICR (range of genomic coordinates).}
#'   \item{Parental.Origin.of.Methylation}{ Indicates whether there is expeiremntal evidence of parental origin of methylation for that specific ICR (presence of "P")}
#'   \item{Nearest.Transcript}{ ID of nearest transcript to the geonimic location of the ICR (if one exists).}
#'   \item{Distance.to.Nearest.Transcript}{ Genomic distance of ICR to nearest transcript.}
#' }
#' @source 16 Alzheimer's brain samples with matching WGBS that used in array design paper
"imprintome_icr_nearest_transcripts"

#' WGBS ratio data from brain tissue from alzheimer study
#'
#' This dataframe lists nearest transcripts for ICRs that have experimental
#' evidence of gametic origin of methylation
#'
#' @name imprintome_icr_gametic_nearest_transcripts
#' @docType data
#' @keywords datasets
#' @format A dataframe with sample ID for wgbs data
#' \describe{
#'   \item{ID}{ ICR ID}
#'   \item{Genomic.Coordinates}{ Genomic location of ICR (range of genomic coordinates).}
#'   \item{Parental.Origin.of.Methylation}{ Indicates whether there is expeiremntal evidence of parental origin of methylation for that specific ICR (presence of "P")}
#'   \item{Nearest.Transcript}{ ID of nearest transcript to the geonimic location of the ICR (if one exists).}
#'   \item{Distance.to.Nearest.Transcript}{ Genomic distance of ICR to nearest transcript.}
#' }
"imprintome_icr_gametic_nearest_transcripts"



#' This dataframe lists the ICRs with evidence of gametic origin of methylation
#'
#' @name imprintome_gametic_icrs
#' @docType data
#' @keywords datasets
#' @format A dataframe with sample ID for wgbs data
#' \describe{
#'   \item{ch}{ Chromosome number.}
#'   \item{start}{ genomic position start of ICR}
#'   \item{end}{ genomic position end of ICR}
#'   \item{icr_id}{ ICR identifier}
#' }
"imprintome_gametic_icrs"


#' List of epicV2 cpg sites that are shared with imprintome
#'
#' @name epicv2_shared_cpg_sites
#' @docType data
#' @keywords datasets
#' @format A dataframe with sample ID for wgbs data
#' \describe{
#'   \item{ch}{ Chromosome number.}
#'   \item{start}{ genomic position start of ICR}
#'   \item{end}{ genomic position end of ICR}
#'   \item{icr_id}{ ICR identifier}
#' }
"epicv2_shared_cpg_sites"


