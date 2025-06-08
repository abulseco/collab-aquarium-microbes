# Extract sequences per ASV
# GOT Project

# SETUP----
## libraries----
library(Biostrings); library(here)

## Files----
here()
seqtab.nochim.16S <- readRDS("seqtab_nochimeras_16S.rds")
seqtab.nochim.18S <- readRDS ("seqtab_nochimeras_18s.rds")

# Export for 16S----
## Extract the ASV sequences (column names)
asv_seqs_16S <- colnames(seqtab.nochim.16S)  # These should be DNA sequences

## Generate ASV IDs
asv_ids_16S <- paste0("ASV", seq_along(asv_seqs_16S))  # e.g., ASV1, ASV2, ...

## Create a mapping table of ASV ID to sequence
asv_map_16S <- data.frame(
  ASV_ID = asv_ids_16S,
  Sequence = asv_seqs_16S,
  stringsAsFactors = FALSE
)

## Write to CSV
write.csv(asv_map_16S, "INTERMEDIATE-FILES/ASV_ID_sequence_table_16S.csv", row.names = FALSE)

# Export for 18S----
## Extract the ASV sequences (column names)
asv_seqs_18S <- colnames(seqtab.nochim.18S)  # These should be DNA sequences

## Generate ASV IDs
asv_ids_18S <- paste0("ASV", seq_along(asv_seqs_18S))  # e.g., ASV1, ASV2, ...

## Create a mapping table of ASV ID to sequence
asv_map_18S <- data.frame(
  ASV_ID = asv_ids_18S,
  Sequence = asv_seqs_18S,
  stringsAsFactors = FALSE
)

## Write to CSV
write.csv(asv_map_18S, "INTERMEDIATE-FILES/ASV_ID_sequence_table_18S.csv", row.names = FALSE)
