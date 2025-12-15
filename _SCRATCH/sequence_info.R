# Extract information about sequence counts
# GOT Project
# For SCRACTH code not final publication to extract statistics
# Also need to check Haley's thesis

# Read in phyloseq objects----
# Access old files
setwd("~/Documents/Projects/r_local/2024_GOT_local/2024_GOT-PROJ_DROPBOX/GOT-PROJ")
phy_16_final <- readRDS("FINAL-PHYLO-OBJ/google_drive/phy-f-16_FINAL.rds")
phy_16_beforeRare <- readRDS("FINAL-PHYLO-OBJ/google_drive/phy-f-16_beforeRare.rds")

# Not sure what the differnece is between final/before rare but will compare
phy_18_final <- readRDS("FINAL-PHYLO-OBJ/google_drive/phy-f-18_FINAL.rds")
phy_18_beforeRare <- readRDS("FINAL-PHYLO-OBJ/google_drive/phy-f-18_beforeRare.rds")

