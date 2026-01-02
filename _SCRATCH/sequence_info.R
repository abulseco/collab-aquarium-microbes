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

data.frame(
  mean_reads = mean(sample_sums(phy_16_beforeRare)),
  sd_reads   = sd(sample_sums(phy_16_beforeRare))
)

# Read in objects----
seqtab.nochim.16S <- readRDS("INTERMEDIATE-FILES/seqtab_nochimeras_16S.rds")
seqtab.nochim.18S <- readRDS("INTERMEDIATE-FILES/seqtab_nochimeras_18s.rds")
phy_16 <- readRDS("INTERMEDIATE-FILES/physeq_16s.rds") # before decontam
phy_18 <- readRDS("INTERMEDIATE-FILES/physeq_18s.rds") # before decontam
phy_16_decontam <- readRDS("INTERMEDIATE-FILES/ps_16S_decontam.rds") # after decontam
phy_18_decontam <- readRDS("INTERMEDIATE-FILES/ps_18s_decontam.rds") # after decontam
phy_16_final_norare <- readRDS("INTERMEDIATE-FILES/phy-f-16_beforeRare.rds")
phy_18_final_norare <- readRDS("INTERMEDIATE-FILES/phy-f-18_beforeRare.rds")
phy_16_final_rare <- readRDS("INTERMEDIATE-FILES/phy-f-16_FINAL.rds")
phy_18_final_rare <- readRDS("INTERMEDIATE-FILES/phy-f-18_FINAL.rds")

# To compare against Ashley's phyloseq object
# (some concern for number of reads per blank)
seqtab.nochim.16S.ashley <- readRDS("INTERMEDIATE-FILES/seqtab_nochimeras_16s_ashley.rds")

# number of reads retained after chimera removal
sum(seqtab.nochim.16S) # 682847
sum(seqtab.nochim.18S) # 551195
sum(seqtab.nochim.16S.ashley) # 705035 why different?

# To export table of reads per sample
reads_per_sample_16 <- rowSums(seqtab.nochim.16S)
View(reads_per_sample_16)
reads_per_sample_18 <- rowSums(seqtab.nochim.18S)
View(reads_per_sample_18)
reads_per_sample_16_ashley <- rowSums(seqtab.nochim.16S.ashley)
View(reads_per_sample_16_ashley) # OK all blanks had high reads! 

write.csv(reads_per_sample_16, "INTERMEDIATE-FILES/raw-reads-per-16.csv")
write.csv(reads_per_sample_18, "INTERMEDIATE-FILES/raw-reads-per-18.csv")
# Remove blanks, calculate average 

# Kerry will fill in numbers before/after taxa removal
# before and after rarefying
# put in final phyloseq objects

# after decontam
sum(otu_table(phy_16_decontam)) # 490891
sum(otu_table(phy_18_decontam)) # 509549
mean(sample_sums(phy_16_decontam))
sd(sample_sums(phy_16_decontam)) # 7791.921 +/- 3920
phy_18_decontam

sample_sums(phy_16_decontam)
sample_sums(phy_18_decontam)

mean(sample_sums(phy_18_decontam))
sd(sample_sums(phy_18_decontam)) # 6980.123 +/- 3194
# Not sure if this is with or without blanks, let's make sure it doesn't include blanks
sample_names(phy_16_decontam) # 16S DOES include blanks
sample_names(phy_18_decontam) 

# after taxa removal?
sum(otu_table(phy_16_final_norare)) # 488959
sum(otu_table(phy_18_final_norare)) # 304227
mean(sample_sums(phy_16_final_norare))
sd(sample_sums(phy_16_final_norare)) 

sample_sums(phy_18_decontam)

# After rarefying
summary(sample_sums(phy_16_final_rare)) # 1184
summary(sample_sums(phy_18_final_rare)) # 331 per sample 
# This seems wrong. Kerry will address on her end. 



