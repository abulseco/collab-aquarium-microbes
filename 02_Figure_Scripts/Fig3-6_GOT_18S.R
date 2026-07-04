library(phyloseq)
library(ggplot2)
library(vegan)
library(plyr)
library(dplyr)
library(ape)
library(RColorBrewer)
library(knitr)
library(randomForest) 
library(dunn.test)
library(pairwiseAdonis)
library(mgcv)
library(emmeans)
library(plotly)
library(psych)

############load phyloseq object----

ps_obj <- readRDS("~/R/GOT_analysis18S/ps_18S_decontam.rds")  # 1526 taxa and 73 samples

# Remove Blanks

ps_obj = subset_samples(ps_obj, Project != "Blank") # 1526 taxa and 67 samples

# Remove "Craniata" which are teleosts, and "Embryophyceae" & Family != "Mollusca" removed since they are known food items
#but leave in families with NA
ps_obj = ps_obj %>% subset_taxa(is.na(Family) | !(Family %in% c("Embryophyceae", "Mollusca", "Craniata"))) #1448 taxa and 67 samples

write.csv(sample_sums(ps_obj), "sample_depths_18S.csv")

saveRDS(ps_obj, "ps_obj_18S_filt.rds")


# time series only
ps_ts = subset_samples(ps_obj, Project == "TimeSeries")  #  1448 taxa and 45 samples

# copper treatment only
ps_cu = subset_samples(ps_obj, Project == "Copper_tx")  # 1448 taxa and 22 samples

# transform relative abundance
ps_obj <- transform_sample_counts(ps_obj, function(x) x / sum(x))
ps_cu <- transform_sample_counts(ps_cu, function(x) x / sum(x))
ps_ts <- transform_sample_counts(ps_ts, function(x) x / sum(x))

# create ASV tables
otus_ps<- otu_table(ps_obj)
write.csv(otus_ps, file='otu_table18S_all.csv')

tax_table <- tax_table(ps_obj)
write.csv(tax_table, file = "tax_table18S_all.csv")


############## Random Forest Time series (salinity) ----

#random forest should be able to pull out the key ASVss that separate these two sites and accurately classify them
# https://rpubs.com/michberr/randomforestmicrobe

library(randomForest)  #version 4.7.1.2

# How many ASVs do we currently have? 
ntaxa(ps_ts)  # 1448

# Make a dataframe of training data with ASVs as column and samples as rows
predictors <- t(otu_table(ps_ts))
dim(predictors)  # 1448  45

# transpose this so # of samples is first, then # of ASVs
predictors <- t(predictors)
dim(predictors)  #45   1448  

#Salinity status
# Make one column for our outcome/response variable 
response_sal <- as.factor(sample_data(ps_ts)$Sal_status)

# Combine them into 1 data frame
rf_sal.data <- data.frame(response_sal, predictors)

#Random forest results
set.seed(2)
data_sal.classify <- randomForest(response_sal~., data = rf_sal.data, ntree = 100)
print(data_sal.classify)

##Call:
#randomForest(formula = response_sal ~ ., data = rf_sal.data,      ntree = 100) 
#Type of random forest: classification
#Number of trees: 100
#No. of variables tried at each split: 38
#
#OOB estimate of  error rate: 4.44%
#Confusion matrix:
#   15 22 31 class.error
#15  8  2  0         0.2
#22  0 27  0         0.0
#31  0  0  8         0.0


# What variables are stored in the output?
names(data_sal.classify)


### Lets make some plots of the most important variables in our model. For a classification tree, variable importance is measured by mean decrease in GINI coefficient (measure of node purity) due to that variable

# Make a data frame with predictor names and their importance
imp_sal <- importance(data_sal.classify)
imp_sal <- data.frame(predictors = rownames(imp_sal), imp_sal)

# Order the predictor levels by importance
imp_sal.sort <- arrange(imp_sal, desc(MeanDecreaseGini))
imp_sal.sort$predictors <- factor(imp_sal.sort$predictors, levels = imp_sal.sort$predictors)

# Select the top 20 predictors
imp_sal.20 <- imp_sal.sort[1:20, ]
# and top 10
imp_sal.10 <- imp_sal.sort[1:10, ]

# ggplot
ggplot(imp_sal.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important ASVs Salinity 18S rRNA") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(imp_sal.10, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important ASVs Salinity 18S rRNA") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

# What are those ASVs?
library(knitr)

otunames_sal <- imp_sal.10$predictors
row10_sal <- rownames(tax_table(ps_ts)) %in% otunames_sal
kable(tax_table(ps_ts)[row10_sal, ])

write.csv((tax_table(ps_ts)[row10_sal, ]), "RF_results_Salinity_18S.csv")


###plot random forest results

rf_sal <- subset_taxa(ps_ts, rownames(tax_table(ps_ts)) %in% c("ASV77","ASV56","ASV6","ASV10","ASV25","ASV58","ASV16","ASV12","ASV3","ASV118"))

plot_bar(rf_sal, "Description", fill = "OTU") + ggtitle("Random Forest ASVs Salinity 18S rRNA") + facet_grid(.~Sal_status, scales = "free") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_fill_brewer(type="qual", palette="PRGn")


#reorder and relabel ASVs
df_sal <- psmelt(rf_sal)

df_sal$OTU <- factor(df_sal$OTU, levels = c("ASV77","ASV56","ASV6","ASV10","ASV25","ASV58","ASV16","ASV12","ASV3","ASV118"),
                     labels = c("ASV77 - Ploeotiidae","ASV56 - NA","ASV6 - Aspidiscidae","ASV10 - Labyrinthulomycetes","ASV25 - Sorodiplophryidae","ASV58 - Aspidiscidae","ASV16 - Chilodonellidae","ASV12 - Aequoreidae","ASV3 - Diplonemidae","ASV118 - Chilodonellidae"))

ggplot(df_sal, aes(x = Description, y = Abundance, fill = OTU)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  ggtitle("B: 18S rRNA - Time Series") +
  facet_grid(. ~ Sal_status, scales = "free", space = "free_x") +
  scale_fill_brewer(type = "qual", palette = "BrBG") +
  xlab("Sample") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.size = (unit(1, "lines"))) +
  scale_y_continuous(expand = c(0, 0))+
  guides(fill = guide_legend(title = NULL))


############## Random Forest All data (DCS) ----
#all data (copper tx and time series)
# How many ASVs do we currently have? 
ntaxa(ps_obj)  # 1448

# Make a dataframe of training data with ASVs as column and samples as rows
predictors_obj <- t(otu_table(ps_obj))
dim(predictors_obj)  # 1448   67

# transpose this so # of samples is first, then # of ASVs
predictors_obj <- t(predictors_obj)
dim(predictors_obj)  # 67   1448

#Disease Condition Score
# Make one column for our outcome/response variable 
response_DCS <- as.factor(sample_data(ps_obj)$DCS)

# Combine them into 1 data frame
rf_DCS.data <- data.frame(response_DCS, predictors_obj)

#Random forest results
set.seed(2)
data_DCS.classify <- randomForest(response_DCS~., data = rf_DCS.data, ntree = 100)
print(data_DCS.classify)

##Call:
#randomForest(formula = response_DCS ~ ., data = rf_DCS.data,      ntree = 100) 
#Type of random forest: classification
#Number of trees: 100
#No. of variables tried at each split: 38
#
#OOB estimate of  error rate: 22.39%
#Confusion matrix:
#   0 1 2 3 class.error
#0 45 0 0 0         0.0
#1  5 0 0 1         1.0
#2  2 0 0 0         1.0
#3  7 0 0 7         0.5


# What variables are stored in the output?
names(data_DCS.classify)


### Lets make some plots of the most important variables in our model. For a classification tree, variable importance is measured by mean decrease in GINI coefficient (measure of node purity) due to that variable

# Make a data frame with predictor names and their importance
imp_DCS <- importance(data_DCS.classify)
imp_DCS <- data.frame(predictors_obj = rownames(imp_DCS), imp_DCS)

# Order the predictor levels by importance
imp_DCS.sort <- arrange(imp_DCS, desc(MeanDecreaseGini))
imp_DCS.sort$predictors_obj <- factor(imp_DCS.sort$predictors_obj, levels = imp_DCS.sort$predictors_obj)

# Select the top 20 predictors
imp_DCS.20 <- imp_DCS.sort[1:20, ]
# and top 10
imp_DCS.10 <- imp_DCS.sort[1:10, ]

# ggplot
ggplot(imp_DCS.20, aes(x = predictors_obj, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important ASVs DCS 18S rRNA") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(imp_DCS.10, aes(x = predictors_obj, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important ASVs DCS 18S rRNA") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

# What are those ASVs?
library(knitr)

otunames_DCS <- imp_DCS.10$predictors_obj
row10_DCS <- rownames(tax_table(ps_obj)) %in% otunames_DCS
kable(tax_table(ps_obj)[row10_DCS, ])

write.csv((tax_table(ps_obj)[row10_DCS, ]), "RF_results_DCS_18S.csv")


###plot random forest results this way

rf_DCS <- subset_taxa(ps_obj, rownames(tax_table(ps_obj)) %in% c("ASV45", "ASV23","ASV73","ASV3","ASV316","ASV90","ASV16","ASV307","ASV11", "ASV26"))

plot_bar(rf_DCS, "Description", fill = "OTU") + ggtitle("Random Forest ASVs DCS 18S rRNA") + facet_grid(.~DCS, scales = "free") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_fill_brewer(type="qual", palette="PRGn")


#reorder and relabel ASVs
df_DCS <- psmelt(rf_DCS)

df_DCS$OTU <- factor(df_DCS$OTU, levels = c("ASV45", "ASV23","ASV73","ASV3","ASV316","ASV90","ASV16","ASV307","ASV11", "ASV26"),
                     labels = c("ASV45 - Arcyriaceae", "ASV23 - Cryptocaryonidae","ASV73 - Diplonemidae","ASV3 - Diplonemidae","ASV316 - NA","ASV90 - Diplonemidae","ASV16 - Chilodonellidae","ASV307 - NA","ASV11 - Paramoebidae", "ASV26 - Gigartinales_X"))

ggplot(df_DCS, aes(x = Description, y = Abundance, fill = OTU)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  ggtitle("D: 18S rRNA - Disease Condition Status") +
  facet_grid(. ~ DCS, scales = "free", space = "free_x") +
  scale_fill_brewer(type = "qual", palette = "PRGn") +
  xlab("Sample") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.size = (unit(1, "lines"))) +
  scale_y_continuous(expand = c(0, 0))+
  guides(fill = guide_legend(title = NULL))


###### TAXA PLOTS ----
custom_col30 = c("#784511","#D2781E","#EAAB6C", "#781122", "#D21E2C","#A51876","#E43FAD","#F098D3","#185EA5","#1E78D2","#3F91E4","#6CABEA","#98C4F0","#117878","#18A5A5","#3FE4E4","#6CEAEA","#98F0F0","#117845","#18A55E","#1ED278","#3FE491","#6CEAAB","#98F0C4","#787811","#A5A518","#D2D21E","#E4E43F","#EAEA6C","#F0F098","#F7F7C5")
custom_col30a <- c("#F0A0FF", "#0075DC", "#993F00", "#4C005C", "#191919","#005C31", "#2BCE48", "#FFCC99", "#808080", "#94FFB5","#8F7C00", "#9DCC00","#E43FAD", "#C20088", "#003380", "#3F91E4", "#FFA405","#FFA8BB", "#426600","#185EA5", "#6CEAAB", "#FF0010", "#5EF1F2", "#00998F","#E0FF66", "#740AFF", "#990000", "#FFFF80", "#EAEA6C","#D2781E")

# 1. Transform the original phyloseq object to relative abundance
# did this earlier: ps_obj <- transform_sample_counts(ps_obj, function(x) x / sum(x))

# 2. Agglomerate at Family level
ps_family <- tax_glom(ps_obj, taxrank = "Family", NArm = TRUE)

# in the OTU table for 18S, family is actually in the column labeled species
ps_species <- tax_glom(ps_obj, taxrank = "Species", NArm = TRUE)

# 3. Get the top 30 families by total abundance (across all samples)
family_sums <- taxa_sums(ps_family)
top30_families <- names(sort(family_sums, decreasing = TRUE))[1:30]

species_sums <- taxa_sums(ps_species)
top30_species <- names(sort(species_sums, decreasing = TRUE))[1:30]

# 4. Prune to only the top 30 families
ps_top30_families <- prune_taxa(top30_families, ps_family)

ps_top30_species <- prune_taxa(top30_species, ps_species)

# 5. Melt for plotting
ps_melt <- psmelt(ps_top30_families)

ps_melt_sp <- psmelt(ps_top30_species)

## make barplot

# top 30 families
ps_melt$new_datetime_column <- as.POSIXct(ps_melt$DateID, format = "%m_%d_%y")
ps_melt2 <- ps_melt[order(ps_melt$new_datetime_column), ]
ps_melt2$DateID <- factor(ps_melt2$DateID, levels = unique(ps_melt2$DateID))

ggplot(ps_melt2, aes(x = DateID, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack",color = "black") +
  facet_grid(.~ Project, scales = "free_x", space = "free_x", labeller = labeller(Project = c("Copper_tx" = "Copper Treatment", "TimeSeries" = "Time Series")))+
  labs(x = "Date", y = "Relative Abundance") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_blank(),
    axis.title.x = element_text(vjust = -1),
    legend.text = element_text(size = 10),
    legend.key.size = (unit(1, "lines"))) +
  scale_y_continuous(expand = c(0, 0))+
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values=custom_col30a)+
  ggtitle("B: 18S rRNA")

# Top 30 species
ps_melt_sp$new_datetime_column <- as.POSIXct(ps_melt_sp$DateID, format = "%m_%d_%y")
ps_melt_sp2 <- ps_melt_sp[order(ps_melt_sp$new_datetime_column), ]
ps_melt_sp2$DateID <- factor(ps_melt_sp2$DateID, levels = unique(ps_melt_sp2$DateID))

plot <- ggplot(ps_melt_sp2, aes(x = DateID, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity", position = "stack",color = "black") +
  facet_grid(.~ Project, scales = "free_x", space = "free_x", labeller = labeller(Project = c("Copper_tx" = "Copper Treatment", "TimeSeries" = "Time Series")))+
  labs(x = "Date", y = "Relative Abundance") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_blank(),
    axis.title.x = element_text(vjust = -1),
    legend.text = element_text(size = 10),
    legend.key.size = (unit(1, "lines"))) +
  scale_y_continuous(expand = c(0, 0))+
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values=custom_col30a)+
  ggtitle("B: 18S rRNA")

ggplotly(plot)

ggsave("Fig_taxa18S.jpg",
       plot = plot,
       width = 10.6, height = 4,
       units = "in",
       dpi = 600)

####
# exporting otu tables to find out actual percentages, etc.
otus_ps_fam<- otu_table(ps_top30_families)
write.csv(otus_ps_fam, file='otu_table_all_fam18S_NA_rar.csv')

tax_table_fam <- tax_table(ps_top30_families)
write.csv(tax_table_fam, file = "tax_table_all_fam18S_NA_rar.csv")


otus_ps_species<- otu_table(ps_top30_species)
write.csv(otus_ps_species, file='otu_table_all_spec18S_NA_rar.csv')

tax_table_species <- tax_table(ps_top30_species)
write.csv(tax_table_species, file = "tax_table_all_spec18S_NA_rar.csv")

metadata_ps_obj <- as(sample_data(ps_obj), "data.frame")
write.csv(metadata_ps_obj, file = "metadata_ps_obj18S.csv")


###### ALPHA DIVERSITY ----
# this was done in a separate file
# rarefied to the lowest sample depth
ps_obj <- rarefy_even_depth(ps_obj, rngseed = 123)  #1134 taxa and 67 samples

min(sample_sums(ps_obj))  # 482

# plot Shannon
sample_data(ps_cu)$Copper <- as.numeric(as.character(sample_data(ps_cu)$Copper))
plot_richness(ps_cu, x = "Copper", measures = "Shannon") +
  ggtitle("B: 18S - Copper Treatment") +
  geom_point(aes(color = Copper), size = 2) +  
  labs(x = "Copper Concentration (ppb)", y = "Shannon Diversity", color = "Copper") +
  theme(
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_blank(),
    axis.line          = element_line(colour = "black"),
    strip.text         = element_blank(),      
    strip.background   = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 0))+
  coord_cartesian(ylim = c(1.0, 5.0))+
  scale_y_continuous(breaks = c(0, 1.0, 2.0, 3.0, 4.0, 5.0))


plot_richness(ps_ts,x = "Sal_status", measures = "Shannon") +             
  geom_boxplot(aes(fill = Sal_status),width = 0.6, alpha = 0.8, colour = "black") +
  ggtitle("D: 18S - Time Series") +
  labs(x = "Treatment", fill = "Treatment", y = "Shannon Diversity") +       
  scale_fill_manual(values = colors) +                  
  theme(
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_blank(),
    axis.line          = element_line(colour = "black"),
    strip.text.x       = element_blank(),               
    strip.background   = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 0)
  ) +
  scale_fill_manual(
    name = "Treatment",          
    values = colors,                   
    labels = c("15ppt Salinity", "22ppt Salinity", "31ppt Salinity/Prazi"))+
  guides(colour = "none")+
  coord_cartesian(ylim = c(1.0, 4.0))+
  scale_y_continuous(breaks = c(0, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0))


############## Alpha stats ----
### copper
shannon_cu <- estimate_richness(ps_cu, measures = "Shannon")

describe(shannon_cu, fast=TRUE)
#         vars  n mean   sd  min  max range   se
#Shannon     1 22 3.25 0.67 1.71 4.17  2.47 0.14


### time series and salinity
shannon_ts <- estimate_richness(ps_ts, measures = "Shannon")

describe(shannon_ts, fast=TRUE)
#        vars  n mean   sd  min  max range   se
#Shannon    1 45 3.21 0.49 1.55 4.03  2.48 0.07

sample_data(ps_ts) <- cbind(sample_data(ps_ts), sample_data(shannon_ts))
ps_ts_df <- as(sample_data(ps_ts), "data.frame")
ps_ts_df %>%
  group_by(Sal_status) %>%
  summarise(across(Shannon, ~ describe(.x, fast = TRUE)))

### GAMM
library(mgcv)

#### Copper - Shannon diversity
sample_data(ps_cu) <- cbind(sample_data(ps_cu), sample_data(shannon_cu))
ps_cu_df <- as(sample_data(ps_cu), "data.frame")
ps_cu_df$Day <- as.numeric(as.character(ps_cu_df$Day))

gamm_mod_cu <- gamm(Shannon ~ s(Copper, k = 9) + s(Day, k = 10), 
                    data = ps_cu_df,  
                    correlation = corAR1(form = ~ Day), # autocorrelation structure 
                    method = "REML")

summary(gamm_mod_cu$lme)
summary(gamm_mod_cu$gam)

# 1) Parametric terms (e.g., Treatment main effects)
summary(gamm_mod_cu$gam)$p.table   # estimates, SE, t, p for parametric part

# 2) Smooth terms (per-treatment time trends)
summary(gamm_mod_cu$gam)$s.table   # edf, F, p for smooths like s(time_num):Treatment


#### Treatment time series - Shannon diversity
sample_data(ps_ts) <- cbind(sample_data(ps_ts), sample_data(shannon_ts))

ps_ts_df <- as(sample_data(ps_ts), "data.frame")
ps_ts_df$Day <- as.numeric(as.character(ps_ts_df$Day))
ps_ts_df$Sal_status <- factor(ps_ts_df$Sal_status)

ps_ts_df$Sal_status <- relevel(factor(ps_ts_df$Sal_status), ref = "31")


gamm_mod_ts <- gamm(Shannon ~ Sal_status + s(Day, by = Sal_status, k = 10),  # smooth per treatment
                    data = ps_ts_df,
                    correlation = corAR1(form = ~ Day),  # autocorrelation structure
                    method = "REML")

summary(gamm_mod_ts$gam)


# 1) Parametric terms (e.g., Treatment main effects)
summary(gamm_mod_ts$gam)$p.table   # estimates, SE, t, p for parametric part

# 2) Smooth terms (per-treatment time trends)
summary(gamm_mod_ts$gam)$s.table   # edf, F, p for smooths like s(time_num):Treatment


###### BETA DIVERSITY ----
#plot pcoa
pcoa_bray = ordinate(ps_obj, "PCoA", "bray", weighted=TRUE)
plot_ordination(ps_obj, pcoa_bray, color = "Sal_status", shape = "Prazi_status") + geom_point(size = 3) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

pcoa_bray_cu = ordinate(ps_cu, "PCoA", "bray", weighted=TRUE)
sample_data(ps_cu)$Copper <- as.numeric(sample_data(ps_cu)$Copper)
plot_ordination(ps_cu, pcoa_bray_cu, color = "Copper") + geom_point(size = 3) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_gradient(name = "Copper") +
  ggtitle("B: 18S rRNA - Copper Treatment")

pcoa_bray_ts = ordinate(ps_ts, "PCoA", "bray", weighted=TRUE)
plot_ordination(ps_ts, pcoa_bray_ts, color = "Sal_status") + geom_point(size = 3) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("D: 18S rRNA - Time Series")+
  scale_color_manual(
    name = "Treatment",          
    values = colors,                   
    labels = c("15 ppt salinity", "22 ppt salinity", "31 ppt salinity/prazi"))

############## Beta stats ----
# Permanova
bray_dist_cu = phyloseq::distance(ps_cu, method="bray", weighted=F)
pcoa_bray_cu = ordinate(ps_cu, "PCoA", "bray", weighted=TRUE)
sample_data(ps_cu)$Copper <- as.numeric(as.character(sample_data(ps_cu)$Copper))
adonis2(bray_dist_cu ~ sample_data(ps_cu)$Copper)
#                          Df SumOfSqs     R2      F Pr(>F)  
#sample_data(ps_cu)$Copper  1   0.5835 0.08058 1.7529  0.008 **
#Residual                  20   6.6572 0.91942                 
#Total                     21   7.2406 1.00000  

bray_dist_ts = phyloseq::distance(ps_ts, method="bray", weighted=F)
pcoa_bray_ts = ordinate(ps_ts, "PCoA", "bray", weighted=TRUE)
adonis2(bray_dist_ts ~ sample_data(ps_ts)$Sal_status)
#                              Df SumOfSqs      R2     F Pr(>F)    
#sample_data(ps_ts)$Sal_status  2   2.8889 0.21043 5.5969  0.001 ***
#Residual                      42  10.8394 0.78957                  
#Total                         44  13.7283 1.00000     

grouping <- sample_data(ps_ts)$Sal_status
pairwise.adonis(bray_dist_ts, grouping, sim.method = "bray", p.adjust.m = "BH")
#     pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
#1 15 vs 22  1  1.453835 6.365881 0.1538921   0.001      0.001  **
#2 15 vs 31  1  1.535917 5.495406 0.2556549   0.001      0.001  **
#3 22 vs 31  1  1.387484 4.969412 0.1308794   0.001      0.001  **


bray_dist_ts_31 = phyloseq::distance(ps_ts_31, method="bray", weighted=F)
pcoa_bray_ts_31 = ordinate(ps_ts_31, "PCoA", "bray", weighted=TRUE)
sample_data(ps_ts_31)$Copper <- as.numeric(as.character(sample_data(ps_ts_31)$Copper))
adonis2(bray_dist_ts_31 ~ sample_data(ps_ts_31)$Copper)
#                             Df SumOfSqs      R2      F Pr(>F)
#sample_data(ps_ts_31)$Copper  1  0.37647 0.13227 0.9146  0.662
#Residual                      6  2.46971 0.86773              
#Total                         7  2.84617 1.00000 


###### DISPERSION/VARIANCE ----
## gamm for beta dispersion
ps_ts_df <- as(sample_data(ps_ts), "data.frame")
ps_ts_df$Day <- as.numeric(as.character(ps_ts_df$Day))
ps_ts_df$Sal_status <- factor(ps_ts_df$Sal_status)

bray_dist_ts = phyloseq::distance(ps_ts, method="bray", weighted=F)
beta_ts <- betadisper(bray_dist_ts, ps_ts_df$Sal_status, type = "centroid")

disp_df <- data.frame(
  Distance = beta_ts$distances,
  Sal_status = ps_ts_df$Sal_status,
  Day = ps_ts_df$Day
)

disp_df$Sal_status <- relevel(factor(disp_df$Sal_status), ref = "31")

gamm_mod_disp <- gamm(Distance ~ Sal_status + s(Day, by = Sal_status, k = 10),  # smooth per treatment
                      data = disp_df,
                      correlation = corAR1(form = ~ Day),  # autocorrelation structure
                      method = "REML")

summary(gamm_mod_disp$gam)

#Parametric coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)   0.44174    0.81710   0.541    0.592
#Sal_status15 -0.05187    0.81880  -0.063    0.950
#Sal_status22  0.14390    0.81992   0.176    0.862
#
#Approximate significance of smooth terms:
#  edf Ref.df     F p-value  
#s(Day):Sal_status31 1.401  1.401 1.648  0.1708  
#s(Day):Sal_status15 2.460  2.460 0.072  0.9761  
#s(Day):Sal_status22 4.964  4.964 3.235  0.0165 *


# 1) Parametric terms (e.g., Treatment main effects)
summary(gamm_mod_disp$gam)$p.table   # estimates, SE, t, p for parametric part

# 2) Smooth terms (per-treatment time trends)
summary(gamm_mod_disp$gam)$s.table   # edf, F, p for smooths like s(time_num):Treatment


##Dispersion describe 

describe(disp_df, fast=TRUE)
write.csv(disp_df, file = "distance_disp.csv")
disp_df %>%
  group_by(Sal_status) %>%
  summarise(across(Distance, ~ describe(.x, fast = TRUE)))
#Sal_status   Distance$vars    $n $mean    $sd $median  $min  $max $range $skew $kurtosis
#1 31                     1     8 0.593 0.0670   0.578 0.504 0.702  0.198 0.346    -1.45 
#2 15                     1    10 0.399 0.0594   0.374 0.333 0.522  0.189 1.04     -0.466
#3 22                     1    27 0.481 0.0662   0.466 0.377 0.596  0.219 0.358    -1.09 


##### Plot beta dispersion ----
disp_df$Sal_status <- factor(disp_df$Sal_status,
                             levels = c("15", "22", "31"))

Fig_Disp18S <- ggplot(disp_df,aes(x = Sal_status, y = Distance)) + geom_point(alpha = 0.8) +            
  geom_boxplot(aes(fill = Sal_status),width = 0.6, alpha = 0.8, colour = "black") +
  geom_jitter(aes(y=Distance), width = 0.15, alpha = 0.3, shape = 16) +
  ggtitle("B: 18S rRNA - Dispersion") +
  labs(x = "Treatment", fill = "Treatment", y = "Distance from centroid") +       
  theme(
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_blank(),
    axis.line          = element_line(colour = "black"),
    strip.text.x       = element_blank(),               
    strip.background   = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(
    name = "Treatment",          
    values = colors,                   
    labels = c("15 ppt salinity", "22 ppt salinity", "31 ppt salinity/prazi"))+
  guides(colour = "none")

