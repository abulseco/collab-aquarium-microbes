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
# this has mitochondria, archea, and chloroplasts removed
# also already has the Blanks (6 samples) removed

ps_obj <- readRDS("~/R/GOT_analysis/phy-f-16_beforeRare.rds")  # 2887 taxa and 57 samples

write.csv(sample_sums(ps_obj), "sample_depths_16S.csv")

# time series only
ps_ts = subset_samples(ps_obj, Project == "TimeSeries")  #2887 taxa and 47 samples

# copper treatment only
ps_cu = subset_samples(ps_obj, Project == "Copper_tx")  #2887 taxa and 10 samples

# transform relative abundance
ps_obj <- transform_sample_counts(ps_obj, function(x) x / sum(x))
ps_cu <- transform_sample_counts(ps_cu, function(x) x / sum(x))
ps_ts <- transform_sample_counts(ps_ts, function(x) x / sum(x))

# create ASV tables
otus_ps<- otu_table(ps_obj)
write.csv(otus_ps, file='otu_table_all.csv')

tax_table <- tax_table(ps_obj)
write.csv(tax_table, file = "tax_table_all.csv")

###### Random Forest Time series (salinity) ----

#random forest should be able to pull out the key ASVss that separate these two sites and accurately classify them
# https://rpubs.com/michberr/randomforestmicrobe

library(randomForest)  #version 4.7.1.2

# Make a dataframe of training data with ASVs as column and samples as rows
predictors <- t(otu_table(ps_ts))
dim(predictors)  # 2887   47

# transpose this so # of samples is first, then # of ASVs
predictors <- t(predictors)
dim(predictors)  # 47  2887  

####
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
# randomForest(formula = response_sal ~ ., data = rf_sal.data,      ntree = 100) 
#Type of random forest: classification
#Number of trees: 100
#No. of variables tried at each split: 53
#
#OOB estimate of  error rate: 17.02%
#Confusion matrix:
#   15 22 31 class.error
#15  6  3  1  0.40000000
#22  1 28  0  0.03448276
#31  0  3  5  0.37500000


# What variables are stored in the output?
names(data_sal.classify)

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
  ggtitle("Most important ASVs Salinity 16S rRNA") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(imp_sal.10, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important ASVs Salinity 16S rRNA") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))


# What are those ASVs?
library(knitr)

otunames_sal <- imp_sal.10$predictors
row10_sal <- rownames(tax_table(ps_ts)) %in% otunames_sal
kable(tax_table(ps_ts)[row10_sal, ])

write.csv((tax_table(ps_ts)[row10_sal, ]), "RF_results_Salinity_16S.csv")


###plot random forest results 

rf_sal <- subset_taxa(ps_ts, rownames(tax_table(ps_ts)) %in% c("ASV26","ASV38","ASV41","ASV103","ASV43","ASV34","ASV5","ASV4","ASV21","ASV108"))

plot_bar(rf_sal, "Description", fill = "OTU") + ggtitle("Random Forest ASVs Salinity 16S rRNA") + facet_grid(.~Sal_status, scales = "free") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_fill_brewer(type="qual", palette="PRGn")


#reorder and relabel ASVs
df_sal <- psmelt(rf_sal)

df_sal$OTU <- factor(df_sal$OTU, levels = c("ASV26","ASV38","ASV41","ASV103","ASV43","ASV34","ASV5","ASV4","ASV21","ASV108"),
                     labels = c("ASV26 - Hyphomonadaceae","ASV38 - Woeseiaceae","ASV41 - Pseudomonadaceae","ASV103 - Hyphomonodaceae","ASV43 - Salinisphaeraceae","ASV34 - Cellvibrionaceae","ASV5 - Alteromonadaceae","ASV4 - Rhodobacteraceae","ASV21 - Nitrospiraceae","ASV108 - Marinicellaceae"))

Fig_RF_TS_16S <- ggplot(df_sal, aes(x = Description, y = Abundance, fill = OTU)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  ggtitle("A: 16S rRNA - Time Series") +
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


###### Random Forest All data (DCS) ----
#all data (copper tx and time series)
# How many ASVs do we currently have? 
ntaxa(ps_obj)  # 2887

# Make a dataframe of training data with ASVs as column and samples as rows
predictors_obj <- t(otu_table(ps_obj))
dim(predictors_obj)  # 2887  57

# transpose this so # of samples is first, then # of ASVs
predictors_obj <- t(predictors_obj)
dim(predictors_obj)  #57  2887

####
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
#randomForest(formula = response_DCS ~ ., data = rf_DCS.data, ntree = 100) 
#Type of random forest: classification
#Number of trees: 100
#No. of variables tried at each split: 45
#
#OOB estimate of  error rate: 21.05%
#Confusion matrix:
#   0 1 2 3 class.error
#0 45 0 0 0           0
#1  5 0 0 0           1
#2  2 0 0 0           1
#3  5 0 0 0           1


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
  ggtitle("Most important ASVs DCS 16S rRNA") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(imp_DCS.10, aes(x = predictors_obj, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important ASVs DCS 16S rRNA") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

# What are those ASVs?
library(knitr)

otunames_DCS <- imp_DCS.10$predictors_obj
row10_DCS <- rownames(tax_table(ps_obj)) %in% otunames_DCS
kable(tax_table(ps_obj)[row10_DCS, ])

write.csv((tax_table(ps_obj)[row10_DCS, ]), "RF_results_DCS_16S.csv")


###plot random forest results this way

rf_DCS <- subset_taxa(ps_obj, rownames(tax_table(ps_obj)) %in% c("ASV18","ASV430","ASV1140","ASV674","ASV401","ASV591","ASV502","ASV428","ASV206","ASV468"))

plot_bar(rf_DCS, "Description", fill = "OTU") + ggtitle("Random Forest ASVs DCS 16S rRNA") + facet_grid(.~DCS, scales = "free") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_fill_brewer(type="qual", palette="PRGn")


#reorder and relabel ASVs
df_DCS <- psmelt(rf_DCS)

df_DCS$OTU <- factor(df_DCS$OTU, levels = c("ASV18","ASV430","ASV1140","ASV674","ASV401","ASV591","ASV502","ASV428","ASV206","ASV468"),
                     labels = c("ASV18 - Kordiimonadaceae","ASV430 - Halieaceae","ASV1140 - Saprospiraceae","ASV674 - Rickettsiales family","ASV401 - Thiohalorhabdaceae","ASV591 - Cellvibrionaceae","ASV502 - Myxococcales family","ASV428 - Micavibrionaceae","ASV206 - Beijerinckiaceae","ASV468 - Cellvibrionaceae"))

ggplot(df_DCS, aes(x = Description, y = Abundance, fill = OTU)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  ggtitle("C: 16S rRNA - Disease Condition Status") +
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

# try to make a break in y axis
library(ggbreak)

Fig_RF_DCS_16S <- ggplot(df_DCS, aes(x = Description, y = Abundance, fill = OTU)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  ggtitle("C: 16S rRNA - Disease Condition Status") +
  facet_grid(. ~ DCS, scales = "free", space = "free_x") +
  scale_fill_brewer(type = "qual", palette = "PRGn") +
  xlab("Sample") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_blank(),
    axis.title.x = element_text(vjust = 8, hjust = 0.3, margin = ggplot2::margin(10,0,0,0)), # centers x axis title
    legend.text = element_text(size = 10),
    legend.key.size = (unit(1, "lines")),
    plot.title = element_text(hjust = 0.12, size = 13.2),
    axis.title.y = element_text(angle = 90, vjust = -2.8),
    plot.margin = ggplot2::margin(5, 5, 5, 5)  # Tightens outer spacing
  ) +
  # Y-axis break, no space, free scaling for lower and upper parts
  scale_y_break(
    breaks = c(0.045, 0.29),
    scales = c(1, 0.2),
    space = 0.6,
    ticklabels = c(0, 0.01, 0.02, 0.03, 0.29, 0.30), expand = c(0,0)
  ) +
  # This removes the space under bars
  coord_cartesian(ylim = c(0, NA), expand = FALSE) +
  guides(fill = guide_legend(title = NULL))


###### TAXA PLOTS ----
# custom color palette to use for taxa plots of top 30 families:
custom_col30 = c("#784511","#D2781E","#EAAB6C", "#781122", "#D21E2C","#A51876","#E43FAD","#F098D3","#185EA5","#1E78D2","#3F91E4","#6CABEA","#98C4F0","#117878","#18A5A5","#3FE4E4","#6CEAEA","#98F0F0","#117845","#18A55E","#1ED278","#3FE491","#6CEAAB","#98F0C4","#787811","#A5A518","#D2D21E","#E4E43F","#EAEA6C","#F0F098","#F7F7C5")

# 1. Transform the original phyloseq object to relative abundance
# did this earlier: ps_obj <- transform_sample_counts(ps_obj, function(x) x / sum(x))

# 2. Agglomerate at Family level
ps_family <- tax_glom(ps_obj, taxrank = "Family", NArm = TRUE)

# 3. Get the top 30 families by total abundance (across all samples)
family_sums <- taxa_sums(ps_family)
top30_families <- names(sort(family_sums, decreasing = TRUE))[1:30]

# 4. Prune to only the top 30 families
ps_top30_families <- prune_taxa(top30_families, ps_family)

# 5. Melt for plotting
ps_melt <- psmelt(ps_top30_families)

## make barplot

ps_melt$new_datetime_column <- as.POSIXct(ps_melt$Sample, format = "%m_%d_%y")
ps_melt2 <- ps_melt[order(ps_melt$new_datetime_column), ]
ps_melt2$Sample <- factor(ps_melt2$Sample, levels = unique(ps_melt2$Sample))


plot <- ggplot(ps_melt2, aes(x = Sample, y = Abundance, fill = Family)) +
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
  scale_fill_manual(values=custom_col30)+
  ggtitle("A: 16S rRNA")

ggplotly(plot)

ggsave("Fig_taxa16S.jpg",
       plot = plot,
       width = 10.5, height = 4,
       units = "in",
       dpi = 600)

##
# exporting otu tables to find out actual percentages, etc.
otus_ps_fam<- otu_table(ps_top30_families)
write.csv(otus_ps_fam, file='otu_table_all_fam_rar.csv')

tax_table_fam <- tax_table(ps_top30_families)
write.csv(tax_table_fam, file = "tax_table_all_fam_rar.csv")

metadata_ps_obj <- as(sample_data(ps_obj), "data.frame")
write.csv(metadata_ps_obj, file = "metadata_ps_obj.csv")


###### ALPHA DIVERSITY ----
# this was done in a separate file
# plot Shannon on rarefied data to minimum sample depth = 1184
ps_obj <- readRDS("~/R/GOT_analysis_rar/phy-f-16_FINAL.rds")  # 2082 taxa and 57 samples

ps_cu$Copper <- as.numeric(ps_cu$Copper)

plot_richness(ps_cu, x = "Copper", measures = "Shannon") +
  ggtitle("A: 16S - Copper Treatment") +
  geom_point(aes(color = Copper), size = 2) +  
  labs(x = "Copper Concentration (ppb)", y = "Shannon Diversity", color = "Copper") +
  theme(
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_blank(),
    axis.line          = element_line(colour = "black"),
    strip.text         = element_blank(),      
    strip.background   = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 0)
  ) +
  coord_cartesian(ylim = c(1.0, 5.0))+
  scale_y_continuous(breaks = c(0, 1.0, 2.0, 3.0, 4.0, 5.0))


plot_richness(ps_ts,x = "Sal_status", measures = "Shannon") +             
  geom_boxplot(aes(fill = Sal_status),width = 0.6, alpha = 0.8, colour = "black") +
  ggtitle("C: 16S - Time Series") +
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
  guides(colour = "none") +
  coord_cartesian(ylim = c(1.0, 4.0))+
  scale_y_continuous(breaks = c(0, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0))


############## Alpha stats ----
### copper
shannon_cu <- estimate_richness(ps_cu, measures = "Shannon")

describe(shannon_cu, fast=TRUE)
#          vars  n mean   sd  min  max range   se
#Shannon     1 10  4.3 0.51 3.31 4.94  1.63 0.16


### time series and salinity
shannon_ts <- estimate_richness(ps_ts, measures = "Shannon")

describe(shannon_ts, fast=TRUE)
#        vars  n mean   sd  min  max range   se
#Shannon    1 47 2.63 0.44 1.19 3.45  2.26 0.06

write.csv(shannon_ts, file = "shannon_ts.csv")

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
                    #correlation = corAR1(form = ~ Day), # autocorrelation structure put out an error
                    method = "REML")

# removed the autocorrelation structure above because > acf(residuals(gamm_mod_cu$gam)) showed no significant autocorrelation which was the reason for the error

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
  scale_color_gradient(name = "Copper (ppb)") +
  ggtitle("A: 16S rRNA - Copper Treatment")

pcoa_bray_ts = ordinate(ps_ts, "PCoA", "bray", weighted=TRUE)
colors = c("#1B9E77", "#D95F02","#7570B3")

plot_ordination(ps_ts, pcoa_bray_ts, color = "Sal_status") + geom_point(size = 3) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("C: 16S rRNA - Time Series")+
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
#                           Df SumOfSqs      R2      F Pr(>F)  
#sample_data(ps_cu)$Copper  1  0.54992 0.18558 1.8229  0.027 *
#Residual                   8  2.41335 0.81442                
#Total                      9  2.96327 1.00000    


bray_dist_ts = phyloseq::distance(ps_ts, method="bray", weighted=F)
pcoa_bray_ts = ordinate(ps_ts, "PCoA", "bray", weighted=TRUE)
adonis2(bray_dist_ts ~ sample_data(ps_ts)$Sal_status)

#                              Df SumOfSqs      R2      F Pr(>F)    
#sample_data(ps_ts)$Sal_status  2   2.0690 0.32091 10.396  0.001 ***
#Residual                      44   4.3783 0.67909                  
#Total                         46   6.4472 1.00000 


grouping <- sample_data(ps_ts)$Sal_status
pairwise.adonis(bray_dist_ts, grouping, sim.method = "bray", p.adjust.m = "BH")
#     pairs Df SumsOfSqs   F.Model        R2 p.value p.adjusted sig
#1 15 vs 22  1 0.5098975  8.490698 0.1866469   0.001      0.001  **
#2 15 vs 31  1 1.1033043  6.395845 0.2855818   0.001      0.001  **
#3 22 vs 31  1 1.5539433 14.409373 0.2916324   0.001      0.001  ** 


bray_dist_ts_31 = phyloseq::distance(ps_ts_31, method="bray", weighted=F)
pcoa_bray_ts_31 = ordinate(ps_ts_31, "PCoA", "bray", weighted=TRUE)
sample_data(ps_ts_31)$Copper <- as.numeric(as.character(sample_data(ps_ts_31)$Copper))
adonis2(bray_dist_ts_31 ~ sample_data(ps_ts_31)$Copper)
#                             Df SumOfSqs      R2      F Pr(>F)
#sample_data(ps_ts_31)$Copper  1  0.38241 0.17735 1.2935  0.207
#Residual                      6  1.77387 0.82265              
#Total                         7  2.15628 1.00000 


###### DISPERSION/VARIANCE ----
## gamm for beta dispersion
ps_ts_df$Day <- as.numeric(as.character(ps_ts_df$Day))
ps_ts_df$Sal_status <- factor(ps_ts_df$Sal_status)

bray_dist_ts = phyloseq::distance(ps_ts, method="bray", weighted=F)
beta_ts <- betadisper(bray_dist_ts, ps_ts_df$Sal_status, type = "centroid")

disp_df <- data.frame(
  Distance = beta_ts$distances,
  Sal_status = ps_ts_df$Sal_status,
  Day = ps_ts_df$Day
)

disp_df$Sal_status <- relevel(factor(ps_ts_df$Sal_status), ref = "31")

gamm_mod_disp <- gamm(Distance ~ Sal_status + s(Day, by = Sal_status, k = 10),  # smooth per treatment
                      data = disp_df,
                      correlation = corAR1(form = ~ Day),  # autocorrelation structure
                      method = "REML")

summary(gamm_mod_disp$gam)

#Parametric coefficients:
#            Estimate Std. Error t value Pr(>|t|)
#(Intercept)   0.32456    2.21787   0.146    0.884
#Sal_status15 -0.07185    2.21830  -0.032    0.974
#Sal_status22 -0.04752    2.21803  -0.021    0.983

#Approximate significance of smooth terms:
#                      edf Ref.df      F  p-value    
#s(Day):Sal_status31 1.832  1.832 16.667 1.95e-05 ***
#s(Day):Sal_status15 1.934  1.934  0.921    0.364    
#s(Day):Sal_status22 2.836  2.836  2.078    0.116 


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
#  Sal_status Distance$vars    $n $mean    $sd $median  $min  $max $range $skew $kurtosis    $se
#1 15                     1    10 0.235 0.0750   0.203 0.147 0.400  0.252 0.898    -0.367 0.0237
#2 22                     1    29 0.229 0.0574   0.212 0.162 0.401  0.239 1.25      1.18  0.0106
#3 31                     1     8 0.503 0.137    0.446 0.369 0.714  0.345 0.582    -1.55  0.0486 


############## Plot beta dispersion ----
disp_df$Sal_status <- factor(disp_df$Sal_status,
                             levels = c("15", "22", "31"))

Fig_Disp16S <- ggplot(disp_df,aes(x = Sal_status, y = Distance)) + geom_point(alpha = 0.8) +            
  geom_boxplot(aes(fill = Sal_status),width = 0.6, alpha = 0.8, colour = "black") +
  geom_jitter(aes(y=Distance), width = 0.15, alpha = 0.3, shape = 16) +
  ggtitle("A: 16S rRNA - Dispersion") +
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
