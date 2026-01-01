# GOT Manuscript
# FAPROTAX analysis by experiment
# Figure 7 - Boxplots, exploratory

# This script assumes you have run "Fig7_Faprotax_byExperiment.R"
# To produce FAPROTAX output from MicroEco 
# And to visualize stats

# Load libraries
library(mgcv)
library(ggplot2)

# Output from MicroEco
faprotax_boxplot <- read.csv("FAPROTAX/function-boxplot-timeseries.csv", header = T)

# Boxplots----
## Time series----
### Chemoheterotrophy----
# Ensure Sal_status is a factor
faprotax_boxplot$Sal_status <- as.factor(faprotax_boxplot$Sal_status)
faprotax_boxplot$Treatment <- as.factor(faprotax_boxplot$Treatment)

chemohet_boxplot <- ggplot(faprotax_boxplot, aes(x = Treatment, y = chemoheterotrophy)) +
  geom_point(alpha = 0.8) +
  geom_boxplot(aes(fill = Sal_status), width = 0.6, alpha = 0.8) +
  scale_fill_manual(values=c("#1B9E77", "#D95F02","#7570B3"))+
  labs(x = "Treatment", y = "Chemoheterotrophy (%)") +
  pretty.theme.noangle() +
  theme(legend.position = "none")
chemohet_boxplot

### Fermentation----
ferment_boxplot <- ggplot(faprotax_boxplot, aes(x = Treatment, y = fermentation)) +
  geom_point(alpha = 0.8) +
  geom_boxplot(aes(fill = Sal_status), width = 0.6, alpha = 0.8) +
  scale_fill_manual(values=c("#1B9E77", "#D95F02","#7570B3"))+
  labs(x = "Treatment", y = "Fermentation (%)") +
  pretty.theme.noangle() +
  theme(legend.position = "none", 
        x.axis.title = "none", 
        x.axis.text = "none")
ferment_boxplot

### Nitrate oxidation----
nitrite_ox_boxplot <- ggplot(faprotax_boxplot, aes(x = Treatment, y = aerobic_nitrite_oxidation)) +
  geom_point(alpha = 0.8) +
  geom_boxplot(aes(fill = Sal_status), width = 0.6, alpha = 0.8) +
  scale_fill_manual(values=c("#1B9E77", "#D95F02","#7570B3"))+
  labs(x = "Treatment", y = "Nitrite Oxidation (%)") +
  pretty.theme.noangle() +
  theme(legend.position = "none")
nitrite_ox_boxplot

### Nitrification----
nitrification_boxplot <- ggplot(faprotax_boxplot, aes(x = Treatment, y = nitrification)) +
  geom_point(alpha = 0.8) +
  geom_boxplot(aes(fill = Sal_status), width = 0.6, alpha = 0.8) +
  scale_fill_manual(values=c("#1B9E77", "#D95F02","#7570B3"))+
  labs(x = "Treatment", y = "Chemoheterotrophy (%)") +
  pretty.theme.noangle() +
  theme(legend.position = "none")
nitrification_boxplot

colors <- c("#1B9E77", "#D95F02","#7570B3")
# With kerry formatting----

## Chemoheterotrophy----
chemo_boxplot <- ggplot(faprotax_boxplot, aes(x = Sal_status, y = chemoheterotrophy)) +
  geom_point(alpha = 0.8) +
  geom_boxplot(aes(fill = Sal_status),width = 0.6, alpha = 0.8, colour = "black") +  
  ggtitle("A: Chemoheterotrophy") +
  labs(x = "Treatment", fill = "Treatment", y = "Chemoheterotrophy (%)") +       
  theme(
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_blank(),
    axis.line          = element_line(colour = "black"),
    strip.text.x       = element_blank(),               
    strip.background   = element_blank(),
    axis.text.x        = element_text(size=10, color = "black", angle = 45, hjust = 1)
  ) +
  scale_fill_manual(
    name = "Treatment",          
    values = colors,                   
    labels = c("15ppt Salinity", "22ppt Salinity", "31ppt Salinity/Prazi"))+
  guides(colour = "none")
chemo_boxplot

## Fermentation----
fermentation_boxplot <- ggplot(faprotax_boxplot, aes(x = Sal_status, y = fermentation)) +
  geom_point(alpha = 0.8) +
  geom_boxplot(aes(fill = Sal_status),width = 0.6, alpha = 0.8, colour = "black") +  
  ggtitle("B: Fermentation") +
  labs(x = "Treatment", fill = "Treatment", y = "Fermentation (%)") +       
  theme(
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_blank(),
    axis.line          = element_line(colour = "black"),
    strip.text.x       = element_blank(),               
    strip.background   = element_blank(),
    axis.text.x        = element_text(size=10, color = "black", angle = 45, hjust = 1)
  ) +
  scale_fill_manual(
    name = "Treatment",          
    values = colors,                   
    labels = c("15ppt Salinity", "22ppt Salinity", "31ppt Salinity/Prazi"))+
  guides(colour = "none")
fermentation_boxplot

## Nitrate oxidation----
nitrite_oxidation_boxplot <- ggplot(faprotax_boxplot, aes(x = Sal_status, y = aerobic_nitrite_oxidation)) +
  geom_point(alpha = 0.8) +
  geom_boxplot(aes(fill = Sal_status),width = 0.6, alpha = 0.8, colour = "black") +  
  ggtitle("C: Nitrite Oxidation") +
  labs(x = "Treatment", fill = "Treatment", y = "Nitrite Oxidation (%)") +       
  theme(
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_blank(),
    axis.line          = element_line(colour = "black"),
    strip.text.x       = element_blank(),               
    strip.background   = element_blank(),
    axis.text.x        = element_text(size=10, color = "black", angle = 45, hjust = 1)
  ) +
  scale_fill_manual(
    name = "Treatment",          
    values = colors,                   
    labels = c("15ppt Salinity", "22ppt Salinity", "31ppt Salinity/Prazi"))+
  guides(colour = "none")
nitrite_oxidation_boxplot

## Nitrification----
nitrification_boxplot <- ggplot(faprotax_boxplot, aes(x = Sal_status, y = nitrification)) +
  geom_point(alpha = 0.8) +
  geom_boxplot(aes(fill = Sal_status),width = 0.6, alpha = 0.8, colour = "black") +  
  ggtitle("D: Nitrification") +
  labs(x = "Treatment", fill = "Treatment", y = "Nitrification") +       
  theme(
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_blank(),
    axis.line          = element_line(colour = "black"),
    strip.text.x       = element_blank(),               
    strip.background   = element_blank(),
    axis.text.x        = element_text(size=10, color = "black", angle = 45, hjust = 1)
  ) +
  scale_fill_manual(
    name = "Treatment",          
    values = colors,                   
    labels = c("15ppt Salinity", "22ppt Salinity", "31ppt Salinity/Prazi"))+
  guides(colour = "none")
nitrification_boxplot

## combine----
all_boxplots <- (chemo_boxplot + 
                   fermentation_boxplot + 
                   nitrite_oxidation_boxplot + 
                   nitrification_boxplot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right",
        legend.justification = "left")

all_boxplots

# Parasites
parasites_boxplot <- ggplot(faprotax_boxplot, aes(x = Treatment, y = animal_parasites_or_symbionts)) +
  geom_point(alpha = 0.8) +
  geom_boxplot(aes(fill = Sal_status), width = 0.6, alpha = 0.8) +
  scale_fill_manual(values=c("#1B9E77", "#D95F02","#7570B3"))+
  labs(x = "Treatment", y = "Parasites (%)") +
  theme(
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_blank(),
    axis.line          = element_line(colour = "black"),
    strip.text.x       = element_blank(),               
    strip.background   = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 0)
  ) +
  theme(legend.position = "none")
parasites_boxplot

# nitrate_red_boxplot <- ggplot(faprotax_boxplot, aes(x = Treatment, y = nitrate_reduction)) +
#   geom_boxplot(aes(fill = Sal_status)) +
#   scale_fill_manual(values=c("#1B9E77", "#D95F02","#7570B3"))+
#   labs(x = "Treatment", y = "Nitrate Reduction (%)") +
#   pretty.theme.noangle() +
#   theme(legend.position = "none")
# nitrate_red_boxplot

functional_plots_timeseries <- chemohet_boxplot + ferment_boxplot + nitrite_ox_boxplot + nitrification_boxplot
functional_plots_timeseries

# With N-Cycling genes


# Stats----
library(lubridate)
faprotax_boxplot$True_date <- mdy(faprotax_boxplot$True_date)
faprotax_boxplot$time_num <- as.numeric(faprotax_boxplot$True_date - min(faprotax_boxplot$True_date, na.rm = TRUE))
# 
# # Test that it worked
# is.data.frame(faprotax_boxplot)      # should be TRUE
# "True_date" %in% names(faprotax_boxplot)
# class(faprotax_boxplot$True_date)    # should be "Date"
# summary(faprotax_boxplot$time_num)   # should be numeric

# # Dummy variable
# faprotax_boxplot$one <- factor(1)
# 
# m_chemohet <- gamm(chemoheterotrophy ~ Treatment + s(time_num, by = Treatment, k = 6),
#   data = faprotax_boxplot,
#   correlation = corCAR1(form = ~ time_num))  # continuous-time AR(1))

# GAMM
# Fit GAMM

# Double check the original metadata for the "days after start" 
library(mgcv)
faprotax_boxplot$Treatment <- as.factor(faprotax_boxplot$Treatment)

gamm_mod <- gamm(nitrification ~ Treatment + s(time_num, by = Treatment, k = 10),  # smooth per treatment
  data = faprotax_boxplot,
  correlation = corAR1(form = ~ time_num)  # autocorrelation structure
)

summary(gamm_mod$gam)
summary(gamm_mod$lme)

library(emmeans)
emm <- emmeans(gamm_mod$gam, ~ Treatment,
               data = faprotax_boxplot,
               weights = "proportional")   # averages over observed covariates
pairs(emm, adjust = "tukey")               # posthoc p-values
confint(emm)

# 1) Parametric terms (e.g., Treatment main effects)
ptab <- summary(gamm_mod$gam)$p.table   # estimates, SE, t, p for parametric part
ptab

# 2) Smooth terms (per-treatment time trends)
stab <- summary(gamm_mod$gam)$s.table   # edf, F, p for smooths like s(time_num):Treatment
stab

# 3) Overall tests via anova.gam (approximate)
anova(gamm_mod$gam)  # parametric vs smooth blocks with F and p

## Kerry's code----
### chemoheterotrophy----
# We need to re-order the factors so that 31 is the reference instead 
# of the default of 15:
faprotax_boxplot$Sal_status <- relevel(factor(faprotax_boxplot$Sal_status), ref = "31")

# faprotax_boxplot$Sal_status <- as.factor(faprotax_boxplot$Sal_status)
# faprotax_boxplot$Treatment <- as.factor(faprotax_boxplot$Treatment)
# faprotax_boxplot$time_num <- as.numeric(faprotax_boxplot$Date)

gamm_fapro_chemohet <- gamm(chemoheterotrophy ~ Sal_status + s(Day, by = Sal_status, k = 10),
                            data = faprotax_boxplot,
                            correlation = corAR1(form = ~ Day),
                            method = "REML",
                            na.action = na.exclude)
summary(gamm_fapro_chemohet$gam)

gamm_fapro_nit <- gamm(nitrification ~ Sal_status + s(Day, by = Sal_status, k = 10),
                            data = faprotax_boxplot,
                            correlation = corAR1(form = ~ Day),
                            method = "REML",
                            na.action = na.exclude)
summary(gamm_fapro_nit$gam)

gamm_fapro_fer <- gamm(fermentation ~ Sal_status + s(Day, by = Sal_status, k = 10),
                       data = faprotax_boxplot,
                       correlation = corAR1(form = ~ Day),
                       method = "REML",
                       na.action = na.exclude)
summary(gamm_fapro_fer$gam)

gamm_fapro_nitrox <- gamm(aerobic_nitrite_oxidation ~ Sal_status + s(Day, by = Sal_status, k = 10),
                       data = faprotax_boxplot,
                       correlation = corAR1(form = ~ Day),
                       method = "REML",
                       na.action = na.exclude)
summary(gamm_fapro_nitrox$gam)


library(gratia)

# Compare 31 vs 15 ppt smooths
diff_31_15 <- difference_smooths(
  gamm_fapro_chemohet$gam,
  select = "s(Day)",
  conditions = data.frame(Sal_status = c("15", "31"))
)

draw(diff_31_15)

### nitrification----
faprotax_boxplot$Sal_status <- as.factor(faprotax_boxplot$Sal_status)
faprotax_boxplot$Treatment <- as.factor(faprotax_boxplot$Treatment)
faprotax_boxplot$time_num <- as.numeric(faprotax_boxplot$Date)

gamm_fapro_nit <- gamm(nitrification ~ Sal_status + s(Day, by = Sal_status, k = 10),
                            data = faprotax_boxplot,
                            correlation = corAR1(form = ~ Day),
                            method = "REML",
                            na.action = na.exclude)
summary(gamm_fapro_nit$gam)

library(gratia)

# Compare 31 vs 15 ppt smooths
diff_31_15 <- difference_smooths(
  gamm_fapro_nit$gam,
  select = "s(Day)",
  conditions = data.frame(Sal_status = c("15", "31"))
)

draw(diff_31_15)

# emm_par <- emmeans(gamm_fapro_chemohet$lme, ~ Sal_status, data = faprotax_boxplot)
# pairs(emm_par)
# 
# library(gratia)
# draw(gamm_fapro_chemohet$gam, select = 1:3)

##### pairwise GAMM -----
# not sure this is correct, but will use in the meantime. It is a manual approach to getting pairwise comparisons from GAMM model predictions

# get predicted means per group

newdata <- data.frame(
  Sal_status = factor(c(15, 22, 31), levels = c(15, 22, 31)),
  Day = mean(faprotax_boxplot$Day)
)

pred <- predict(gamm_fapro_chemohet$gam, newdata = newdata, se.fit = TRUE) # if there were no smooth terms, you would use gamm_mod_ts$lme, but we have smooth terms
pred_df <- data.frame(newdata, fit = pred$fit, se = pred$se.fit)
pred_df

#Compute pairwise comparisons manually
#Compute pairwise contrasts (differences between predicted means) and test whether they are significant using standard errors and t-tests.

pairs <- combn(1:nrow(pred_df), 2, simplify = FALSE)

pairwise_results <- do.call(rbind, lapply(pairs, function(idx) {
  i <- idx[1]
  j <- idx[2]
  diff <- pred_df$fit[i] - pred_df$fit[j]
  se_diff <- sqrt(pred_df$se[i]^2 + pred_df$se[j]^2)
  t_val <- diff / se_diff
  p_val <- 2 * pt(-abs(t_val), df = nrow(faprotax_boxplot) - length(unique(faprotax_boxplot$Sal_status)))
  data.frame(
    Comparison = paste(pred_df$Sal_status[i], "vs", pred_df$Sal_status[j]),
    Estimate = diff,
    SE = se_diff,
    t_value = t_val,
    p_value = p_val
  )
}))

pairwise_results


# adjust p-values for multiple comparisons
pairwise_results$p_adj <- p.adjust(pairwise_results$p_value, method = "BH")
pairwise_results

## Copper----
# Need to make sure you are using the correct dataset here
# spe func results --> is this really what you want? double check! 
faprotax_copper <- read.csv("FAPROTAX/function-scatter-copper.csv", header = TRUE)

# Change dates
# Change date to the proper format
# For some reason this is changing the format of the data incorrectly?? 
# Need to make sure "True_date" has full length year (2016 not 16)
faprotax_copper $True_date <- as.Date(faprotax_copper $True_date, format = "%m/%d/%Y")
faprotax_copper  <- faprotax_copper [order(faprotax_copper $True_date), ]

### chemohet----
# chemohet_scatter <- ggplot(faprotax_copper, aes(x = True_date, y = chemoheterotrophy)) +
#   geom_point() +
#   geom_line() +
#   scale_x_date(date_labels = "%b %Y") +
#   pretty.theme() +
#   labs(x = "Date", y = "Chemoheterotrophy (%)") 
# # +
# #   theme(axis.title.x = element_blank(),     # Removes x-axis title
# #         axis.text.x = element_blank())
# chemohet_scatter

chemohet_scatter <- ggplot(faprotax_copper, aes(x = Copper, y = chemoheterotrophy)) +
  geom_point() +
  geom_line() +
  pretty.theme() +
  labs(x = "Copper (ppb)", y = "Chemoheterotrophy (%)") 
chemohet_scatter

ferment_scatter <- ggplot(faprotax_copper, aes(x = Copper, y = fermentation)) +
  geom_point() +
  geom_line() +
  pretty.theme() +
  labs(x = "Copper (ppb)", y = "Fermentation (%)") 
ferment_scatter

nitrite_ox_scatter <- ggplot(faprotax_copper, aes(x = Copper, y = aerobic_nitrite_oxidation)) +
  geom_point() +
  geom_line() +
  pretty.theme() +
  labs(x = "Copper (ppb)", y = "Nitrate Oxidation (%)") 
nitrite_ox_scatter

nitrification_scatter <- ggplot(faprotax_copper, aes(x = Copper, y = nitrification)) +
  geom_point() +
  geom_line() +
  pretty.theme() +
  labs(x = "Copper (ppb)", y = "Nitrification (%)") 
nitrification_scatter

# Combine the plots 
# fapro_combined <- (chemohet_boxplot | chemohet_scatter) / (ferment_boxplot | ferment_scatter) /
#   (nitrite_ox_boxplot | nitrite_ox_scatter) / (nitrification_boxplot | nitrification_scatter)
# fapro_combined

fapro_time_combined <- (chemohet_boxplot | ferment_boxplot) / (nitrite_ox_boxplot | nitrification_boxplot) 
fapro_time_combined
