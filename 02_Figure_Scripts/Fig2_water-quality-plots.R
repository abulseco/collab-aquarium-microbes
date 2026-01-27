# GOT Manuscript
# Figure 2

# Setup----
## Load libraries----
library(ggplot2); library(patchwork)

## Figure formatting----
pretty.theme <- function(){
  theme_bw()+
    theme(axis.text.x=element_text(size=14, color = "black", angle = 45, hjust = 1),
          axis.text.y=element_text(size=14, color = "black"),
          axis.title.x=element_text(size=14, color = "black"),             
          axis.title.y=element_text(size=14, color = "black"),             
          panel.grid.major.x=element_blank(),                                          
          panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank(),  
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
          plot.title = element_text(size=20),
          legend.text = element_text(size=12, face="italic"),          
          legend.title = element_blank(),                              
          legend.position="none")
}

pretty.theme.noaxes <- function(){
  theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),             
          panel.grid.major.x=element_blank(),                                          
          panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank(),  
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
          legend.title = element_blank(),                              
          legend.position="none")
}

## Import data----
# Using the 16S map for water quality data
water_qual <- read.csv("GOT_16S_map_forWaterQuality.csv", header = T)

# Plots----
## Salinity----
# Filter out NAs
water_qual <- water_qual[!is.na(water_qual$Salinity) & !is.na(water_qual$True_date), ]

# Convert to a continuous variable
water_qual$Salinity <- as.numeric(as.character(water_qual$Salinity))
water_qual$Copper <- as.numeric(as.character(water_qual$Copper))
water_qual$NO3 <- as.numeric(as.character(water_qual$NO3))
water_qual$NO2 <- as.numeric(as.character(water_qual$NO2))
water_qual$NH3 <- as.numeric(as.character(water_qual$NH3))
str(water_qual)

# Change date to the proper format
# For some reason this is changing the format of the data incorrectly?? 
# Need to make sure "True_date" has full length year (2016 not 16)
water_qual$True_date <- as.Date(water_qual$True_date, format = "%m/%d/%Y")
water_qual <- water_qual[order(water_qual$True_date), ]

# ggplot(water_qual, aes(x = True_date, y = Salinity)) +
#   geom_line() +
#   geom_point() +
#   labs(x = "Date", y = "Salinity", title = "Salinity Over Time") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# df$Date <- as.Date(df$Date, format = "%m/%d/%Y")
 
# Example shaded date ranges
# DCS = 1
shaded_regions1 <- data.frame(
  xmin = as.Date(c("2016-03-06", "2016-05-01", "2016-09-13", "2017-10-27")),
  xmax = as.Date(c("2016-03-06", "2016-05-08","2016-09-28", "2017-10-27")),
  ymin = -Inf,
  ymax = Inf
)

# DCS = 2
shaded_regions2 <- data.frame(
  xmin = as.Date(c("2016-10-26", "2017-09-13")),
  xmax = as.Date(c("2016-10-28", "2017-09-15")),
  ymin = -Inf,
  ymax = Inf
)

# DCS = 3
shaded_regions3 <- data.frame(
  xmin = as.Date(c("2016-02-17", "2016-08-30","2017-09-27")),
  xmax = as.Date(c("2016-02-28", "2016-08-31","2017-10-11")),
  ymin = -Inf,
  ymax = Inf
)

# Prazi administered
prazi_dates <- as.Date(c("2016-08-30", "2016-09-13", "2016-09-28", "2016-10-12", "2016-10-26", "2016-11-09", "2016-11-22", "2016-12-07"))
prazi_dates

## Salinity----
salinity_plot <- ggplot(water_qual, aes(x = True_date, y = Salinity)) +
  geom_vline(xintercept = as.numeric(prazi_dates), linetype = "dashed", color = "gray40", linewidth = 0.3, alpha = 0.4) +
  geom_rect(data = shaded_regions1,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "#ffffcc", alpha = 0.8) +
  geom_rect(data = shaded_regions2,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "#ffcc66", alpha = 0.8) +
  geom_rect(data = shaded_regions3,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "#ff9933", alpha = 0.8) +
  geom_point() +
  geom_line() +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  pretty.theme() +
  labs(x = "Date", y = "Salinity (ppt)") +
  theme(axis.title.x = element_blank(),     # Removes x-axis title
        axis.text.x = element_blank())
salinity_plot 

# For a record of prior color regions
# geom_rect(data = shaded_regions1,
#           aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
#           inherit.aes = FALSE,
#           fill = "#21C25C", alpha = 0.5) +
#   geom_rect(data = shaded_regions2,
#             aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
#             inherit.aes = FALSE,
#             fill = "#EEBF61", alpha = 0.5) +
#   geom_rect(data = shaded_regions3,
#             aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
#             inherit.aes = FALSE,
#             fill = "#F59B4D", alpha = 0.5) +

## Copper----
cop_plot <- ggplot(water_qual, aes(x = True_date, y = Copper)) +
  geom_rect(data = shaded_regions1,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "#ffffcc", alpha = 0.5) +
  geom_rect(data = shaded_regions2,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "#ffcc66", alpha = 0.5) +
  geom_rect(data = shaded_regions3,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "#ff9933", alpha = 0.5) +
  geom_point() +
  geom_line() +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  pretty.theme() +
  labs(x = "Date", y = "Copper (ppb)") +
  theme(axis.title.x = element_blank(),     # Removes x-axis title
        axis.text.x = element_blank())
cop_plot 

## Nitrate----
nitrate_plot <- ggplot(water_qual, aes(x = True_date, y = NO3)) +
  geom_rect(data = shaded_regions1,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "#ffffcc", alpha = 0.5) +
  geom_rect(data = shaded_regions2,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "#ffcc66", alpha = 0.5) +
  geom_rect(data = shaded_regions3,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "#ff9933", alpha = 0.5) +
  geom_point() +
  geom_line() +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  pretty.theme() +
  labs(x = "Date", y = "Nitrate (ppm)") +
  theme(axis.title.x = element_blank(),     # Removes x-axis title
        axis.text.x = element_blank())
nitrate_plot 

### Boxplot----
nitrate_boxplot <- ggplot(water_qual, aes(x = True_date, y = NO3)) +
  geom_boxplot(aes(x = Sal_status, y = NO3)) +
  pretty.theme() +
  labs(x = "Date", y = "Nitrate (ppm)") 
nitrate_boxplot # This one is interesting

## Nitrite----
nitrite_plot <- ggplot(water_qual, aes(x = True_date, y = NO2)) +
  geom_rect(data = shaded_regions1,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "#ffffcc", alpha = 0.5) +
  geom_rect(data = shaded_regions2,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "#ffcc66", alpha = 0.5) +
  geom_rect(data = shaded_regions3,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "#ff9933", alpha = 0.5) +
  geom_point() +
  geom_line() +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  pretty.theme() +
  labs(x = "Date", y = "Nitrite (ppm)") +
  theme(axis.title.x = element_blank(),     # Removes x-axis title
        axis.text.x = element_blank())
nitrite_plot 

### Boxplot----
nitrite_boxplot <- ggplot(water_qual, aes(x = True_date, y = NO3)) +
  geom_boxplot(aes(x = Sal_status, y = NO2)) +
  pretty.theme() +
  labs(x = "Date", y = "Nitrite (ppm)") +
  theme(axis.title.x = element_blank(),     # Removes x-axis title
        axis.text.x = element_blank())
nitrite_boxplot 

## Ammonia----
amm_plot <- ggplot(water_qual, aes(x = True_date, y = NH3)) +
  geom_rect(data = shaded_regions1,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "#ffffcc", alpha = 0.5) +
  geom_rect(data = shaded_regions2,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "#ffcc66", alpha = 0.5) +
  geom_rect(data = shaded_regions3,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "#ff9933", alpha = 0.5) +
  geom_point() +
  geom_line() +
  scale_x_date(date_breaks = "3 month", date_labels = "%b %Y") +
  pretty.theme() +
  labs(x = "Date", y = "Ammonia (ppm)")
amm_plot 

## Feed----
# Removing feed from the final plot
feed_plot <- ggplot(water_qual, aes(x = True_date, y = Food_total)) +
  geom_rect(data = shaded_regions1,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "#21C25C", alpha = 0.5) +
  geom_rect(data = shaded_regions2,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "#EEBF61", alpha = 0.5) +
  geom_rect(data = shaded_regions3,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "#F59B4D", alpha = 0.5) +
  geom_point() +
  geom_line() +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  pretty.theme() +
  labs(x = "Date", y = "Food Total (lbs)") 

+
  theme(axis.title.x = element_blank(),     # Removes x-axis title
        axis.text.x = element_blank())
feed_plot

## Crypto plot----
# These data are not rarefied, need to decide what to cut off at if you take this route
## Salinity----
crypto_plot <- ggplot(water_qual, aes(x = True_date, y = Crypto)) +
  geom_rect(data = shaded_regions1,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "#ffffcc", alpha = 0.5) +
  geom_rect(data = shaded_regions2,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "#ffcc66", alpha = 0.5) +
  geom_rect(data = shaded_regions3,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "#ff9933", alpha = 0.5) +
  geom_point() +
  geom_line() +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
  pretty.theme() +
  labs(x = "Date", y = expression(italic("C. irritans") ~ "(%)")) 
crypto_plot

# water_qual_plot <- salinity_plot / cop_plot / nitrate_plot / crypto_plot
# water_qual_plot

# GROUP PLOTS----
# Group plots together differently
nitrogen_plots <- nitrate_plot / nitrite_plot / amm_plot
nitrogen_plots

water_qual_plots <- salinity_plot / cop_plot / crypto_plot
water_qual_plots

combined_plots <- water_qual_plots | nitrogen_plots
combined_plots


