library(readr)
library(dplyr)
library(ggplot2)

# Read CellProfiler output
df <- read_csv("sample_cellprofiler_output.csv")
sec_per_frame <- 3 # seconds between frames

# Calculate per-cell mean (and SEM) fluorescence across frames
summary_df <- df %>%
  group_by(ImageNumber) %>%
  mutate(AVG_FL = mean(Intensity_IntegratedIntensity),
         SE_FL = sd(Intensity_IntegratedIntensity) / sqrt(n()),
         Time = (ImageNumber - 1)*sec_per_frame/60) %>%
  dplyr::select(ImageNumber, AVG_FL, SE_FL, Time) %>%
  distinct()

# Normalize fluorescence measurements to initial baseline and plot
norm_factor <- as.double(summary_df['AVG_FL'][1,])

summary_df %>%
  mutate(NORM_FL = AVG_FL / norm_factor,
         SE_NORM = SE_FL / norm_factor) %>%
  ggplot() + 
  geom_ribbon(aes(x = Time, ymin = NORM_FL - SE_NORM, ymax = NORM_FL + SE_NORM), alpha = 0.2) +
  geom_line(aes(x = Time, y = NORM_FL), size = 1) +
  theme_linedraw()
