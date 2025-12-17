library(ggplot2)
library(readr)
library(dplyr)
library(ggpubr) # Required for ggarrange

# Define consistent colors and names for both plots
# This ensures the legend merges correctly
custom_colors <- c("Calibrated CPP" = "black", "Heled & Drummond" = "red")

# ==========================================
# PART 1: Benchmark Plot (Time)
# ==========================================

# 1. Read and Clean Data
# (Assuming benchmark_results.csv exists in the WD)
setwd("~/code/CalibratedCPP/calibratedcpp-beast/validation/calibratedcpp/likelihood_benchmark/")
data <- read_csv("benchmark_results.csv")

clean_data <- data %>%
  select(Benchmark, `Param: nCalibrations`, Score, Unit) %>%
  rename(
    Method = Benchmark,
    Calibrations = `Param: nCalibrations`,
    Time = Score
  ) %>%
  mutate(
    # Standardize names to match the Color Palette keys exactly
    Method = ifelse(grepl("measureCPP", Method), "Calibrated CPP", "Heled & Drummond"),
    Calibrations = as.numeric(Calibrations)
  )

# 2. Create Plot 1
plot_time <- ggplot(clean_data, aes(x = Calibrations, y = Time, color = Method, group = Method)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  
  # Log scale for Time
  scale_y_log10(labels = scales::comma) +
  
  # Force specific colors and legend title
  scale_color_manual(name = "Method", values = custom_colors) +
  
  labs(
    x = "Number of non-nested calibrations (k)",
    y = expression(paste("Time (", mu, "s)"))
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 12)
  )

# ==========================================
# PART 2: Likelihood Plot (Validation)
# ==========================================

setwd("~/code/CalibratedCPP/calibratedcpp-beast/validation/calibratedcpp/likelihood_comparison/")

lik_vals <- read.csv("comparison_results.csv")

plot_lik <- ggplot() + 
  # Move color inside aes() to create a mapping for the legend
  geom_line(aes(x = lik_vals$BirthRate, y = lik_vals$HeledAndDrummond_LogLikelihood, color = "Heled and Drummond"), linewidth = 1.1) + 
  geom_point(aes(x = lik_vals$BirthRate, y = lik_vals$CPP_LogLikelihood, color = "CPP"), size = 1.5) + 
  
  # Manually define the colors for the labels created above
  scale_color_manual(name = "Method", 
                     values = c("Heled and Drummond" = "red", "CPP" = "black")) +
  
  ylab("Log-likelihoods") + 
  xlab("Birth-rate") + 
  theme_bw() +
  theme(legend.position = "bottom") # Optional: moves legend to bottom

# ==========================================
# PART 3: Combine and Save
# ==========================================

# ggarrange will now detect that the legends are identical and merge them
combined_plot <- ggarrange(plot_time, plot_lik, 
                           ncol = 2, nrow = 1, 
                           common.legend = TRUE, 
                           legend = "bottom", labels = "AUTO")

print(combined_plot)

ggsave("combined_benchmark_validation.pdf", 
       combined_plot, device = "pdf", width = 9, height = 4)
