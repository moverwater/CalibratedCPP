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

# (Reusing your data loading logic)
heled_and_drummond_logs <- list()
cpp_logs <- list()
birthRates <- seq(1.1,5,0.1)

# Loop to read logs
for(i in 1:length(birthRates)){
  heled_and_drummond_logs[[i]] <- readLog(paste0("heled_and_drummond_likelihood_",birthRates[i],".log"),as.mcmc=F,burnin = 0.5)
  cpp_logs[[i]] <- readLog(paste0("cpp_likelihood_",birthRates[i],".log"),as.mcmc=F,burnin = 0.5)
}

# Extract values
heled_and_drummond_lik <- sapply(heled_and_drummond_logs, function(x) x$birthDeath[1])
cpp_lik <- sapply(cpp_logs, function(x) x$birthDeath[1])

# Create Plot 2
# Note: We use aes(color = "Name") to manually map the strings to the legend
plot_lik <- ggplot() + 
  geom_line(aes(x = birthRates, y = heled_and_drummond_lik, color = "Heled & Drummond"), linewidth = 1.1) + 
  geom_point(aes(x = birthRates, y = cpp_lik - lfactorial(100), color = "Calibrated CPP"), size = 2) + 
  
  # ESSENTIAL: Use the SAME name and values as Plot 1
  scale_color_manual(name = "Method", values = custom_colors) +
  
  ylab("Log-likelihoods") + 
  xlab("Birth-rate") + 
  theme_bw() +
  theme(legend.position = "bottom")

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
