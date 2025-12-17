library(ggplot2)

setwd("~/code/CalibratedCPP/calibratedcpp-beast/validation/calibratedcpp/likelihood_comparison/")

lik_vals <- read.csv("comparison_results.csv")

liks <- ggplot() + 
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

ggsave("likelihood_comparison.pdf", liks, device="pdf",width=6,height=4)

