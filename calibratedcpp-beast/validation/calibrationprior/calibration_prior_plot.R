setwd("~/code/CalibratedCPP/calibratedcpp-beast/validation/calibrationprior/")

log <- readLog("test_calibration_prior.log", as.mcmc = F)

colnames(log) <- c("Sample", "tree.height", "tree.treeLength", 
                   "clade.10", "clade.9", "clade.8", "clade.7", 
                   "clade.5", "clade.4", "clade.3", "clade.6", 
                   "clade.2", "clade.1")

mcmc_cov <- coverage <- sapply(1:N, function(i) mean(log[[paste0("clade.",i)]] >= t_lo[i] & log[[paste0("clade.",i)]] <= t_hi[i]))

plots <- list()
for (i in 1:N) {
  clade_name <- paste0("clade.", i)
  
  # separate data frames for different-length series
  df_log  <- data.frame(value = log[[clade_name]])
  df_wsim <- data.frame(value = Wsim[, i])
  
  plots[[i]] <- ggplot() +
    geom_density(data = df_log, aes(x = value), fill = "skyblue", alpha = 0.5) +
    geom_density(data = df_wsim, aes(x = value), fill = "grey", alpha = 0.5) +
    geom_vline(xintercept = t_lo[i], color = "red", linetype = "dashed") +
    geom_vline(xintercept = t_hi[i], color = "red", linetype = "dashed") +
    labs(
      title = paste0("Clade ", i, " (Coverage=", round(mcmc_cov[i], 3), ")"),
      x = TeX(paste0("$T_{", i, "}$"))
    ) + ylab("") +
    theme_minimal() +
    theme(plot.title = element_text(size = 10, face = "bold"))
}
plots[[1]] <- plots[[1]] + ylab("Density")
plots[[3]] <- plots[[3]] + ylab("Density")
plots[[5]] <- plots[[5]] + ylab("Density")
plots[[7]] <- plots[[7]] + ylab("Density")
plots[[9]] <- plots[[9]] + ylab("Density")

panel <- ggarrange(plotlist = plots, ncol = 2, nrow = ceiling(N / 2))
print(panel)
ggsave("mcmc_calibrations.pdf",panel,device="pdf")
