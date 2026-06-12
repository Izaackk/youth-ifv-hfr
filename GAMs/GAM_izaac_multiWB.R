library(R.matlab)
library(mgcv)
library(gratia)
library(visreg)
library(ggplot2)

rm(list = ls())

outdir <- ".../HFR_output/MatchRate1/multiWB"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

datasets <- list(
  list(
    name = "vSC_WB",
    input_mat = ".../HFR_output/MatchRate1/SC_DWI_new/WB/data_bin20_step10.mat",
    measure_var = "var.bins",
    sex_var = "sex.bins"
  ),
  list(
    name = "vMSN_WB",
    input_mat = ".../HFR_output/MatchRate1/MBC/diff_MBC_minmaxnorm_novetices/Inter/data_bin20_step10.mat",
    measure_var = "mvar.all",
    sex_var = "sex.bins"
  ),
  list(
    name = "vMD_WB",
    input_mat = ".../HFR_output/MatchRate1/Myelin/WB/data_bin20_step10.mat",
    measure_var = "H",
    sex_var = "sex.de.bins"
  )
)

GAM_results <- data.frame()
all_fit_data <- data.frame()

for (cfg in datasets) {
  cat("Processing ", cfg$name, "\n", sep = "")
  data_all <- readMat(cfg$input_mat)

  Data <- list()
  Data$Age <- as.numeric(data_all$age.bins)
  Data$y <- as.numeric(data_all[[cfg$measure_var]])
  Data$Sex <- as.numeric(data_all[[cfg$sex_var]])
  Data$mFD <- as.numeric(data_all$mFD.bins)

  mod_gam <- gam(y ~ s(Age, bs = "cs", k = 3) + Sex,
                 data = Data, method = "REML", na.action = "na.omit")
  gam_summ <- summary(mod_gam)

  p_value <- gam_summ$s.table[1, "p-value"]
  F_value <- gam_summ$s.table[1, "F"]
  R_sq <- gam_summ$r.sq

  gam_null <- gam(y ~ Sex, data = Data, method = "REML", na.action = "na.omit")
  r2_null <- summary(gam_null)$r.sq
  delta_temp <- gam_summ$r.sq - r2_null

  derv <- derivatives(mod_gam, term = "s(Age)")
  mean_deriv <- mean(derv$.derivative, na.rm = TRUE)
  delta_R_sq <- ifelse(mean_deriv < 0, -delta_temp, delta_temp)

  anova_pvalue <- anova.gam(gam_null, mod_gam, test = "Chisq")$`Pr(>Chi)`[2]
  sig <- ifelse(anova_pvalue < 0.05, 1, 0)

  lm_mod <- lm(y ~ Age + Sex + mFD, data = Data)
  lm_summ <- summary(lm_mod)
  age_lm_c <- lm_summ$coefficients["Age", "Estimate"]
  lm_t <- lm_summ$coefficients["Age", "t value"]

  Gam_Z_Age <- qnorm(anova_pvalue / 2, lower.tail = FALSE)
  if (!is.na(lm_t) && lm_t < 0) {
    Gam_Z_Age <- -Gam_Z_Age
  }

  GAM_results <- rbind(
    GAM_results,
    data.frame(
      measure = cfg$name,
      p_value = p_value,
      F_value = F_value,
      R_sq = R_sq,
      delta_R_sq = delta_R_sq,
      anova_pvalue = anova_pvalue,
      anova_p0.05 = sig,
      age_lm_c = age_lm_c
    )
  )

  v <- visreg(mod_gam, "Age", plot = FALSE)
  current_data <- v$fit
  current_data$measure <- cfg$name
  current_data$anova_p <- anova_pvalue
  current_data$p_value <- p_value
  current_data$sig <- sig
  current_data$Gam_Z_Age <- Gam_Z_Age
  current_data$delta_R_sq <- delta_R_sq
  all_fit_data <- rbind(all_fit_data, current_data)

  myplot <- visreg(
    mod_gam,
    "Age",
    gg = TRUE,
    type = "conditional",
    scale = "response",
    col.point = "black",
    overlay = TRUE,
    partial = TRUE,
    rug = FALSE,
    ylab = cfg$name
  )
  myplot <- myplot +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme_classic() +
    scale_x_continuous(
      breaks = c(6, 8, 10, 12, 14, 16, 18, 20, 22),
      expand = c(0, 0.45)
    )

  print(myplot)
  ggsave(filename = paste0(cfg$name, "_age.jpg"), plot = myplot,
         path = outdir, width = 60, height = 60, units = "mm", device = "jpeg")
  ggsave(filename = paste0(cfg$name, "_age.pdf"), plot = myplot,
         path = outdir, width = 60, height = 60, units = "mm", device = "pdf")
}

write.csv(GAM_results, file.path(outdir, "GAM_Result.csv"), row.names = FALSE, quote = FALSE)
write.csv(all_fit_data, file.path(outdir, "GAM_Fitted_AgeCurves.csv"), row.names = FALSE, quote = FALSE)

cat("Saved GAM results: ", file.path(outdir, "GAM_Result.csv"), "\n", sep = "")
