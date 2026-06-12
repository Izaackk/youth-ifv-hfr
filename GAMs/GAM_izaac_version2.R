library(R.matlab)
library(mgcv)
library(gratia)
library(visreg)
library(ggplot2)
library(paletteer)
library(scales)
library(dplyr)
library(pals)
library(tidyr)

rm(list = ls())

sa_mat <- ".../HFR_output/MatchRate1/SAaxis/SA_HFR87.mat"
input_mat <- ".../HFR_output/MatchRate1/Myelin/ROI/data_bin20_step10.mat"
outdir <- ".../HFR_output/MatchRate1/Myelin/ROI/GAM/k3"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

SA <- readMat(sa_mat)
schaefer400_SAaxis <- as.numeric(SA$SA.HFR)

data_all <- readMat(input_mat)
Data <- list()
Data$age.bins <- as.numeric(data_all$age.bins)
Data$dice.bins <- data_all$H
Data$sex.de.bins <- as.numeric(data_all$sex.de.bins)
Data$mFD.bins <- as.numeric(data_all$mFD.bins)

Age <- Data$age.bins
Sex <- Data$sex.de.bins
mFD <- Data$mFD.bins

n_region <- 87
p_value <- numeric(n_region)
F_value <- numeric(n_region)
R_sq <- numeric(n_region)
delta_R_sq <- numeric(n_region)
anova_pvalue <- numeric(n_region)
sig <- numeric(n_region)
lm_c <- numeric(n_region)
all_data <- data.frame()
all_derv_data <- data.frame()

Age_points <- data.frame(Age = seq(8, 22, by = 1))
Age_points$Sex <- mean(Data$sex.de.bins, na.rm = TRUE)
Age_points$mFD <- mean(Data$mFD.bins, na.rm = TRUE)

for (i in seq_len(n_region)) {
  cat("Processing ROI: ", i, "/", n_region, "\n", sep = "")
  y <- Data$dice.bins[, i]

  mod_gam <- gam(y ~ s(Age, bs = "cs", k = 3) + Sex + mFD,
                 data = Data, method = "REML", na.action = "na.omit")
  gam_summ <- summary(mod_gam)

  p_value[i] <- gam_summ$s.table[1, "p-value"]
  F_value[i] <- gam_summ$s.table[1, "F"]
  R_sq[i] <- gam_summ$r.sq

  gam_null <- gam(y ~ Sex + mFD, data = Data, method = "REML", na.action = "na.omit")
  r2_null <- summary(gam_null)$r.sq
  delta_temp <- gam_summ$r.sq - r2_null

  derv <- derivatives(mod_gam, term = "s(Age)", type = "central", data = Age_points)
  derv$region <- i
  derv$SA.axis <- schaefer400_SAaxis[i]
  all_derv_data <- rbind(all_derv_data, derv)

  mean_deriv <- mean(derv$.derivative, na.rm = TRUE)
  delta_R_sq[i] <- ifelse(mean_deriv < 0, -delta_temp, delta_temp)

  anova_p <- anova.gam(gam_null, mod_gam, test = "Chisq")$`Pr(>Chi)`[2]
  anova_pvalue[i] <- anova_p
  sig[i] <- ifelse(anova_p < 0.05, 1, 0)

  Gam_Z_Age <- qnorm(anova_p / 2, lower.tail = FALSE)
  lm_mod <- lm(y ~ Age + Sex + mFD, data = Data)
  lm_t <- summary(lm_mod)$coefficients["Age", "t value"]
  lm_c[i] <- summary(lm_mod)$coefficients["Age", "Estimate"]
  if (!is.na(lm_t) && lm_t < 0) {
    Gam_Z_Age <- -Gam_Z_Age
  }

  v <- visreg(mod_gam, "Age", plot = FALSE)
  current_data <- v$fit
  current_data$SA.axis <- schaefer400_SAaxis[i]
  current_data$region <- i
  current_data$anova_p <- anova_p
  current_data$p_value <- p_value[i]
  current_data$sig <- sig[i]
  current_data$Gam_Z_Age <- Gam_Z_Age
  current_data$delta_R_sq <- delta_R_sq[i]
  all_data <- rbind(all_data, current_data)
}

GAM_Result <- data.frame(
  p_value,
  F_value,
  R_sq,
  delta_R_sq,
  anova_pvalue,
  sig,
  age_lm_c = lm_c
)
colnames(GAM_Result) <- c(
  "p_value",
  "F_value",
  "R_sq",
  "delta_R_sq",
  "anova_pvalue",
  "anova_p0.05",
  "age_lm_c"
)

write.csv(GAM_Result, file.path(outdir, "GAM_Result.csv"), row.names = FALSE, quote = FALSE)
cat("Saved GAM results: ", file.path(outdir, "GAM_Result.csv"), "\n", sep = "")

if (".value" %in% colnames(all_derv_data)) {
  all_derv_data <- all_derv_data %>% rename(Age = .value)
} else if ("data" %in% colnames(all_derv_data)) {
  all_derv_data <- all_derv_data %>% rename(Age = data)
} else if ("cond" %in% colnames(all_derv_data)) {
  all_derv_data <- all_derv_data %>% rename(Age = cond)
}

write.csv(all_derv_data, file.path(outdir, "GAM_Derivatives_Age.csv"),
          row.names = FALSE, quote = FALSE)

all_derv_data_slim <- all_derv_data %>%
  dplyr::select(Age, region, .derivative) %>%
  filter(!is.na(.derivative))

deriv_wide <- tidyr::pivot_wider(
  all_derv_data_slim,
  names_from = region,
  values_from = .derivative
)
deriv_wide <- deriv_wide[order(deriv_wide$Age), ]

write.csv(deriv_wide, file.path(outdir, "GAM_Derivatives_AgeMatrix.csv"),
          row.names = FALSE, quote = FALSE)
cat("Saved derivative matrix: ", file.path(outdir, "GAM_Derivatives_AgeMatrix.csv"), "\n", sep = "")

all_data <- all_data[all_data$SA.axis != 0, ]
colormap_limit <- range(schaefer400_SAaxis, na.rm = TRUE)

All_net_figs <- ggplot(all_data, aes(x = Age, y = visregFit, group = SA.axis)) +
  geom_line(linewidth = 1, alpha = 0.7, aes(color = as.numeric(SA.axis))) +
  paletteer::scale_color_paletteer_c(
    "grDevices::PuOr",
    direction = 1,
    limits = colormap_limit,
    oob = squish
  ) +
  ylim(0, 80) +
  xlim(7, 22) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    panel.background = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 0.1, r = 0, b = 0, l = 1, unit = "cm")
  )

print(All_net_figs)
ggsave(file.path(outdir, "GAM_age_develop_color_SA.pdf"), plot = All_net_figs)
ggsave(file.path(outdir, "GAM_age_develop_color_SA.tiff"), plot = All_net_figs)

derv_plot <- ggplot(all_derv_data, aes(x = Age, y = .derivative, group = region, color = SA.axis)) +
  geom_line(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = 2) +
  paletteer::scale_color_paletteer_c("grDevices::PuOr", limits = colormap_limit, oob = squish) +
  theme_classic() +
  labs(x = "Age", y = "d(Measure)/dAge", title = "GAM Estimated Derivative by Age")

print(derv_plot)
ggsave(file.path(outdir, "GAM_age_derivative_SA.pdf"), plot = derv_plot, width = 7, height = 5)
ggsave(file.path(outdir, "GAM_age_derivative_SA.tiff"), plot = derv_plot, width = 7, height = 5)
