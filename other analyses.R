# ============================================================
#  EXTENDED HFR ANALYSIS — Novel Mathematical & Topological
#  Tests on Protein RMSF Signals
# ============================================================
# Original: HFR (FFT High-Frequency Ratio) vs Protein Length
# Added:
#  A. Spectral Entropy            (information theory)
#  B. Higuchi Fractal Dimension   (fractal geometry)
#  C. Hurst Exponent              (long-range memory / self-similarity)
#  D. Sample Entropy              (non-linear dynamics / regularity)
#  E. Wavelet Energy Ratio        (multi-resolution HFR analogue)
#  F. Persistent Homology (TDA)   (topological data analysis)
#     - H0 total persistence (connected component lifetimes)
#     - H1 max persistence   (loop/cycle lifetime)
#  G. Spectral Edge Frequency     (novel frequency-mass centre)
#  H. Multi-metric regression     (all predictors → length)
# ============================================================

library(dplyr)
library(signal)       # FFT
library(tibble)
library(ggplot2)
library(tidyr)
library(patchwork)    # plot layout
library(mgcv)         # GAM
library(scales)       # colour scales

# ── Optional: TDA (install if absent) ─────────────────────
if (!requireNamespace("TDA", quietly = TRUE)) {
  install.packages("TDA", repos = "https://cloud.r-project.org")
}
library(TDA)

# ── Optional: wavelets ────────────────────────────────────
if (!requireNamespace("wavelets", quietly = TRUE)) {
  install.packages("wavelets", repos = "https://cloud.r-project.org")
}
library(wavelets)


# ============================================================
#  SECTION 0 — Existing HFR function (unchanged)
# ============================================================
compute_HFR <- function(rmsf_vector, cutoff_ratio) {
  n <- length(rmsf_vector)
  if (n < 2) return(NA)
  fft_result  <- fft(rmsf_vector)
  magnitude   <- Mod(fft_result)[1:floor(n / 2)]
  if (length(magnitude) == 0) return(NA)
  freqs        <- seq(0, 1, length.out = length(magnitude))
  cutoff_index <- floor(length(freqs) * cutoff_ratio)
  if (cutoff_index >= length(magnitude)) cutoff_index <- length(magnitude) - 1
  high_freq_sum <- sum(magnitude[(cutoff_index + 1):length(magnitude)], na.rm = TRUE)
  total_sum     <- sum(magnitude, na.rm = TRUE)
  if (total_sum == 0) return(0)
  high_freq_sum / total_sum
}


# ============================================================
#  SECTION 1 — Novel Metric Functions
# ============================================================

# ── A. SPECTRAL ENTROPY ───────────────────────────────────
# Measures how "spread" or "peaked" the power spectrum is.
# A flat spectrum → high entropy (many frequencies contribute equally).
# A peaked spectrum → low entropy (dominated by a few frequencies).
# Novel link to biology: highly disordered proteins are expected to
# show flatter spectra (higher spectral entropy) than compact globular
# proteins, independently of simple HFR.
compute_spectral_entropy <- function(rmsf_vector) {
  n <- length(rmsf_vector)
  if (n < 4) return(NA)
  fft_result <- fft(rmsf_vector)
  power      <- (Mod(fft_result)[1:floor(n / 2)])^2
  power      <- power / sum(power)          # normalise to probability
  power      <- power[power > 0]
  -sum(power * log2(power))                 # Shannon entropy in bits
}


# ── B. HIGUCHI FRACTAL DIMENSION ─────────────────────────
# Measures the self-similarity / complexity of the RMSF signal.
# D ≈ 1 → smooth signal; D ≈ 2 → maximally rough (space-filling).
# Hypothesis: longer proteins may show higher-dimensional (more
# complex) flexibility patterns due to multi-domain dynamics.
compute_higuchi_fd <- function(x, kmax = 6) {
  n <- length(x)
  if (n < 2 * kmax) return(NA)
  Lk <- numeric(kmax)
  for (k in 1:kmax) {
    Lmk <- numeric(k)
    for (m in 1:k) {
      idx   <- seq(m, n, by = k)
      if (length(idx) < 2) { Lmk[m] <- NA; next }
      Nm    <- floor((n - m) / k)
      Lmk[m] <- (sum(abs(diff(x[idx]))) * (n - 1)) / (k * Nm * k)
    }
    Lk[k] <- mean(Lmk, na.rm = TRUE)
  }
  valid <- which(Lk > 0)
  if (length(valid) < 2) return(NA)
  fit <- lm(log(Lk[valid]) ~ log(1 / valid))
  coef(fit)[2]   # slope = fractal dimension
}


# ── C. HURST EXPONENT (R/S method) ──────────────────────
# H > 0.5 → persistent (trending) fluctuations
# H = 0.5 → random walk (Brownian motion)
# H < 0.5 → anti-persistent (mean-reverting) fluctuations
# Biologically, H captures whether local flexibility is correlated
# across the sequence — a property invisible to FFT-based HFR.
compute_hurst <- function(x) {
  n <- length(x)
  if (n < 20) return(NA)
  sizes <- unique(floor(exp(seq(log(10), log(floor(n / 2)), length.out = 8))))
  sizes <- sizes[sizes >= 8]
  if (length(sizes) < 3) return(NA)
  rs_vals <- sapply(sizes, function(s) {
    starts <- seq(1, n - s + 1, by = floor(s / 2))
    rs_sub <- sapply(starts, function(i) {
      seg   <- x[i:(i + s - 1)]
      seg_d <- seg - cumsum(rep(mean(seg), s)) / seq_along(seg)
      # deviation from mean
      seg_d <- cumsum(seg - mean(seg))
      R     <- max(seg_d) - min(seg_d)
      S     <- sd(seg)
      if (S == 0) return(NA)
      R / S
    })
    mean(rs_sub, na.rm = TRUE)
  })
  valid <- which(!is.na(rs_vals) & rs_vals > 0)
  if (length(valid) < 2) return(NA)
  fit <- lm(log(rs_vals[valid]) ~ log(sizes[valid]))
  coef(fit)[2]
}


# ── D. SAMPLE ENTROPY (SampEn) ───────────────────────────
# Measures the unpredictability of the RMSF time series.
# Low SampEn → regular, predictable fluctuations (rigid protein).
# High SampEn → irregular, complex fluctuations (flexible protein).
# Unlike Higuchi FD, SampEn is template-match based (nonlinear).
compute_sample_entropy <- function(x, m = 2, r_frac = 0.2) {
  n   <- length(x)
  if (n < (m + 2) * 2) return(NA)
  r   <- r_frac * sd(x, na.rm = TRUE)
  if (is.na(r) || r == 0) return(NA)
  
  count_templates <- function(m_len) {
    count <- 0
    for (i in 1:(n - m_len)) {
      for (j in (i + 1):(n - m_len)) {
        if (max(abs(x[i:(i + m_len - 1)] - x[j:(j + m_len - 1)])) < r)
          count <- count + 1
      }
    }
    count
  }
  
  # For speed, subsample if very long
  if (n > 300) {
    set.seed(42)
    x <- x[sort(sample(n, 300))]
    n <- 300
  }
  
  A <- count_templates(m + 1)
  B <- count_templates(m)
  if (B == 0) return(NA)
  -log(A / B)
}


# ── E. WAVELET ENERGY RATIO ──────────────────────────────
# Analogous to HFR but using the Discrete Wavelet Transform (DWT).
# Each decomposition level corresponds to a different frequency band.
# The ratio of detail (high-frequency) to approximation (low-freq)
# energy provides a multi-resolution view beyond FFT.
compute_wavelet_energy_ratio <- function(x, wf = "la8", n_levels = 4) {
  n <- length(x)
  if (n < 2^n_levels * 2) return(NA)
  # Pad to power of 2
  n_pad <- 2^ceiling(log2(n))
  x_pad <- c(x, rep(0, n_pad - n))
  tryCatch({
    wt      <- wavelets::dwt(x_pad, filter = wf, n.levels = n_levels, boundary = "periodic")
    detail_energy <- sum(sapply(wt@W, function(d) sum(d^2)))
    approx_energy <- sum(wt@V[[n_levels]]^2)
    total         <- detail_energy + approx_energy
    if (total == 0) return(NA)
    detail_energy / total
  }, error = function(e) NA)
}


# ── F. TOPOLOGICAL DATA ANALYSIS (Persistent Homology) ───
# The RMSF signal is embedded in 2D delay coordinates [x(t), x(t+τ)],
# forming a point cloud in ℝ². We compute the Vietoris-Rips filtration
# and extract:
#   • H0 total persistence: total "lifespan" of connected components
#     (captures how many distinct flexibility clusters exist)
#   • H1 max persistence:   maximum loop lifetime
#     (captures cyclical / periodic structure in the RMSF pattern)
# This is a direct application of topological data analysis (TDA)
# to protein flexibility — a genuinely novel approach.
compute_tda_features <- function(x, tau = 3, max_pts = 200) {
  n <- length(x)
  if (n < tau + 5) return(list(H0_total = NA, H1_max = NA, H0_n = NA))
  
  # Delay embedding
  idx1 <- 1:(n - tau)
  idx2 <- (1 + tau):n
  pts  <- cbind(x[idx1], x[idx2])
  
  # Subsample for speed
  if (nrow(pts) > max_pts) {
    set.seed(42)
    pts <- pts[sample(nrow(pts), max_pts), ]
  }
  
  tryCatch({
    diag <- TDA::ripsDiag(pts, maxdimension = 1, maxscale = diff(range(pts)) * 1.5,
                          library = "GUDHI", printProgress = FALSE)$diagram
    
    # H0 features
    h0    <- diag[diag[, "dimension"] == 0, , drop = FALSE]
    # Remove the one class that lives forever (born=0, death=Inf)
    h0    <- h0[is.finite(h0[, "Death"]), , drop = FALSE]
    H0_total <- if (nrow(h0) > 0) sum(h0[, "Death"] - h0[, "Birth"]) else 0
    H0_n     <- nrow(h0)   # number of transient components (flexibility clusters)
    
    # H1 features
    h1    <- diag[diag[, "dimension"] == 1, , drop = FALSE]
    H1_max <- if (nrow(h1) > 0) max(h1[, "Death"] - h1[, "Birth"]) else 0
    
    list(H0_total = H0_total, H1_max = H1_max, H0_n = H0_n)
  }, error = function(e) list(H0_total = NA, H1_max = NA, H0_n = NA))
}


# ── G. SPECTRAL EDGE FREQUENCY ───────────────────────────
# The frequency below which X% of the total spectral power lies.
# Borrowed from EEG analysis — applied here to RMSF for the first time.
# A protein with lower spectral edge → most power at low frequencies
# (slow, large-scale correlated motions).
compute_spectral_edge <- function(x, edge_frac = 0.95) {
  n   <- length(x)
  if (n < 4) return(NA)
  fft_result <- fft(x)
  power      <- (Mod(fft_result)[1:floor(n / 2)])^2
  total      <- sum(power)
  if (total == 0) return(NA)
  cumsum_power <- cumsum(power) / total
  idx <- which(cumsum_power >= edge_frac)[1]
  if (is.na(idx)) return(1)
  idx / length(power)   # normalised frequency (0–1)
}


# ============================================================
#  SECTION 2 — File Processing
# ============================================================
folder_path  <- "/Users/alperkaragol/Desktop/RMSF"   # ← adjust path
file_list    <- list.files(path = folder_path, pattern = "\\.tsv$", full.names = TRUE)
cutoff_ratios <- seq(0, 0.99, by = 0.01)

# Storage for HFR sweep (same as original)
results_hfr <- tibble(file = character(), HFR = numeric(),
                      ProteinLength = numeric(), Cutoff = numeric())

# Storage for single-value metrics per protein
results_novel <- tibble(
  file            = character(),
  ProteinLength   = numeric(),
  SpectralEntropy = numeric(),
  HiguchiFD       = numeric(),
  HurstExp        = numeric(),
  SampleEnt       = numeric(),
  WaveletER       = numeric(),
  SpectralEdge95  = numeric(),
  TDA_H0_total    = numeric(),
  TDA_H1_max      = numeric(),
  TDA_H0_n        = numeric()
)

cat(sprintf("Processing %d files...\n", length(file_list)))

for (file in file_list) {
  df <- tryCatch(read.delim(file, header = TRUE), error = function(e) NULL)
  if (is.null(df)) next
  
  rmsf_cols <- df[, grepl("RMSF", names(df)), drop = FALSE]
  if (ncol(rmsf_cols) == 0) next
  
  df$RMSF_median <- apply(rmsf_cols, 1, median)
  x   <- df$RMSF_median
  n   <- length(x)
  fn  <- basename(file)
  
  # ── HFR across 100 cutoffs (original) ──
  for (ratio in cutoff_ratios) {
    results_hfr <- results_hfr %>%
      add_row(file = fn, HFR = compute_HFR(x, ratio),
              ProteinLength = n, Cutoff = ratio)
  }
  
  # ── Novel single metrics ──
  cat(sprintf("  → Novel metrics: %s\n", fn))
  tda  <- compute_tda_features(x)
  results_novel <- results_novel %>%
    add_row(
      file            = fn,
      ProteinLength   = n,
      SpectralEntropy = compute_spectral_entropy(x),
      HiguchiFD       = compute_higuchi_fd(x),
      HurstExp        = compute_hurst(x),
      SampleEnt       = compute_sample_entropy(x),
      WaveletER       = compute_wavelet_energy_ratio(x),
      SpectralEdge95  = compute_spectral_edge(x, 0.95),
      TDA_H0_total    = tda$H0_total,
      TDA_H1_max      = tda$H1_max,
      TDA_H0_n        = tda$H0_n
    )
}


# ============================================================
#  SECTION 3 — Correlation Analysis of Novel Metrics
# ============================================================
novel_metrics <- c("SpectralEntropy","HiguchiFD","HurstExp","SampleEnt",
                   "WaveletER","SpectralEdge95","TDA_H0_total","TDA_H1_max","TDA_H0_n")

metric_labels <- c(
  SpectralEntropy = "Spectral Entropy\n(bits)",
  HiguchiFD       = "Higuchi Fractal\nDimension",
  HurstExp        = "Hurst Exponent\n(R/S method)",
  SampleEnt       = "Sample Entropy\n(m=2)",
  WaveletER       = "Wavelet Energy\nRatio (DWT)",
  SpectralEdge95  = "Spectral Edge\nFrequency (95%)",
  TDA_H0_total    = "TDA H₀\nTotal Persistence",
  TDA_H1_max      = "TDA H₁\nMax Persistence",
  TDA_H0_n        = "TDA H₀\nCluster Count"
)

# Pearson + Spearman + Kendall for each metric
cor_table <- purrr::map_dfr(novel_metrics, function(m) {
  sub <- results_novel %>% select(ProteinLength, all_of(m)) %>% drop_na()
  if (nrow(sub) < 5) return(NULL)
  y <- sub[[m]]
  x <- sub$ProteinLength
  p <- cor.test(x, y, method = "pearson")
  s <- cor.test(x, y, method = "spearman")
  k <- cor.test(x, y, method = "kendall")
  tibble(
    Metric   = m,
    Label    = metric_labels[m],
    Pearson_R  = round(p$estimate, 3),
    Pearson_p  = round(p$p.value, 4),
    Spearman_r = round(s$estimate, 3),
    Spearman_p = round(s$p.value, 4),
    Kendall_tau= round(k$estimate, 3),
    Kendall_p  = round(k$p.value, 4)
  )
})

cat("\n===== NOVEL METRIC CORRELATIONS WITH PROTEIN LENGTH =====\n")
print(cor_table, n = Inf)
write.csv(cor_table, "~/Desktop/novel_metric_correlations.csv", row.names = FALSE)


# ============================================================
#  SECTION 4 — GAM Non-Linear Regressions
# ============================================================
# Generalised Additive Models allow detection of non-linear
# relationships that simple Pearson R would miss.
gam_results <- purrr::map(novel_metrics, function(m) {
  sub <- results_novel %>% select(ProteinLength, all_of(m)) %>% drop_na()
  if (nrow(sub) < 10) return(NULL)
  fml <- as.formula(paste(m, "~ s(ProteinLength, k=5)"))
  tryCatch(mgcv::gam(fml, data = sub, method = "REML"), error = function(e) NULL)
})
names(gam_results) <- novel_metrics

cat("\n===== GAM SUMMARIES (Non-linear effects) =====\n")
for (m in novel_metrics) {
  if (!is.null(gam_results[[m]])) {
    cat("\n---", m, "---\n")
    print(summary(gam_results[[m]])$s.table)
  }
}


# ============================================================
#  SECTION 5 — Multiple Regression: Do novel metrics predict
#              protein length better together than HFR alone?
# ============================================================
# Use the HFR at the two most significant cutoffs (0.0 and 0.1)
hfr_wide <- results_hfr %>%
  filter(Cutoff %in% c(0.00, 0.10)) %>%
  mutate(varname = paste0("HFR_cut", Cutoff)) %>%
  select(file, varname, HFR) %>%
  pivot_wider(names_from = varname, values_from = HFR)

combined <- results_novel %>%
  left_join(hfr_wide, by = "file") %>%
  drop_na()

if (nrow(combined) >= 10) {
  mod_hfr_only  <- lm(ProteinLength ~ HFR_cut0 + HFR_cut0.1, data = combined)
  mod_full      <- lm(ProteinLength ~ HFR_cut0 + HFR_cut0.1 +
                        SpectralEntropy + HiguchiFD + HurstExp +
                        SampleEnt + WaveletER + SpectralEdge95 +
                        TDA_H0_total + TDA_H1_max + TDA_H0_n,
                      data = combined)
  
  cat("\n===== REGRESSION: HFR only =====\n")
  print(summary(mod_hfr_only))
  cat("\n===== REGRESSION: Full (HFR + Novel) =====\n")
  print(summary(mod_full))
  cat("\nModel comparison (ANOVA):\n")
  print(anova(mod_hfr_only, mod_full))
}


# ============================================================
#  SECTION 6 — Topological Phase Plot
#              H0 total persistence vs H1 max persistence,
#              coloured by protein length — a topological map
#              of the RMSF signal landscape
# ============================================================
p_topo <- ggplot(results_novel %>% drop_na(TDA_H0_total, TDA_H1_max),
                 aes(x = TDA_H0_total, y = TDA_H1_max, colour = ProteinLength)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_colour_viridis_c(name = "Protein\nLength", option = "plasma") +
  geom_smooth(method = "lm", colour = "black", se = TRUE, linewidth = 0.7) +
  labs(
    title    = "Topological Phase Portrait of RMSF Signals",
    subtitle = "Vietoris-Rips persistent homology on delay-embedded RMSF\nH₀ = connected component persistence | H₁ = loop persistence",
    x        = "H₀ Total Persistence (flexibility cluster lifetimes)",
    y        = "H₁ Max Persistence (cyclic structure strength)"
  ) +
  theme_bw(base_size = 13) +
  theme(plot.subtitle = element_text(size = 9, colour = "grey40"))


# ============================================================
#  SECTION 7 — Composite 9-Panel Figure: Novel Metrics vs Length
# ============================================================
make_panel <- function(metric, label) {
  sub <- results_novel %>% select(ProteinLength, y = all_of(metric)) %>% drop_na()
  cr  <- cor.test(sub$ProteinLength, sub$y, method = "pearson")
  ann <- sprintf("R=%.2f\np=%s", cr$estimate, format.pval(cr$p.value, digits=1, eps=0.001))
  
  ggplot(sub, aes(x = ProteinLength, y = y)) +
    geom_point(alpha = 0.4, colour = "steelblue4", size = 1.8) +
    geom_smooth(method = "gam", formula = y ~ s(x, k=5),
                colour = "firebrick", fill = "firebrick", alpha = 0.15, linewidth = 0.8) +
    annotate("text", x = Inf, y = Inf, label = ann,
             hjust = 1.1, vjust = 1.3, size = 3.5, fontface = "bold") +
    labs(title = label, x = "Protein Length", y = NULL) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(size = 10, face = "bold"))
}

panels <- mapply(make_panel, novel_metrics, metric_labels, SIMPLIFY = FALSE)
fig_novel <- patchwork::wrap_plots(panels, ncol = 3) +
  patchwork::plot_annotation(
    title    = "Novel Signal & Topological Metrics vs Protein Length",
    subtitle = "GAM smooths shown (red). Pearson R annotated per panel.",
    theme    = theme(plot.title    = element_text(size = 15, face = "bold"),
                     plot.subtitle = element_text(size = 10, colour = "grey40"))
  )


# ============================================================
#  SECTION 8 — HFR Signature Heatmap: R-value across 100 cutoffs
#              (enhanced version of the original PDF figure)
# ============================================================
cor_hfr <- results_hfr %>%
  group_by(Cutoff) %>%
  summarize(
    R_pearson  = cor(ProteinLength, HFR, use = "complete.obs"),
    p_val      = cor.test(ProteinLength, HFR)$p.value,
    R_spearman = cor(ProteinLength, HFR, method = "spearman", use = "complete.obs"),
    .groups = "drop"
  ) %>%
  mutate(
    sig        = case_when(p_val < 0.001 ~ "***",
                           p_val < 0.01  ~ "**",
                           p_val < 0.05  ~ "*",
                           TRUE          ~ ""),
    row_band   = factor(floor(Cutoff * 10) / 10),
    col_pos    = round((Cutoff %% 0.1) * 100) / 10
  )

p_heatmap <- ggplot(cor_hfr, aes(x = Cutoff, y = R_pearson)) +
  geom_col(aes(fill = R_pearson), width = 0.009) +
  geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed") +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 0, name = "Pearson R") +
  geom_text(aes(label = sig, y = R_pearson - 0.005),
            size = 2.5, vjust = 1, colour = "black") +
  labs(title    = "Pearson R (HFR vs Protein Length) across 100 Cutoffs",
       subtitle = "* p<0.05  ** p<0.01  *** p<0.001",
       x = "FFT Cutoff Ratio", y = "Pearson R") +
  theme_bw(base_size = 12)

p_spearman <- ggplot(cor_hfr, aes(x = Cutoff, y = R_spearman)) +
  geom_col(aes(fill = R_spearman), width = 0.009) +
  geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed") +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 0, name = "Spearman ρ") +
  labs(title    = "Spearman ρ (HFR vs Protein Length) across 100 Cutoffs",
       subtitle = "Non-parametric rank correlation — robust to outliers",
       x = "FFT Cutoff Ratio", y = "Spearman ρ") +
  theme_bw(base_size = 12)


# ============================================================
#  SECTION 9 — Mathematical Novelty: Scale-Free Behaviour Test
#  Tests whether HFR ~ Length^(-α) (power-law decay),
#  which would imply scale-free dynamics of protein flexibility.
#  We fit a log-log linear model at each cutoff and extract α.
# ============================================================
powerlaw_fits <- results_hfr %>%
  filter(HFR > 0, ProteinLength > 0) %>%
  group_by(Cutoff) %>%
  summarize(
    alpha  = tryCatch({
      fit <- lm(log(HFR) ~ log(ProteinLength))
      coef(fit)[2]
    }, error = function(e) NA),
    R2     = tryCatch({
      fit <- lm(log(HFR) ~ log(ProteinLength))
      summary(fit)$r.squared
    }, error = function(e) NA),
    .groups = "drop"
  )

p_powerlaw <- ggplot(powerlaw_fits, aes(x = Cutoff, y = alpha)) +
  geom_line(colour = "darkviolet", linewidth = 1) +
  geom_ribbon(aes(ymin = alpha - 0.05, ymax = alpha + 0.05),
              fill = "orchid", alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  labs(
    title    = "Power-Law Exponent α: HFR ~ Length^(−α)",
    subtitle = "If α ≠ 0 and R² is high → scale-free flexibility dynamics",
    x        = "FFT Cutoff Ratio",
    y        = expression(paste("Exponent ", alpha))
  ) +
  theme_bw(base_size = 12)

p_powerlaw_r2 <- ggplot(powerlaw_fits, aes(x = Cutoff, y = R2)) +
  geom_line(colour = "darkorange", linewidth = 1) +
  geom_area(fill = "orange", alpha = 0.25) +
  labs(title    = "Power-Law Fit Quality (R²) across Cutoffs",
       subtitle = "High R² → log-log linearity supports scale-free hypothesis",
       x = "FFT Cutoff Ratio", y = "R² of log-log fit") +
  ylim(0, 1) +
  theme_bw(base_size = 12)


# ============================================================
#  SECTION 10 — Save All Outputs
# ============================================================
ggsave("~/Desktop/HFR_novel_9panel.pdf",   fig_novel,    width = 18, height = 18, limitsize = FALSE)
ggsave("~/Desktop/HFR_topo_phase.pdf",     p_topo,       width = 9,  height = 7)
ggsave("~/Desktop/HFR_heatmap_R.pdf",
       p_heatmap / p_spearman,             width = 14, height = 10)
ggsave("~/Desktop/HFR_powerlaw.pdf",
       p_powerlaw / p_powerlaw_r2,         width = 12, height = 10)

# Original 10x10 facet plot (unchanged)
cor_hfr_label <- cor_hfr %>%
  mutate(label_text = paste0("R=", round(R_pearson, 2), "\np=",
                             format.pval(p_val, digits = 1, eps = 0.001)))

p_orig <- ggplot(results_hfr, aes(x = ProteinLength, y = HFR)) +
  geom_point(alpha = 0.3, colour = "darkblue", size = 0.5) +
  geom_smooth(method = "lm", colour = "red", se = FALSE, linewidth = 0.5) +
  facet_wrap(~Cutoff, scales = "free_y", ncol = 10, labeller = label_both) +
  geom_text(data = cor_hfr_label,
            aes(label = label_text, x = Inf, y = Inf),
            hjust = 1.1, vjust = 1.2, size = 1.8, fontface = "bold", inherit.aes = FALSE) +
  theme_bw() +
  labs(title = "HFR vs Protein Length across 100 Cutoff Values (0.00 – 0.99)",
       x = "Length", y = "HFR") +
  theme(strip.text = element_text(size = 5), axis.text = element_text(size = 4),
        axis.title = element_text(size = 8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave("~/Desktop/HFR_100_Cutoffs.pdf", p_orig, width = 20, height = 20, limitsize = FALSE)

cat("\n✅  All outputs saved to ~/Desktop/\n")
cat("   • HFR_100_Cutoffs.pdf           — original 10×10 grid (unchanged)\n")
cat("   • HFR_novel_9panel.pdf          — 9 novel metrics vs protein length\n")
cat("   • HFR_topo_phase.pdf            — topological phase portrait (TDA)\n")
cat("   • HFR_heatmap_R.pdf             — Pearson & Spearman across 100 cutoffs\n")
cat("   • HFR_powerlaw.pdf              — scale-free exponent analysis\n")
cat("   • novel_metric_correlations.csv — full correlation table\n")