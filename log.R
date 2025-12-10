# Load required packages
library(dplyr)
library(signal)
library(tibble)
library(ggplot2)
library(tidyr)

# 1. Define HFR Function
compute_HFR <- function(rmsf_vector, cutoff_ratio) {
  n <- length(rmsf_vector)
  if(n < 2) return(NA)
  fft_result <- fft(rmsf_vector)
  magnitude <- Mod(fft_result)[1:floor(n/2)]
  if(length(magnitude) == 0) return(NA)
  freqs <- seq(0, 1, length.out = length(magnitude))
  cutoff_index <- floor(length(freqs) * cutoff_ratio)
  if (cutoff_index >= length(magnitude)) cutoff_index <- length(magnitude) - 1
  high_freq_sum <- sum(magnitude[(cutoff_index + 1):length(magnitude)], na.rm = TRUE)
  total_sum <- sum(magnitude, na.rm = TRUE)
  if (total_sum == 0) return(0)
  return(high_freq_sum / total_sum)
}

# 2. Setup
cutoff_ratios <- seq(0, 0.99, by = 0.01) # 100 Cutoffs
folder_path <- "/Users/alperkaragol/Desktop/RMSF" # Your path
file_list <- list.files(path = folder_path, pattern = "\\.tsv$", full.names = TRUE)
results <- tibble(file = character(), HFR = numeric(), ProteinLength = numeric(), Cutoff = numeric())

# 3. Process Files
print("Processing data...")
for (file in file_list) {
  df <- tryCatch(read.delim(file, header = TRUE), error=function(e) NULL)
  if (!is.null(df)) {
    rmsf_cols <- df[, grepl("RMSF", names(df)), drop = FALSE]
    if(ncol(rmsf_cols) > 0) {
      rmsf_vector <- apply(rmsf_cols, 1, median)
      p_len <- length(rmsf_vector)
      for (ratio in cutoff_ratios) {
        hfr_value <- compute_HFR(rmsf_vector, cutoff_ratio = ratio)
        results <- results %>% add_row(file = basename(file), HFR = hfr_value, ProteinLength = p_len, Cutoff = ratio)
      }
    }
  }
}

# 4. ADVANCED TRANSFORMATION: Log-Log Scaling
# We use dplyr::filter to avoid conflicts with the signal package
results_log <- results %>%
  dplyr::filter(HFR > 0) %>%
  mutate(
    Log_Length = log10(ProteinLength),
    Log_HFR = log10(HFR)
  )

# 5. Calculate Correlations on TRANSFORMED Data
cor_data <- results_log %>%
  group_by(Cutoff) %>%
  summarize(
    # Compute correlation on Log-Log data
    cor_val = cor(Log_Length, Log_HFR, use = "complete.obs"),
    p_val = cor.test(Log_Length, Log_HFR)$p.value,
    .groups = 'drop'
  ) %>%
  mutate(
    # Label now reflects the improved correlation
    label_text = paste0("R=", round(cor_val, 2), "\np=", format.pval(p_val, digits = 1, eps = 0.001))
  )

# 6. Plotting the Power Law (Log-Log) Relationship
p <- ggplot(results_log, aes(x = Log_Length, y = Log_HFR)) +
  geom_point(alpha = 0.3, color = "#006400", size = 0.5) + # Dark Green for transformed
  geom_smooth(method = "lm", color = "magenta", se = FALSE, size = 0.5) + # Linear fit on Log data
  
  # 10x10 Grid
  facet_wrap(~Cutoff, scales = "free_y", ncol = 10, labeller = label_both) +
  
  # Add Stats
  geom_text(data = cor_data, aes(label = label_text, x = Inf, y = Inf), 
            hjust = 1.1, vjust = 1.2, size = 1.8, fontface = "bold", inherit.aes = FALSE) +
  
  theme_bw() +
  labs(title = "Power Law Scaling: Log(Length) vs Log(HFR) across 100 Cutoffs",
       subtitle = "Data transformed using log10 to linearize the decay relationship",
       x = "Log10(Protein Length)", 
       y = "Log10(HFR)") +
  
  theme(
    strip.text = element_text(size = 5),
    axis.text = element_text(size = 4),
    axis.title = element_text(size = 9, face="bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p)

# Save high-res version
ggsave("~/Desktop/HFR_LogLog_Analysis.pdf", p, width = 20, height = 20, limitsize = FALSE)
print("Saved transformed analysis to Desktop/HFR_LogLog_Analysis.pdf")
