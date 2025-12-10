# Load required packages
library(dplyr)
library(signal)  # for FFT
library(tibble)
library(ggplot2)
library(tidyr)

# 1. Define the function to compute HFR
compute_HFR <- function(rmsf_vector, cutoff_ratio) {
  n <- length(rmsf_vector)
  if(n < 2) return(NA)
  
  fft_result <- fft(rmsf_vector)
  magnitude <- Mod(fft_result)[1:floor(n/2)]  # keep positive frequencies
  
  if(length(magnitude) == 0) return(NA)
  
  freqs <- seq(0, 1, length.out = length(magnitude))
  
  # Cutoff: retain high-frequency part above the given ratio
  cutoff_index <- floor(length(freqs) * cutoff_ratio)
  
  if (cutoff_index >= length(magnitude)) {
    cutoff_index <- length(magnitude) - 1
  }
  
  # Sum magnitudes
  high_freq_sum <- sum(magnitude[(cutoff_index + 1):length(magnitude)], na.rm = TRUE)
  total_sum <- sum(magnitude, na.rm = TRUE)
  
  if (total_sum == 0) return(0)
  
  HFR <- high_freq_sum / total_sum
  return(HFR)
}

# 2. Setup Parameters and Path
# Generate 100 cutoff ratios from 0.00 to 0.99
cutoff_ratios <- seq(0, 0.99, by = 0.01)

# Set folder path
folder_path <- "/Users/alperkaragol/Desktop/RMSF"  # Replace with your path
file_list <- list.files(path = folder_path, pattern = "\\.tsv$", full.names = TRUE)

# Initialize result dataframe
results <- tibble(file = character(), 
                  HFR = numeric(), 
                  ProteinLength = numeric(), 
                  Cutoff = numeric())

# 3. Process each file
print(paste("Processing", length(file_list), "files across", length(cutoff_ratios), "cutoffs..."))

for (file in file_list) {
  df <- tryCatch(read.delim(file, header = TRUE), error=function(e) NULL)
  
  if (!is.null(df)) {
    rmsf_cols <- df[, grepl("RMSF", names(df)), drop = FALSE]
    
    if(ncol(rmsf_cols) > 0) {
      df$RMSF_median <- apply(rmsf_cols, 1, median)
      rmsf_vector <- df$RMSF_median
      p_len <- length(rmsf_vector)
      
      # Loop through all 100 cutoffs
      for (ratio in cutoff_ratios) {
        hfr_value <- compute_HFR(rmsf_vector, cutoff_ratio = ratio)
        
        results <- results %>%
          add_row(file = basename(file), 
                  HFR = hfr_value, 
                  ProteinLength = p_len, 
                  Cutoff = ratio)
      }
    }
  }
}

# 4. Calculate Correlations
cor_data <- results %>%
  group_by(Cutoff) %>%
  summarize(
    cor_val = cor(ProteinLength, HFR, use = "complete.obs"),
    p_val = cor.test(ProteinLength, HFR)$p.value,
    .groups = 'drop'
  ) %>%
  mutate(
    # Create compact label: R=0.XX, p=0.XX
    label_text = paste0("R=", round(cor_val, 2), "\np=", format.pval(p_val, digits = 1, eps = 0.001))
  )

# 5. Plotting (10x10 Grid)
p <- ggplot(results, aes(x = ProteinLength, y = HFR)) +
  geom_point(alpha = 0.3, color = "darkblue", size = 0.5) + # Smaller points
  geom_smooth(method = "lm", color = "red", se = FALSE, size = 0.5) + # Thinner line
  
  # 10x10 Facet Grid
  facet_wrap(~Cutoff, scales = "free_y", ncol = 10, labeller = label_both) +
  
  # Add correlation text (Very small size to fit)
  geom_text(data = cor_data, aes(label = label_text, x = Inf, y = Inf), 
            hjust = 1.1, vjust = 1.2, size = 1.8, fontface = "bold", inherit.aes = FALSE) +
  
  theme_bw() +
  labs(title = "HFR vs Protein Length across 100 Cutoff Values (0.00 - 0.99)",
       x = "Length", y = "HFR") +
  
  # Theme adjustments for crowded plot
  theme(
    strip.text = element_text(size = 5),        # Smaller title for each subplot
    axis.text = element_text(size = 4),         # Smaller axis numbers
    axis.title = element_text(size = 8),
    panel.grid.major = element_blank(),         # Remove grid lines for clarity
    panel.grid.minor = element_blank()
  )

# Print plot (Might be messy in small window)
print(p)

# RECOMMENDED: Save as a large image file to zoom in
ggsave("~/Desktop/HFR_100_Cutoffs.pdf", p, width = 20, height = 20, limitsize = FALSE)
print("Plot saved to Desktop as HFR_100_Cutoffs.pdf")
