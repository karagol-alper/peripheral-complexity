# Load required packages
library(dplyr)
library(signal)  # for FFT
library(tibble)

# Define the function to compute HFR
compute_HFR <- function(rmsf_vector, cutoff_ratio = 0.99) {
  n <- length(rmsf_vector)
  fft_result <- fft(rmsf_vector)
  magnitude <- Mod(fft_result)[1:(n/2)]  # keep positive frequencies
  freqs <- seq(0, 1, length.out = n)[1:(n/2)]

  # Cutoff: retain high-frequency part above the given ratio
  cutoff_index <- floor(length(freqs) * cutoff_ratio)
  high_freq_sum <- sum(magnitude[(cutoff_index + 1):length(magnitude)])
  total_sum <- sum(magnitude)

  HFR <- high_freq_sum / total_sum
  return(HFR)
}

# Set folder path containing the .tsv files
folder_path <- "/Users/alperkaragol/Desktop/RMSF"  # replace with your path
file_list <- list.files(path = folder_path, pattern = "\\.tsv$", full.names = TRUE)

# Initialize result dataframe
results <- tibble(file = character(), HFR = numeric(), ProteinLength = numeric())

# Process each file
for (file in file_list) {
  df <- read.delim(file, header = TRUE)
  
  # Compute median RMSF per residue
  df$RMSF_median <- apply(df[, grepl("RMSF", names(df))], 1, median)

  # Calculate HFR
  rmsf_vector <- df$RMSF_median
  hfr_value <- compute_HFR(rmsf_vector)

  # Save results
  results <- results %>%
    add_row(file = basename(file), HFR = hfr_value, ProteinLength = length(rmsf_vector))
}

# Sort and inspect
results <- results %>% arrange(HFR)
print(results)

# Optional: plot HFR vs protein length
plot(results$ProteinLength, results$HFR, pch = 16,
     xlab = "Protein Length", ylab = "High-Frequency Ratio (HFR)",
     main = "Compactness vs Peripheral Flexibility (HFR)")
abline(lm(HFR ~ ProteinLength, data = results), col = "red")
