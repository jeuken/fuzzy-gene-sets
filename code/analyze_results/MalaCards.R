# Load necessary libraries
library(readr)
library(AnnotationDbi)
library(org.Hs.eg.db)  # Adjust if using a different organism
library(dplyr)  # For data manipulation

# Path to your CSV file
file_path <- '/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/MalaCards - Genes associated with Alzheimers Disease.csv'

# Read the CSV file, specifying the separator and skipping the header rows
df <- read_delim(file_path, delim = ",", skip = 3)

# Extract only the Symbol and Score columns
df_filtered <- df[c("Symbol", "Score")]

# Convert gene symbols to Ensembl IDs
df_filtered$Ensembl_ID <- mapIds(org.Hs.eg.db, 
                                 keys = df_filtered$Symbol, 
                                 column = "ENSEMBL", 
                                 keytype = "SYMBOL", 
                                 multiVals = "first")

# Filter out rows with NA values in the Ensembl_ID column
df_filtered <- na.omit(df_filtered)  # Option 1: Using na.omit()

# Alternatively, using dplyr:
# df_filtered <- df_filtered %>% filter(!is.na(Ensembl_ID))  # Option 2: Using dplyr

# Check the result
print(df_filtered)

# Save the resulting DataFrame to a new CSV file
output_file_path <- '/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/MalaCards - Genes associated with Alzheimers Disease - Ensembl.csv'
write.csv(df_filtered, output_file_path, row.names = FALSE)
