# Load necessary libraries
library(org.Hs.eg.db)
library(AnnotationDbi)

# Define the file paths
gmt_path <- "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/KEGG_2022_entrez.gmt"
output_path <- "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/IDs/entrez2ensembl_kegg.txt"

# Function to extract all Entrez gene IDs from the GMT file
extract_entrez_ids <- function(gmt_path) {
  lines <- readLines(gmt_path)
  entrez_ids <- unique(unlist(lapply(lines, function(line) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) > 2) {
      return(parts[-(1:2)])  # Remove pathway name and description, keep gene names
    } else {
      return(NULL)
    }
  })))
  return(entrez_ids)
}

# Extract Entrez IDs from the GMT file
entrez_ids <- extract_entrez_ids(gmt_path)

# Debugging step: Check if the Entrez IDs were extracted correctly
cat("Number of Entrez IDs extracted:", length(entrez_ids), "\n")

# Convert Entrez IDs to Ensembl IDs using org.Hs.eg.db
ensembl_ids <- mapIds(org.Hs.eg.db,
                      keys = entrez_ids,
                      column = "ENSEMBL",
                      keytype = "ENTREZID",
                      multiVals = "first")

# Create a mapping data frame from the mapIds result
mapping_df <- data.frame(Entrez_ID = names(ensembl_ids),
                         Ensembl_ID = as.character(ensembl_ids),
                         stringsAsFactors = FALSE)

# Manually add missing mappings as a data frame
manual_mappings <- data.frame(
  Entrez_ID = c('4578', '1564', '4574', '4511', '4576', '4558', '4555', '4577', '4570', '4553',
                '102724788', '4566', '4550', '4563', '4568', '780851', '4579', '4575', '4573',
                '392133', '6315', '26851', '4572', '4567', '4571', '4565', '4556', '780853', 
                '4569', '4549', '107987479', '4564', '9103', '780852'),
  Ensembl_ID = c('ENSG00000210117', 'ENSG00000205702', 'ENSG00000210151', 'ENSG00000210140', 
                 'ENSG00000210195', 'ENSG00000210049', 'ENSG00000210154', 'ENSG00000210077',
                 'ENSG00000210135', 'ENSG00000210127', 'ENSG00000100033', 'ENSG00000210156',
                 'ENSG00000210082', 'ENSG00000210164', 'ENSG00000210191', 'ENSG00000263934',
                 'ENSG00000210144', 'ENSG00000210184', 'ENSG00000210174', 'ENSG00000176510',
                 'ENSG00000230223', 'ENSG00000265185', 'ENSG00000210107', 'ENSG00000209082',
                 'ENSG00000210196', 'ENSG00000210100', 'ENSG00000210194', 'ENSG00000264940',
                 'ENSG00000210112', 'ENSG00000211459', 'ENSG00000100197', 'ENSG00000210176',
                 'ENSG00000244682', 'ENSG00000262074'),
  stringsAsFactors = FALSE
)


# Append the manual mappings to the original mapping dataframe
mapping_df <- rbind(mapping_df, manual_mappings)

# Optional: Remove rows where Ensembl_ID is NA
mapping_df <- mapping_df[!is.na(mapping_df$Ensembl_ID), ]

# Save the mapping to a TSV file
write.table(mapping_df, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)

# Output the number of mapped and unmapped Entrez IDs
cat("Total Entrez IDs:", length(entrez_ids), "\n")
cat("Mapped Entrez IDs:", nrow(mapping_df), "\n")
cat("Unmapped Entrez IDs:", length(entrez_ids) - nrow(mapping_df), "\n")
