# Load the necessary libraries
library(org.Hs.eg.db)
library(AnnotationDbi)

# Define the function to map IDs and save the results
gene_id_mapping <- function(input_id_type, output_id_type, output_dir) {
  # Retrieve the mappings from the specified input ID type to the output ID type
  id_map <- select(org.Hs.eg.db, 
                   keys = keys(org.Hs.eg.db, keytype = input_id_type), 
                   columns = c(output_id_type), 
                   keytype = input_id_type)
  
  # Remove rows with NA values (if any)
  id_map <- na.omit(id_map)
  
  # Rename the columns
  colnames(id_map) <- c(input_id_type, output_id_type)
  
  # View the resulting table (optional)
  head(id_map)
  
  # Save the table as a txt file with tab-separated values in the specified directory
  write.table(id_map, 
              file.path(output_dir, paste0(input_id_type, "_to_", output_id_type, ".txt")), 
              sep = "\t", 
              row.names = FALSE, 
              quote = FALSE)
  
  return(id_map)  # Return the complete mapping
}

# Example usage of the function
gene_id_mapping("ENTREZID", "ENSEMBL", "../../data/IDs")
gene_id_mapping("SYMBOL", "ENSEMBL", "../../data/IDs")

