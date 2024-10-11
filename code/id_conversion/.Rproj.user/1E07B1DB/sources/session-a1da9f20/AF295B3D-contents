# Ensure necessary libraries are loaded
library(graph)
library(graphite)
library(Rgraphviz)

# Load KEGG pathways for humans
kegg <- pathways("hsapiens", "kegg")

# Convert the pathway identifiers to ENTREZID format
kegg <- convertIdentifiers(kegg, "ENTREZID")

# Retrieve the "Dilated cardiomyopathy" pathway
p <- kegg[["Dilated cardiomyopathy"]]

# Create a graph from the pathway
graph <- pathwayGraph(p)

# Define the file path where you want to save the plot
output_file <- "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/dilated_cardiomyopathy_graph.png"

# Save the plot as a PNG file with high resolution (e.g., 300 dpi)
png(output_file, width = 2000, height = 2000, res = 300)

# Plot the graph with labels
plot(graph, main = pathwayTitle(p), edge.arrow.size = 0.5)  # Adjust edge.arrow.size as needed

# Turn off the device to save the plot
dev.off()

# Print message to confirm save
cat("Graph saved successfully to:\n", output_file, "\n")

# Save p@protEdges as a TSV file
prot_edges <- p@protEdges

# Convert to a data frame
prot_edges_df <- as.data.frame(prot_edges)

# Define the output file path for the TSV file
output_prot_edges_file <- "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/dilated_cardiomyopathy_prot_edges.tsv"

# Save the data frame as a TSV file
write.table(prot_edges_df, output_prot_edges_file, sep = "\t", row.names = FALSE, quote = FALSE)

# Print message to confirm save
cat("Protein edges saved successfully to:\n", output_prot_edges_file, "\n")

