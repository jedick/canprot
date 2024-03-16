# canprot/extdata/aa/prepare.R
# Prepare amino amino acid compositions for ribosomal proteins
# in Nitrososphaeria genomes analyzed by Luo et al., 2024
# doi:10.1093/ismejo/wrad031
# 20240307

# Read Table S1 (existing) and S3 (new) genome data
# Note: These were exported from the SI spreadsheets
# and columns 1-2 were renamed Order and Family
d1 <- read.csv("Table_S1.csv")
d3 <- read.csv("Table_S3.csv")

# For each table, keep the genomes with both habitat and respiration data
d1 <- d1[d1$Respiration.typed != "-" & d1$Habitate.typee != "-", ]
d3 <- d3[d3$Respiration.typed != "-" & d3$Habitate.typee != "-", ]

# Extract columns and normalize names
d1 <- d1[, c(1:21)]
d3 <- d3[, c(1:3, 5, 6:12, 14:23)]
newnames <- c("Order", "Family", "Accession", "Organism", "Source description", "Scaffolds", "Genome size (kb)",
"Completeness", "Contamination", "Heterogeneity", "Estimated size", "Rep", "tRNA",
"5S rRNA", "23S rRNA", "16S rRNA", "Predicted OGT (°C)",
"Respiration type", "Habitat type", "Sample pH", "Sample Temp (°C)"
)
colnames(d1) <- newnames
colnames(d3) <- newnames

# Merge the tables
nitrososphaeria_MAGs <- rbind(d1, d3)
# Remove trailing space
nitrososphaeria_MAGs$"Source description" <- gsub("\ $", "", nitrososphaeria_MAGs$"Source description")
# Remove trailing character (this not an ASCII space)
nitrososphaeria_MAGs$Accession <- gsub(" $", "", nitrososphaeria_MAGs$Accession)
write.csv(nitrososphaeria_MAGs, "nitrososphaeria_MAGs.csv", quote = c(4, 5), row.names = FALSE)

# List unique combinations of respiration and habitat types
sort(table(paste(nitrososphaeria_MAGs$Respiration, nitrososphaeria_MAGs$Habitat)), decreasing = TRUE)
# Aerobic Nonthermal: 133
# Anaerobic Nonthermal: 24
# Anaerobic Thermal: 15
# Aerobic Thermal: 10
# Microaerobic Nonthermal: 2
# Microaerobic Thermal: 2

# Download genomes
for(accession in nitrososphaeria_MAGs$Accession) {
  outfile <- file.path(".", paste0(accession, ".zip"))
  if(!file.exists(outfile)) {
    print(outfile)
    # 'datasets' command from NCBI's command line tools
    cmd <- paste("datasets download genome accession", accession, "--include gff3,rna,cds,protein,genome,seq-report")
    system(cmd)
    file.rename("ncbi_dataset.zip", outfile)
  }
}

# Get protein FASTA
aalist <- lapply(nitrososphaeria_MAGs$Accession, function(accession) {
  infile <- file.path(".", paste0(accession, ".zip"))
  print(infile)
  cmd <- paste("unzip", infile)
  system(cmd)
  faafile <- file.path("ncbi_dataset/data/", accession, "protein.faa")
  outfile <- paste0(accession, ".faa")
  myaa <- NULL
  if(file.exists(faafile)) {
    myaa <- read_fasta(faafile)
    myaa <- sum_aa(myaa)
    myaa$protein <- "sum"
    myaa$organism <- accession
  }
  system("rm -rf ncbi_dataset/ README.md")
  myaa
})
aa <- do.call(rbind, aalist)
write.csv(aa, "nitrososphaeria_aa.csv", row.names = FALSE, quote = FALSE)
