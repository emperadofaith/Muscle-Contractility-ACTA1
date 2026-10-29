# Download packages
install.packages(c("ggplot2", "rentrez"))

library(rentrez)
library(ggplot2)

# Set assignment details
gene_symbol <- "ACTA1"
transcript_accession <- "NM_001100.4"
clinvar_id <- "18292"
variant_protein <- "p.Val165Met"
variant_position <- 165

# Create output folders
for (folder in c("data", "results", "figures")) {
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
}

# STEP 2: Retrieve the human ACTA1 gene record

hits <- entrez_search(
  db = "gene",
  term = paste0(gene_symbol, "[Gene Name] AND txid9606[Organism]"),
  retmax = 20
)

if (length(hits$ids) == 0) {
  stop("ACTA1 gene was not found in NCBI Gene.")
}

gene_id <- NA_character_
gene_info <- NULL

for (id in hits$ids) {
  rec <- entrez_summary(db = "gene", id = id)
  
  if (toupper(rec$name) == gene_symbol) {
    gene_id <- id
    gene_info <- rec
    break
  }
}

if (is.na(gene_id)) {
  stop("Exact human ACTA1 gene record was not found.")
}

cat("Gene:", gene_symbol, "\n")
cat("NCBI Gene ID:", gene_id, "\n")
cat("Official name:", gene_info$description, "\n")
cat("Chromosome:", gene_info$chromosome, "\n")

writeLines(
  c(
    paste("Gene symbol:", gene_symbol),
    paste("NCBI Gene ID:", gene_id),
    paste("Official name:", gene_info$description),
    paste("Chromosome:", gene_info$chromosome),
    paste("Transcript:", transcript_accession),
    paste("ClinVar Variation ID:", clinvar_id)
  ),
  "results/gene_summary.txt"
)

# STEP 3: Retrieve the protein matching NM_001100.4

gb <- entrez_fetch(
  db = "nuccore",
  id = transcript_accession,
  rettype = "gb",
  retmode = "text"
)

# Extract the RefSeq protein accession
protein_pattern <- '/protein_id="NP_[0-9]+\\.[0-9]+"'
protein_match <- regmatches(
  gb,
  regexpr(protein_pattern, gb)
)

if (length(protein_match) == 0 || protein_match == "") {
  stop("Protein accession not found. Check the transcript record.")
}

protein_accession <- gsub(
  '^/protein_id="|"$',
  "",
  protein_match
)

print(protein_accession)

cat("Matched protein:", protein_accession, "\n")

# Download the matching protein sequence
protein_fasta <- entrez_fetch(
  db = "protein",
  id = protein_accession,
  rettype = "fasta",
  retmode = "text"
)

if (!startsWith(protein_fasta, ">")) {
  stop("Protein FASTA retrieval failed.")
}

fasta_file <- file.path(
  "data",
  paste0(gene_symbol, "_protein.fasta")
)

writeLines(protein_fasta, fasta_file)

writeLines(
  c(
    paste("Transcript:", transcript_accession),
    paste("Protein accession:", protein_accession)
  ),
  "results/protein_accession.txt"
)

# STEP 4: Calculate protein length and amino-acid composition

fasta_lines <- readLines(fasta_file)

sequence <- paste(
  fasta_lines[!grepl("^>", fasta_lines)],
  collapse = ""
)

sequence <- toupper(gsub("\\s+", "", sequence))
aa <- strsplit(sequence, "")[[1]]

protein_length <- length(aa)

# Count the 20 standard amino acids
standard_aa <- strsplit(
  "ACDEFGHIKLMNPQRSTVWY",
  ""
)[[1]]

aa_counts <- table(
  factor(aa, levels = standard_aa)
)

aa_results <- data.frame(
  AminoAcid = standard_aa,
  Count = as.integer(aa_counts)
)

write.csv(
  aa_results,
  "results/amino_acid_composition.csv",
  row.names = FALSE
)

cat("Protein length:", protein_length, "amino acids\n")
cat("Three most abundant amino acids:\n")

top_three <- head(
  aa_results[order(-aa_results$Count), ],
  3
)

print(top_three)

# STEP 5: Create the amino-acid composition graph

fig1 <- ggplot(
  aa_results,
  aes(x = AminoAcid, y = Count)
) +
  geom_col(fill = "grey45") +
  theme_bw() +
  labs(
    title = "ACTA1 Amino-Acid Composition",
    x = "Amino acid",
    y = "Number of residues"
  )

print(fig1)

ggsave(
  "figures/01_amino_acid_composition.png",
  plot = fig1,
  width = 8,
  height = 5,
  dpi = 300
)

# STEP 6: Check the assigned variant position

cat("Assigned variant:", variant_protein, "\n")
cat("Variant position:", variant_position, "\n")

# Check that the position exists in the protein
if (protein_length < variant_position) {
  stop("Variant position exceeds the retrieved protein length.")
}

# Correctly retrieve the amino acid at position 165
reference_residue <- aa[variant_position]

cat(
  "Reference residue at position 165:",
  reference_residue,
  "\n"
)

# Val165Met means valine (V) is replaced by methionine (M)
if (is.na(reference_residue)) {
  stop("Could not retrieve the amino acid at position 165.")
} else if (reference_residue != "V") {
  warning(
    "Position 165 is not valine. Verify the reference protein sequence."
  )
} else {
  cat("Reference residue verified: Valine (V)\n")
  cat("Variant change: Valine (V) to Methionine (M)\n")
}

# STEP 7: Plot the variant position on the reference protein

fig2 <- ggplot() +
  geom_segment(
    aes(
      x = 1,
      xend = protein_length,
      y = 0,
      yend = 0
    ),
    linewidth = 2
  ) +
  geom_point(
    aes(x = variant_position, y = 0),
    size = 4
  ) +
  annotate(
    "text",
    x = variant_position,
    y = 0.25,
    label = variant_protein,
    size = 4
  ) +
  scale_x_continuous(
    limits = c(0, protein_length),
    expand = expansion(mult = 0.02)
  ) +
  coord_cartesian(ylim = c(-0.6, 0.6)) +
  theme_void() +
  labs(
    title = "ACTA1 Reference Protein and Variant Position",
    x = "Amino-acid position"
  ) +
  theme(
    axis.title.x = element_text(),
    axis.text.x = element_text(),
    axis.ticks.x = element_line()
  )

print(fig2)

ggsave(
  "figures/02_variant_position.png",
  plot = fig2,
  width = 8,
  height = 2.5,
  dpi = 300
)

# STEP 8: Save the assigned ClinVar reference

clinvar_url <- paste0(
  "https://www.ncbi.nlm.nih.gov/clinvar/variation/",
  clinvar_id,
  "/"
)

writeLines(
  c(
    paste("Gene:", gene_symbol),
    paste("Transcript:", transcript_accession),
    paste("DNA variant: c.493G>A"),
    paste("Protein variant:", variant_protein),
    paste("ClinVar Variation ID:", clinvar_id),
    paste("ClinVar URL:", clinvar_url)
  ),
  "results/variant_summary.txt"
)

cat("\nAnalysis completed.\n")
cat("Check the data, results, and figures folders.\n")
