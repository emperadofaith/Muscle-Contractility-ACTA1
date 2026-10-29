# Molecular Mechanism of Muscle-Contractility: ACTA1

## STUDENT INFORMATION ##
**Student name:** Faith Denielle Emperado

**Student ID:** 202303451

## PROJECT INFORMATION ##
**Assignment number:** 01

**Gene:** ACTA1

**NCBI Gene ID:** 58 

**Assigned transcript:** NM_001100.4

**Protein accession:** NP_001091.1

**ClinVar Variation ID:** 18292

**Variant type:** missense

**Assigned variant** :c.493G>A (p.Val165Met)

**Clinvar condition:** pathogenic

## RESEARCH QUESTION


## METHODS

### 1. R Project Setup
A project folder named `Muscle_Contractility_ACTA1` was created in RStudio. An R script named `analysis.R` and three subfolders (`data`, `results`, and `figures`) were prepared to organize the code and output files. The `rentrez` and `ggplot2` packages were installed and loaded.

```r
install.packages(c("rentrez", "ggplot2"))

library(rentrez)
library(ggplot2)

gene_symbol <- "ACTA1"
transcript_accession <- "NM_001100.4"
clinvar_id <- "18292"
variant_protein <- "p.Val165Met"

for (folder in c("data", "results", "figures")) {
  dir.create(folder, showWarnings = FALSE)
}
```

### 2. Retrieval of the ACTA1 Gene Record
The NCBI Gene database was searched using the `rentrez` package to retrieve the human ACTA1 gene record. The exact gene symbol was checked, and the NCBI Gene ID, official gene name, and chromosome information were recorded in `results/gene_summary.txt`.

```r
hits <- entrez_search(
  db = "gene",
  term = paste0(
    gene_symbol,
    "[Gene Name] AND txid9606[Organism]"
  ),
  retmax = 20
)

gene_info <- entrez_summary(db = "gene", id = "58")
```

### 3. Retrieval of the Reference Protein Sequence
The assigned RefSeq transcript, `NM_001100.4`, was retrieved from NCBI Nucleotide. Its GenBank record was examined to identify the corresponding RefSeq protein accession. The protein sequence was then downloaded in FASTA format and saved as `data/ACTA1_protein.fasta`. The protein accession was recorded in `results/protein_accession.txt`.

```r
gb <- entrez_fetch(
  db = "nuccore",
  id = transcript_accession,
  rettype = "gb",
  retmode = "text"
)

protein_fasta <- entrez_fetch(
  db = "protein",
  id = protein_accession,
  rettype = "fasta",
  retmode = "text"
)

writeLines(
  protein_fasta,
  "data/ACTA1_protein.fasta"
)
```

### 4. Analysis of Protein Length and Amino-Acid Composition
The FASTA file was read in R, and the amino-acid sequence was extracted. The total protein length was calculated, and the frequencies of the 20 standard amino acids were counted. The results were saved as `results/amino_acid_composition.csv`, and the three most abundant amino acids were identified.

```r
fasta_lines <- readLines(
  "data/ACTA1_protein.fasta"
)

sequence <- paste(
  fasta_lines[!grepl("^>", fasta_lines)],
  collapse = ""
)

aa <- strsplit(toupper(sequence), "")[[1]]
protein_length <- length(aa)

standard <- strsplit(
  "ACDEFGHIKLMNPQRSTVWY", ""
)[[1]]

counts <- table(factor(aa, levels = standard))

results <- data.frame(
  AminoAcid = standard,
  Count = as.integer(counts)
)

write.csv(
  results,
  "results/amino_acid_composition.csv",
  row.names = FALSE
)

head(results[order(-results$Count), ], 3)
```

### 5. Visualization of Amino-Acid Composition
A bar graph was generated using `ggplot2` to show the number of occurrences of each standard amino acid in the ACTA1 reference protein. The graph was saved as `figures/01_amino_acid_composition.png` at 300 dpi.

```r
fig1 <- ggplot(
  results,
  aes(x = AminoAcid, y = Count)
) +
  geom_col(fill = "grey45") +
  theme_bw() +
  labs(
    title = "ACTA1 Amino-Acid Composition",
    x = "Amino acid",
    y = "Number of residues"
  )

ggsave(
  "figures/01_amino_acid_composition.png",
  fig1,
  width = 8,
  height = 5,
  dpi = 300
)
```

### 6. Analysis of the Assigned ClinVar Variant
The assigned variant, `NM_001100.4:c.493G>A (p.Val165Met)`, was examined using its ClinVar Variation ID, 18292. The ClinVar record was used to document the variant description, associated condition, current classification, review status, and supporting evidence. The missense change indicates a substitution of valine with methionine at amino-acid position 165. The database classification and available evidence were recorded without assuming that the predicted molecular effect was experimentally confirmed.

[NCBI ClinVar: Variation 18292](https://www.ncbi.nlm.nih.gov/clinvar/variation/18292/)

### 7. Visualization of the Variant Position
Option A from the laboratory manual was selected because the assigned variant is a missense variant with a defined protein position. A linear protein map was generated to show the total reference-protein length and the location of the Val165Met substitution. The figure was saved as `figures/02_variant_position.png`.

```r
variant_position <- 165

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
    label = "p.Val165Met",
    size = 4
  ) +
  theme_void() +
  labs(
    title = "ACTA1 Reference Protein and Variant Position",
    x = "Amino-acid position"
  )

ggsave(
  "figures/02_variant_position.png",
  fig2,
  width = 8,
  height = 2.5,
  dpi = 300
)
```

### 8. Interpretation and Documentation
The retrieved sequence, amino-acid composition, and variant-position plot were used to describe the ACTA1 protein and the location of the assigned variant. The normal role of ACTA1 in skeletal muscle, the available evidence for the variant, and its possible effects on muscle function were discussed. Experimental findings, clinical associations, and computational predictions were distinguished. The analysis files, figures, database links, access dates, references, limitations, and a proposed validation experiment were documented in the GitHub repository.

## RESULTS

**ACTA1 reference protein**
  * Transcript Accession: NM_001100.4
  * Retrieved Protein Accession: NP_001091.1 (Length: 377 amino acids)

# Figure 1. amino acid composition
<img width="2400" height="1500" alt="image" src="https://github.com/user-attachments/assets/b6dc27c4-8748-41e8-af80-2a7cd545c2db" />

  The graph illustrates the amino-acid composition of the ACTA1 protein, showing that each of the 20 standard amino acids occurs in different amounts. Isoleucine (I) has the highest count, with 30 residues, followed by alanine (A) with 29, while tryptophan (W) has the lowest count, with only 4 residues. These differences reflect the protein’s amino-acid composition, which contributes to its chemical properties, folding, and overall structure. However, the graph alone cannot establish the protein’s function or determine how the assigned Val165Met mutation affects its structure and activity.

# Variant interpretation

# Figure 2. variant position
<img width="1278" height="398" alt="image" src="https://github.com/user-attachments/assets/52456bde-62c0-45f5-b99e-95efc81149d1" />

  The figure shows the position of the Val165Met (p.Val165Met) variant in the ACTA1 reference protein. The black line represents the protein's amino-acid sequence, while the black dot marks position 165, where valine (Val) is replaced by methionine (Met). The marker lies within the protein sequence, indicating that the mutation affects a specific amino-acid residue rather than the end of the protein. Since ACTA1 encodes skeletal muscle alpha-actin, this substitution may influence the protein's structure or interactions involved in muscle contraction.

Overall, the figures indicate that the ACTA1 protein contains varying amounts of the 20 standard amino acids and that the Val165Met variant occurs at a specific position within its sequence. Since ACTA1 encodes skeletal muscle alpha-actin, this substitution is relevant to studying the protein's role in muscle contraction. The amino-acid composition and variant position provide a basis for further investigation of how the substitution may influence protein structure and function.

# Variant inpterpretation


