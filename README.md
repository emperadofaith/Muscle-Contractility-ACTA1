# Molecular Mechanism of Muscle-Contractility: ACTA1

## STUDENT INFORMATION ##
**Student name:** Faith Denielle Emperado

**Student ID:** 202303451

## ASSIGNMENT INFORMATION ##
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
How does ACTA1 contribute to muscle function, and what is the scientifically supported effect of its assigned ClinVar variant on the encoded protein and muscle activity?

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

# Figure 2. variant position
<img width="1278" height="398" alt="image" src="https://github.com/user-attachments/assets/52456bde-62c0-45f5-b99e-95efc81149d1" />

  The figure shows the position of the Val165Met (p.Val165Met) variant in the ACTA1 reference protein. The black line represents the protein's amino-acid sequence, while the black dot marks position 165, where valine (Val) is replaced by methionine (Met). The marker lies within the protein sequence, indicating that the mutation affects a specific amino-acid residue rather than the end of the protein. Since ACTA1 encodes skeletal muscle alpha-actin, this substitution may influence the protein's structure or interactions involved in muscle contraction.

Overall, the figures indicate that the ACTA1 protein contains varying amounts of the 20 standard amino acids and that the Val165Met variant occurs at a specific position within its sequence. Since ACTA1 encodes skeletal muscle alpha-actin, this substitution is relevant to studying the protein's role in muscle contraction. The amino-acid composition and variant position provide a basis for further investigation of how the substitution may influence protein structure and function.

# ClinVar RECORD OF NM_001100.4(ACTA1):c.493G>A 

### Table 1. ACTA1 Variant Annotation from ClinVar

| Field | Annotation |
|---|---|
| Gene | ACTA1 |
| Variant Name | NM_001100.4(ACTA1):c.493G>A (p.Val165Met) |
| Transcript Accession | NM_001100.4 |
| Coding Change | c.493G>A |
| Protein Change | p.Val165Met (V165M) |
| Variant Type | Single-nucleotide variant (SNV) |
| Sequence Effect | Missense |
| Associated Condition | Alpha-actinopathy |
| Clinical Significance | Pathogenic |
| Review Status | Reviewed by expert panel |
| ClinVar Variation ID | 18292 |
| Supporting Publications | PMID: 10508519, 15198992, 16427282, 17705262 |
| Date Accessed | October 29, 2026 |

The variant is a missense mutation, meaning that a single nucleotide substitution changes one amino acid in the protein. In ACTA1, the coding DNA change from guanine (G) to adenine (A) at position c.493 changes the codon so that valine (Val) is replaced by methionine (Met) at amino-acid position 165. 
NCBI+1

ACTA1 encodes skeletal muscle alpha-actin, a protein important for muscle contraction. ClinVar reports experimental evidence that this variant affects ACTA1 function. The expert-panel assessment also describes evidence of altered actin behavior and protein aggregation in experimental models. The variant is classified as pathogenic for alpha-actinopathy.

# DISCUSSION

The ACTA1 gene encodes skeletal muscle alpha-actin, a protein that forms the thin filaments of muscle fibers. It interacts with myosin to produce muscle contraction and helps maintain the structure and stability of sarcomeres, the basic contractile units of muscle fibers. Calcium signaling regulates the interaction between actin and myosin, while alpha-actin provides an essential structural component for effective contraction. Therefore, changes in alpha-actin may affect muscle contraction and the structural integrity of muscle fibers.

The assigned variant, NM_001100.4:c.493G>A (p.Val165Met), is a missense variant in which valine is replaced by methionine at amino acid position 165. ClinVar classifies this variant as pathogenic for alpha-actinopathy. Experimental studies in C2C12 muscle cells observed actin aggregates inside the nucleus and smaller amounts in the cytoplasm. Studies using fruit flies also found abnormalities in Z-disc assembly and sarcomere stability, suggesting that the variant can disrupt normal muscle structure. Computational analysis using REVEL produced a score of 0.961, predicting a potentially damaging effect on protein function. However, this score is a prediction rather than direct experimental proof. The variant has also been reported in individuals with muscle disorders, providing clinical evidence of its association with disease. Although these findings support its pathogenicity, the exact mechanism by which the substitution affects muscle function is not fully established. 

**The proposed molecular mechanism can be summarized as follows:**

**ACTA1 gene → skeletal muscle alpha-actin → thin filament formation and sarcomere stability → p.Val165Met substitution → abnormal actin aggregation and disrupted sarcomere organization → possible impairment of muscle contraction and structural integrity.**

*This chain connects the genetic change to the observed structural abnormalities and their possible effects on muscle function. However, a direct effect on calcium signaling or muscle relaxation has not been established by the evidence discussed.*

To test this interpretation, normal and p.Val165Met ACTA1-expressing muscle cells could be examined using microscopy. The distribution of actin and the organization of muscle structures could be compared between the two groups by recording the presence and location of actin aggregates. If the variant-expressing cells show more aggregates or disrupted actin organization, this would support the idea that p.Val165Met affects protein organization and muscle structure. Additional experiments, such as measuring muscle force, would be needed to determine whether these structural changes directly impair contraction.

# LIMITATION
* Cell model limitations: Results from cultured muscle cells may not fully represent how the variant affects actual human muscle tissue.
* Uncertain cause and effect: Observing actin aggregates does not prove that these aggregates directly cause muscle weakness. ACTA1 variants may affect muscle function through different mechanisms.

# Possible Validation
The experiment could be validated by comparing normal ACTA1-expressing muscle cells with cells carrying the p.Val165Met variant under the same conditions. Microscopy images should be collected from multiple independent experiments, and actin aggregation and sarcomere organization should be measured using consistent methods. Muscle force measurement could also be added to determine whether structural abnormalities are associated with reduced contraction. If the variant-expressing cells consistently show abnormal actin organization and lower force production than the control cells, this would provide stronger evidence that the variant affects muscle structure and function.

# QUESTIONS TO ANSWER
**1. What does your assigned gene encode, and how does this protein contribute to muscle function?**

The ACTA1 gene encodes skeletal muscle alpha-actin, a protein that forms the thin filaments in muscle fibers. It interacts with myosin to produce muscle contraction and helps maintain sarcomere structure and stability. Calcium regulates the actin–myosin interaction, while alpha-actin provides an essential part of the contractile machinery.

**2. What was the exact protein accession and length retrieved with R?**
* Protein accession: NP_001091.1
* Protein length: 377 amino acids.

**3. What did the amino-acid composition chart show? What can it NOT tell you?**

The chart showed the frequency of each of the 20 amino acids in the ACTA1 protein. Isoleucine (I) was the most abundant, with 30 residues, while tryptophan (W) was the least abundant, with 4 residues. However, the chart cannot show the protein's three-dimensional structure, directly measure its activity, or determine how the assigned variant affects muscle function.
4. What is your fixed ClinVar variant, and what is its current classification and evidence strength?
The assigned variant is NM_001100.4:c.493G>A (NP_001091.1:p.Val165Met), also known as V165M or, under older numbering, V163M. It is a missense variant that replaces valine with methionine at amino acid position 165. ClinVar classifies it as Pathogenic for autosomal dominant alpha-actinopathy, with review by an expert panel. Supporting evidence includes experimental studies showing altered ACTA1 function, reports of the variant in affected individuals, segregation with disease in a family, its absence from the cited population database, and computational predictions. The REVEL score is 0.961, supporting a potentially damaging effect, although this is computational evidence rather than direct experimental proof. 

**5. What is the scientifically supported mechanism linking this variant to muscle function, and what remains uncertain?**

Experimental studies found that the variant was associated with actin aggregates in the nucleus and, at lower levels, in the cytoplasm of C2C12 muscle cells. Studies in fruit flies also showed disrupted Z-disc assembly and reduced sarcomere stability. These findings support a mechanism in which altered actin organization affects muscle structure and may interfere with normal contraction. However, the exact molecular process causing these changes is not fully established, and a direct effect on calcium signaling or muscle relaxation has not been demonstrated by the evidence discussed here. 

**6. What additional experiment would help test the proposed mechanism?**

Microscopic examination combined with muscle-force measurement could help test the proposed mechanism. Normal and p.Val165Met ACTA1-expressing muscle cells could be compared for actin aggregation and sarcomere organization using microscopy. Their contractile force could then be measured to determine whether structural abnormalities are associated with reduced contraction. Consistent structural abnormalities and lower force in variant-expressing cells would provide further support for the proposed mechanism.

# DATABASE AND REFERENCES
1. National Center for Biotechnology Information (NCBI).
   ACTA1 Gene. Gene ID: 58.
   
   https://www.ncbi.nlm.nih.gov/gene/58

3. National Center for Biotechnology Information (NCBI).
   ACTA1 Reference Transcript. NM_001100.4.
   
   https://www.ncbi.nlm.nih.gov/nuccore/NM_001100.4

5. National Center for Biotechnology Information (NCBI).
   ACTA1 Reference Protein. NP_001091.1.

   https://www.ncbi.nlm.nih.gov/protein/NP_001091.1

7. National Center for Biotechnology Information (NCBI).
   ClinVar: ACTA1 c.493G>A (p.Val165Met).
   Variation ID: 18292.
   
   https://www.ncbi.nlm.nih.gov/clinvar/variation/18292/

**Supoorting Publications**
   
Cooper, S. T., et al. (2004). Evidence for a dominant-negative effect in ACTA1 nemaline myopathy caused by abnormal folding, aggregation and altered polymerization of mutant actin isoforms. Human Molecular Genetics, 13(16), 1727–1743.

https://doi.org/10.1093/hmg/ddh185

Clinical Genome Resource (ClinGen), Congenital Myopathies Variant Curation Expert Panel. ACTA1 c.493G>A (p.Val165Met) variant interpretation.

https://erepo.clinicalgenome.org/evrepo/ui/interpretation/dafe07a8-2578-4074-8f38-bad29ef7c164
