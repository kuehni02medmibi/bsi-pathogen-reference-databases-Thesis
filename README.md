# BSI Pathogen Reference Databases

Custom Minimap2 and Kraken2 reference databases for bloodstream infection (BSI) pathogen detection using Oxford Nanopore sequencing.

---

## Overview

This repository contains the reference database resources used for taxonomic classification of microbial cell-free DNA (mcfDNA) in bloodstream infection (BSI) studies using Oxford Nanopore sequencing.

A curated pathogen panel was developed to include the most clinically relevant Gram-positive and Gram-negative bacterial species associated with bloodstream infections while simultaneously accounting for host-derived DNA through inclusion of the human reference genome (GRCh38).

The same reference sequence collection was used to generate both the Minimap2 and Kraken2 databases. This ensures that any differences in taxonomic classification performance arise from the underlying classification algorithms rather than differences in database composition.

---

## Pathogen Panel

### Gram-positive bacteria

- Staphylococcus aureus
- Staphylococcus epidermidis
- Staphylococcus haemolyticus
- Staphylococcus hominis
- Staphylococcus capitis
- Staphylococcus lugdunensis
- Enterococcus faecalis
- Enterococcus faecium
- Streptococcus pneumoniae
- Streptococcus pyogenes
- Streptococcus agalactiae
- Streptococcus dysgalactiae

### Gram-negative bacteria

- Escherichia coli
- Klebsiella pneumoniae
- Klebsiella oxytoca
- Enterobacter cloacae
- Citrobacter freundii
- Serratia marcescens
- Proteus mirabilis
- Morganella morganii
- Pseudomonas aeruginosa
- Acinetobacter baumannii
- Aggregatibacter aphrophilus

### Host Reference

- Homo sapiens (GRCh38)

---

## Database Construction Strategy

### Reference Selection

For each bacterial species, a representative genome assembly was obtained from the NCBI RefSeq database using the following criteria:

1. RefSeq reference genome preferred
2. Representative genome if no reference genome was available
3. Complete genome assembly preferred
4. Latest assembly version selected

The human reference genome GRCh38 was included to represent host-derived cfDNA background.

### Taxonomic Annotation

Each reference sequence was assigned its corresponding NCBI TaxID.

Taxonomic resources were derived from the NCBI Taxonomy database and include:

- names.dmp
- nodes.dmp
- merged.dmp
- delnodes.dmp

A reference-to-taxonomy mapping file (`ref2taxid.tsv`) was generated to ensure unambiguous taxonomic assignment.

---

## Database Types

### Minimap2 Database

Components:

- Combined reference FASTA
- Minimap2 index (`.mmi`)
- Taxonomy files
- `ref2taxid.tsv`

Purpose:

- Alignment-based classification
- Genome coverage analysis
- Read-length analysis
- Mapping statistics

### Kraken2 Database

Components:

- Kraken2 hash tables (`hash.k2d`)
- Taxonomy database (`taxo.k2d`)
- Database options (`opts.k2d`)

Purpose:

- Exact k-mer-based taxonomic classification
- Rapid taxonomic profiling
- Classifier comparison experiments

---

## Repository Contents

```text
selected_assemblies_v2.tsv
build_bsi_v2.py
taxonomy/
ref2taxid.tsv
README.md
```

### Important Files

| File | Description |
|--------|--------|
| selected_assemblies_v2.tsv | Complete list of reference assemblies used |
| build_bsi_v2.py | Database generation workflow |
| taxonomy/ | NCBI taxonomy resources |
| ref2taxid.tsv | Sequence-to-TaxID mapping |

---

## Reproducibility

Both Minimap2 and Kraken2 databases were generated from the same curated reference collection.

This design allows direct comparison of classifier performance while controlling for differences in database composition.

---

## Citation

If you use these databases, please cite the associated thesis.

DOI: *to be added after Zenodo deposition*

---

## Author

Niklas Kühn

Master's Thesis – Novel Approaches for Growth-independent Detection of Bloodstream Pathogens: From Culture-associated Sequencing to Cell-free DNA Analysis


