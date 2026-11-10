# Bioinformatics Algorithms on Human Gene Sequences

---

## Overview

This assignment applies three fundamental algorithmic paradigms — **basic sequence analysis**, **greedy algorithms**, and **divide & conquer** — to real human gene sequences retrieved from NCBI. All sequences are stored locally as FASTA files and processed entirely in pure Python (no external bioinformatics libraries).

---

## Project Structure

```
Data_Scicence_for_Bioinformatics/
│
├── q1_insulin_analysis.py      # Question 1 – DNA sequence analysis of human INS gene
├── q2_greedy_motif.py          # Question 2 – Greedy Motif Search across 4 disease genes
├── q3_divide_and_conquer.py    # Question 3 – Hirschberg's D&C global sequence alignment
│
├── INS_sequence.fasta          # Human Insulin gene (NM_000207.3) — used in Q1
├── HBB.fasta                   # Hemoglobin subunit beta (NM_000518.5) — used in Q2
├── CFTR.fasta                  # Cystic Fibrosis Transmembrane Regulator — used in Q2
├── HTT.fasta                   # Huntingtin gene — used in Q2
├── TP53.fasta                  # Tumour Protein p53 — used in Q2
├── BRCA1.fasta                 # Breast Cancer gene 1 — used in Q3
└── BRCA2.fasta                 # Breast Cancer gene 2 — used in Q3
```

---

## Questions & Algorithms

### Question 1 – Human Insulin (INS) Gene Analysis

**File:** `q1_insulin_analysis.py`  
**Input:** `INS_sequence.fasta` (NM_000207.3, 446 bp mRNA)

Performs basic DNA sequence analysis on the human insulin gene:

| Sub-part | Task |
|----------|------|
| (a) | Read and parse FASTA file; extract header & sequence |
| (b) | Count total sequence length in base pairs |
| (c) | Count individual nucleotide frequencies (A, T, G, C) |
| (d) | Calculate **GC content** (%) |
| (e) | Generate the **5′→3′ reverse complement** |

**Key functions:**
- `read_fasta(filepath)` — Parses FASTA format, strips headers, upper-cases sequence
- `calculate_gc_content(seq)` — Computes G, C counts and GC%
- `get_reverse_complement(seq)` — Builds complement using a lookup map, then reverses

---

### Question 2 – Greedy Motif Search

**File:** `q2_greedy_motif.py`  
**Input:** First 500 bp of `HBB.fasta`, `CFTR.fasta`, `HTT.fasta`, `TP53.fasta`

Implements the **Greedy Motif Search** algorithm to find a conserved 10-mer motif across four clinically significant human disease genes.

**Algorithm steps:**
1. Initialise best motifs with the first k-mer of each sequence
2. Iterate through all k-mers in the first sequence
3. For each candidate k-mer, greedily extend to subsequent sequences by selecting the **profile-most-probable k-mer** using a **4×k profile matrix** (built with Laplace pseudocounts to avoid zero probabilities)
4. Score the motif set by summing mismatches from the consensus column-by-column; keep the lowest-scoring set

**Key functions:**
- `build_profile_matrix(motifs, k)` — Constructs probability matrix with Laplace's Rule of Succession (+1 pseudocount)
- `profile_most_probable_kmer(dna_seq, k, profile)` — Finds highest-probability k-mer under the current profile
- `score_motifs(motifs, k)` — Derives consensus string and mismatch score
- `greedy_motif_search(dna_sequences, k)` — Full greedy search returning best motifs, consensus, and score

**Output includes:** Best motif per gene, consensus motif, mismatch score, and positional conservation percentage.

---

### Question 3 – Divide & Conquer Global Sequence Alignment (Hirschberg's Algorithm)

**File:** `q3_divide_and_conquer.py`  
**Input:** First 1000 bp of `BRCA1.fasta` and `BRCA2.fasta`

Implements **Hirschberg's Algorithm** — a space-efficient divide-and-conquer approach to global pairwise sequence alignment. Combines the correctness of Needleman-Wunsch with **O(n) space complexity** instead of O(mn).

**Scoring scheme:**

| Event | Score |
|-------|-------|
| Match | +2 |
| Mismatch | −1 |
| Gap (indel) | −2 |

**Algorithm steps:**
1. **Base case** (`len(X) ≤ 1` or `len(Y) ≤ 1`): Run standard Needleman-Wunsch with traceback
2. **Divide**: Split sequence X at midpoint
3. **Conquer**: Compute forward DP scores for `X[:mid]` vs `Y` and backward DP scores for the reversed second half
4. **Combine**: Find split point `k*` in Y that maximises `forward_score[k] + backward_score[k]`
5. **Recurse**: Align left halves and right halves independently and concatenate

**Key functions:**
- `linear_space_score(X, Y)` — O(|Y|) memory DP that returns only the final score row
- `base_case_align(X, Y)` — Full Needleman-Wunsch with traceback for base cases
- `divide_and_conquer_align(X, Y)` — Recursive Hirschberg alignment

**Output includes:** Alignment length, similarity pairs, percent identity, mismatches, gaps, and a 75-column alignment preview.

---

## How to Run

> **Requirements:** Python 3.7+ — no third-party packages needed.

All scripts must be run from the project directory (where the `.fasta` files reside):

```bash
cd Data_Scicence_for_Bioinformatics

# Question 1 – INS gene analysis
python q1_insulin_analysis.py

# Question 2 – Greedy motif search
python q2_greedy_motif.py

# Question 3 – Divide & conquer alignment
python q3_divide_and_conquer.py
```

---

## Gene Sequences Used

| Gene | Full Name | NCBI Accession | Clinical Relevance |
|------|-----------|----------------|--------------------|
| INS | Insulin | NM_000207.3 | Type 1 & Type 2 Diabetes |
| HBB | Hemoglobin Subunit Beta | NM_000518.5 | Sickle Cell Anaemia, Thalassaemia |
| CFTR | CF Transmembrane Regulator | — | Cystic Fibrosis |
| HTT | Huntingtin | — | Huntington's Disease |
| TP53 | Tumour Protein p53 | — | Multiple Cancers |
| BRCA1 | Breast Cancer Gene 1 | — | Hereditary Breast & Ovarian Cancer |
| BRCA2 | Breast Cancer Gene 2 | — | Hereditary Breast & Ovarian Cancer |

---

## Algorithmic Complexity

| Question | Algorithm | Time Complexity | Space Complexity |
|----------|-----------|-----------------|------------------|
| Q1 | Linear scan | O(n) | O(n) |
| Q2 | Greedy Motif Search | O(t · n · k · t) | O(t · k) |
| Q3 | Hirschberg's Algorithm | O(m · n) | **O(m + n)** |

> Q2: t = number of sequences, n = sequence length, k = motif length  
> Q3: m, n = lengths of the two sequences being aligned

---

## Notes

- All FASTA files were sourced from [NCBI](https://www.ncbi.nlm.nih.gov/).
- Q2 uses only the first **500 bp** of each gene (5′ regulatory/coding region) to keep runtime tractable.
- Q3 uses only the first **1000 bp** of BRCA1 and BRCA2 coding regions for demonstration.
- Hirschberg's algorithm is preferred over standard Needleman-Wunsch for long sequences because its space usage is **linear** rather than quadratic.
