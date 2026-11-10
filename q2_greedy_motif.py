def read_fasta(filepath):
    """Reads a FASTA file and returns the header and uppercase DNA sequence."""
    header = ""
    seq_lines = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:]
            else:
                seq_lines.append(line.upper())
    return header, "".join(seq_lines)

def build_profile_matrix(motifs, k):
    """
    Constructs a 4 x k profile probability matrix from selected motifs
    using Laplace's Rule of Succession (+1 pseudocount) to avoid zero probabilities.
    """
    t = len(motifs)
    profile = {base: [0.0] * k for base in "ACGT"}
    for col in range(k):
        counts = {base: 1 for base in "ACGT"}  # Pseudocount = 1
        for row in range(t):
            counts[motifs[row][col]] += 1
        for base in "ACGT":
            profile[base][col] = counts[base] / (t + 4)
    return profile

def profile_most_probable_kmer(dna_seq, k, profile):
    """Finds the k-mer in dna_seq that has the highest probability given the profile matrix."""
    max_prob = -1.0
    best_kmer = dna_seq[0:k]
    for i in range(len(dna_seq) - k + 1):
        kmer = dna_seq[i:i + k]
        prob = 1.0
        for col, base in enumerate(kmer):
            prob *= profile[base][col]
        if prob > max_prob:
            max_prob = prob
            best_kmer = kmer
    return best_kmer

def score_motifs(motifs, k):
    """
    Calculates the consensus motif and the mismatch score (number of bases
    differing from the most frequent nucleotide in each column). Lower score = better.
    """
    t = len(motifs)
    consensus = []
    total_score = 0
    for col in range(k):
        counts = {base: 0 for base in "ACGT"}
        for row in range(t):
            counts[motifs[row][col]] += 1
        most_frequent_base = max(counts, key=counts.get)
        consensus.append(most_frequent_base)
        total_score += (t - counts[most_frequent_base])
    return "".join(consensus), total_score

def greedy_motif_search(dna_sequences, k):
    """
    Greedy Motif Search Algorithm across t DNA sequences for a motif of length k.
    """
    t = len(dna_sequences)
    # Initialize best motifs with the first k-mer of each sequence
    best_motifs = [seq[0:k] for seq in dna_sequences]
    _, best_score = score_motifs(best_motifs, k)

    # Iterate through every k-mer in the first DNA sequence
    first_seq = dna_sequences[0]
    for i in range(len(first_seq) - k + 1):
        current_motifs = [first_seq[i:i + k]]
        # Greedily pick the profile-most probable k-mer from sequences 2 to t
        for j in range(1, t):
            profile = build_profile_matrix(current_motifs, k)
            next_motif = profile_most_probable_kmer(dna_sequences[j], k, profile)
            current_motifs.append(next_motif)

        _, current_score = score_motifs(current_motifs, k)
        if current_score < best_score:
            best_motifs = current_motifs
            best_score = current_score

    consensus, final_score = score_motifs(best_motifs, k)
    return best_motifs, consensus, final_score

if __name__ == "__main__":
    fasta_files = ["HBB.fasta", "CFTR.fasta", "HTT.fasta", "TP53.fasta"]
    headers = []
    dna_list = []

    for file in fasta_files:
        header, seq = read_fasta(file)
        headers.append(header.split()[0])
        # Use first 500 bp (5' regulatory/coding region) of each gene
        dna_list.append(seq[:500])

    k = 10  # Motif length (10 bp)
    best_motifs, consensus, final_score = greedy_motif_search(dna_list, k)

    total_positions = k * len(dna_list)
    conservation_pct = ((total_positions - final_score) / total_positions) * 100

    print("==========================================================")
    print("       QUESTION 2: GREEDY MOTIF SEARCH RESULTS            ")
    print("==========================================================")
    for gene_file, acc, motif in zip(fasta_files, headers, best_motifs):
        print(f"{gene_file:12s} ({acc:12s}) -> Motif: {motif}")
    print("----------------------------------------------------------")
    print(f"Selected Consensus Motif : {consensus}")
    print(f"Motif Mismatch Score     : {final_score} mismatches (out of {total_positions} bases)")
    print(f"Positional Conservation  : {conservation_pct:.2f}%")
    print("==========================================================")