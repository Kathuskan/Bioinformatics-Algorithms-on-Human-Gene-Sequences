
def read_fasta(filepath):
    """Reads a FASTA file and extracts the header and uppercase sequence."""
    header = ""
    seq_lines = []
    with open(filepath, "r") as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line
            else:
                seq_lines.append(line.upper())
    return header, "".join(seq_lines)

def calculate_gc_content(seq):
    """Calculates G, C counts and GC percentage."""
    g_count = seq.count("G")
    c_count = seq.count("C")
    gc_percent = ((g_count + c_count) / len(seq)) * 100
    return g_count, c_count, gc_percent

def get_reverse_complement(seq):
    """Generates the 5' -> 3' reverse complement of a DNA sequence."""
    complement_map = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    complement = "".join(complement_map.get(base, base) for base in seq)
    return complement[::-1]  # Reverse the complement string

if __name__ == "__main__":
    # Make sure 'sequence.fasta' is in the same folder as this script
    fasta_file = "INS_sequence.fasta"
    header, dna_seq = read_fasta(fasta_file)

    total_len = len(dna_seq)
    a_count = dna_seq.count("A")
    t_count = dna_seq.count("T")
    g_count, c_count, gc_pct = calculate_gc_content(dna_seq)
    rev_comp_seq = get_reverse_complement(dna_seq)

    print("==========================================================")
    print("       QUESTION 1: HUMAN INSULIN (INS) GENE ANALYSIS      ")
    print("==========================================================")
    print(f"FASTA Header      : {header}")
    print(f"Accession Number  : {header.split()[0][1:]}")
    print(f"Total Length      : {total_len} bp")
    print(f"Base Counts       : A={a_count}, T={t_count}, G={g_count}, C={c_count}")
    print(f"Total G + C       : {g_count + c_count} bp")
    print(f"1(d) GC Content   : {gc_pct:.2f}%")
    print("----------------------------------------------------------")
    print("1(e) Reverse Complement (5' -> 3', 70 bp per line):")
    for i in range(0, len(rev_comp_seq), 70):
        print(rev_comp_seq[i:i+70])
    print("==========================================================")