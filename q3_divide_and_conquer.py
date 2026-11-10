MATCH = 2
MISMATCH = -1
GAP = -2

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

def linear_space_score(X, Y):
    """Computes the final row of DP alignment scores using only O(len(Y)) memory."""
    n = len(Y)
    prev_row = [j * GAP for j in range(n + 1)]
    for x_char in X:
        curr_row = [prev_row[0] + GAP]
        for j in range(1, n + 1):
            diag = prev_row[j - 1] + (MATCH if x_char == Y[j - 1] else MISMATCH)
            up = prev_row[j] + GAP
            left = curr_row[j - 1] + GAP
            curr_row.append(max(diag, up, left))
        prev_row = curr_row
    return prev_row

def base_case_align(X, Y):
    """Standard Needleman-Wunsch alignment for base cases where len(X) <= 1 or len(Y) <= 1."""
    m, n = len(X), len(Y)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(m + 1):
        dp[i][0] = i * GAP
    for j in range(n + 1):
        dp[0][j] = j * GAP

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diag = dp[i - 1][j - 1] + (MATCH if X[i - 1] == Y[j - 1] else MISMATCH)
            up = dp[i - 1][j] + GAP
            left = dp[i][j - 1] + GAP
            dp[i][j] = max(diag, up, left)

    align_X, align_Y = [], []
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0 and dp[i][j] == dp[i - 1][j - 1] + (MATCH if X[i - 1] == Y[j - 1] else MISMATCH):
            align_X.append(X[i - 1])
            align_Y.append(Y[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and dp[i][j] == dp[i - 1][j] + GAP:
            align_X.append(X[i - 1])
            align_Y.append("-")
            i -= 1
        else:
            align_X.append("-")
            align_Y.append(Y[j - 1])
            j -= 1
    return "".join(reversed(align_X)), "".join(reversed(align_Y))

def divide_and_conquer_align(X, Y):
    """Hirschberg's Divide and Conquer Global Sequence Alignment Algorithm."""
    if len(X) == 0:
        return "-" * len(Y), Y
    elif len(Y) == 0:
        return X, "-" * len(X)
    elif len(X) == 1 or len(Y) == 1:
        return base_case_align(X, Y)
    else:
        mid = len(X) // 2

        # Conquer: Forward DP on first half, Backward DP on reversed second half
        score_left = linear_space_score(X[:mid], Y)
        score_right = linear_space_score(X[mid:][::-1], Y[::-1])
        score_right.reverse()

        # Combine: Find index k_star maximizing forward + backward score
        combined = [l + r for l, r in zip(score_left, score_right)]
        k_star = combined.index(max(combined))

        # Recursively align subproblems
        left_X, left_Y = divide_and_conquer_align(X[:mid], Y[:k_star])
        right_X, right_Y = divide_and_conquer_align(X[mid:], Y[k_star:])

        return left_X + right_X, left_Y + right_Y

if __name__ == "__main__":
    h1, brca1_full = read_fasta("BRCA1.fasta")
    h2, brca2_full = read_fasta("BRCA2.fasta")

    # Aligning the first 1000 bp coding region of BRCA1 and BRCA2
    seq1 = brca1_full[:1000]
    seq2 = brca2_full[:1000]

    aligned_brca1, aligned_brca2 = divide_and_conquer_align(seq1, seq2)

    similarity_pairs = sum(1 for a, b in zip(aligned_brca1, aligned_brca2) if a == b and a != "-")
    mismatches = sum(1 for a, b in zip(aligned_brca1, aligned_brca2) if a != b and a != "-" and b != "-")
    gaps = sum(1 for a, b in zip(aligned_brca1, aligned_brca2) if a == "-" or b == "-")
    total_len = len(aligned_brca1)
    identity_pct = (similarity_pairs / total_len) * 100

    print("==========================================================")
    print("  QUESTION 3: DIVIDE & CONQUER ALIGNMENT (BRCA1 vs BRCA2) ")
    print("==========================================================")
    print(f"BRCA1 Sequence          : {h1.split()[0]} (Length analyzed: {len(seq1)} bp)")
    print(f"BRCA2 Sequence          : {h2.split()[0]} (Length analyzed: {len(seq2)} bp)")
    print(f"Total Alignment Length  : {total_len} columns")
    print(f"Total Similarity Pairs  : {similarity_pairs} exact nucleotide matches ({identity_pct:.2f}%)")
    print(f"Total Mismatches        : {mismatches}")
    print(f"Total Gaps (Indels)     : {gaps}")
    print("----------------------------------------------------------")
    print("Alignment Preview (First 75 positions):")
    print("BRCA1 :", aligned_brca1[:75])
    print("Match :", "".join("|" if a == b and a != "-" else " " for a, b in zip(aligned_brca1[:75], aligned_brca2[:75])))
    print("BRCA2 :", aligned_brca2[:75])
    print("==========================================================")