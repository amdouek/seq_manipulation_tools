"""
CpG-Reduced Codon Optimiser
===========================
Optimises coding sequences to minimise CpG dinucleotides while preserving 
the original amino acid sequence. Designed for studying methylation dynamics of 
transgenic ORFs.

Author: Alon M Douek
Usage:
    python cpg_optimiser.py -s ATGCGTTCGACG...
    python cpg_optimiser.py -f input.fasta -o output.fasta
    python cpg_optimiser.py -s SEQUENCE --method dp --verbose
    python cpg_optimiser.py -s SEQUENCE --reduce-gpc
"""

import argparse
import sys
from typing import Dict, List, Tuple, Optional
from collections import defaultdict

# =============================================================================
# CODON TABLES
# =============================================================================

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F',
    'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
    'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y',
    'TAA': '*', 'TAG': '*', 'TGA': '*',
    'CAT': 'H', 'CAC': 'H',
    'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D',
    'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C',
    'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# Build reverse lookup: amino acid -> list of codons
AA_TO_CODONS: Dict[str, List[str]] = defaultdict(list)
for codon, aa in CODON_TABLE.items():
    AA_TO_CODONS[aa].append(codon)

# Human codon usage frequencies (per thousand) - optional weighting
# Source: Kazusa codon usage database
HUMAN_CODON_FREQ = {
    'TTT': 17.6, 'TTC': 20.3, 'TTA': 7.7, 'TTG': 12.9,
    'CTT': 13.2, 'CTC': 19.6, 'CTA': 7.2, 'CTG': 39.6,
    'ATT': 16.0, 'ATC': 20.8, 'ATA': 7.5, 'ATG': 22.0,
    'GTT': 11.0, 'GTC': 14.5, 'GTA': 7.1, 'GTG': 28.1,
    'TCT': 15.2, 'TCC': 17.7, 'TCA': 12.2, 'TCG': 4.4,
    'CCT': 17.5, 'CCC': 19.8, 'CCA': 16.9, 'CCG': 6.9,
    'ACT': 13.1, 'ACC': 18.9, 'ACA': 15.1, 'ACG': 6.1,
    'GCT': 18.4, 'GCC': 27.7, 'GCA': 15.8, 'GCG': 7.4,
    'TAT': 12.2, 'TAC': 15.3, 'TAA': 1.0, 'TAG': 0.8,
    'CAT': 10.9, 'CAC': 15.1, 'CAA': 12.3, 'CAG': 34.2,
    'AAT': 17.0, 'AAC': 19.1, 'AAA': 24.4, 'AAG': 31.9,
    'GAT': 21.8, 'GAC': 25.1, 'GAA': 29.0, 'GAG': 39.6,
    'TGT': 10.6, 'TGC': 12.6, 'TGA': 1.6, 'TGG': 13.2,
    'CGT': 4.5, 'CGC': 10.4, 'CGA': 6.2, 'CGG': 11.4,
    'AGT': 12.1, 'AGC': 19.5, 'AGA': 12.2, 'AGG': 12.0,
    'GGT': 10.8, 'GGC': 22.2, 'GGA': 16.5, 'GGG': 16.5,
}

# =============================================================================
# CORE FUNCTIONS
# =============================================================================

def clean_sequence(sequence: str) -> str:
    """Remove whitespace and convert to uppercase."""
    return ''.join(sequence.upper().split())


def translate(sequence: str) -> str:
    """Translate DNA sequence to amino acid sequence."""
    sequence = clean_sequence(sequence)
    protein = []
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        aa = CODON_TABLE.get(codon, 'X')
        protein.append(aa)
    return ''.join(protein)


def count_cpg(sequence: str) -> int:
    """Count CpG dinucleotides in a sequence."""
    return sequence.upper().count('CG')


def count_gpc(sequence: str) -> int:
    """Count GpC dinucleotides in a sequence."""
    return sequence.upper().count('GC')


def get_codon_cpg_score(codon: str, prev_codon: Optional[str] = None, 
                        reduce_gpc: bool = False) -> int:
    """
    Calculate methylation target score for a codon.
    
    Args:
        codon: Current codon (3 nucleotides)
        prev_codon: Previous codon (for boundary checking)
        reduce_gpc: If True, also penalise GpC dinucleotides
    
    Returns:
        Score (lower is better)
    
    Scoring:
        - CpG (5'-CG-3'): Canonical mammalian methylation target (always counted)
        - GpC (5'-GC-3'): Non-canonical, optional (only if reduce_gpc=True)
    """
    # Always count CpG (canonical methylation target)
    score = codon.count('CG')
    
    # Check for cross-boundary CpG: previous ends with C, current starts with G
    if prev_codon and prev_codon[-1] == 'C' and codon[0] == 'G':
        score += 1
    
    # Optionally count GpC (non-canonical methylation)
    if reduce_gpc:
        score += codon.count('GC')
        # Cross-boundary GpC: previous ends with G, current starts with C
        if prev_codon and prev_codon[-1] == 'G' and codon[0] == 'C':
            score += 1
    
    return score


# =============================================================================
# OPTIMISATION ALGORITHMS
# =============================================================================

def optimise_greedy(sequence: str, use_freq: bool = False, 
                    reduce_gpc: bool = False) -> str:
    """
    Greedy codon optimisation to minimise CpG (and optionally GpC).
    
    For each position, selects the codon with minimum methylation target score,
    considering the previous codon for boundary effects.
    
    Args:
        sequence: Input DNA sequence (must be divisible by 3)
        use_freq: Weight by human codon usage frequency
        reduce_gpc: Also minimise GpC dinucleotides
    
    Returns:
        Optimised DNA sequence
    """
    sequence = clean_sequence(sequence)
    
    if len(sequence) % 3 != 0:
        raise ValueError(f"Sequence length ({len(sequence)}) must be divisible by 3")
    
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
    amino_acids = [CODON_TABLE.get(c, 'X') for c in codons]
    
    optimised = []
    
    for i, aa in enumerate(amino_acids):
        candidates = AA_TO_CODONS.get(aa, [codons[i]])
        prev_codon = optimised[-1] if optimised else None
        
        # Score each candidate
        scored = []
        for c in candidates:
            # Primary score: CpG (always minimised)
            cpg_score = c.count('CG')
            if prev_codon and prev_codon[-1] == 'C' and c[0] == 'G':
                cpg_score += 1
            
            # Secondary score: GpC (only used for tiebreaking if reduce_gpc=True)
            if reduce_gpc:
                gpc_score = c.count('GC')
                if prev_codon and prev_codon[-1] == 'G' and c[0] == 'C':
                    gpc_score += 1
            else:
                gpc_score = 0  # Don't consider GpC in tiebreaking
            
            # Tertiary score: codon frequency (negative because higher is better)
            freq_score = -HUMAN_CODON_FREQ.get(c, 10) if use_freq else 0
            
            # Tuple for sorting: (CpG, GpC, -frequency, codon)
            scored.append((cpg_score, gpc_score, freq_score, c))
        
        # Sort by CpG (primary), then GpC (secondary), then frequency (tertiary)
        scored.sort()
        optimised.append(scored[0][3])
    
    return ''.join(optimised)


def optimise_dp(sequence: str, use_freq: bool = False, 
                reduce_gpc: bool = False) -> str:
    """
    Dynamic programming optimisation for globally minimal CpG count.
    
    Considers all possible codon combinations and finds the sequence
    with the absolute minimum number of CpG dinucleotides (including
    cross-codon boundaries). When reduce_gpc=True, GpC is used as a
    secondary optimisation target (tiebreaker).
    
    Args:
        sequence: Input DNA sequence (must be divisible by 3)
        use_freq: Add penalty for rare codons
        reduce_gpc: Also minimise GpC dinucleotides (secondary to CpG)
    
    Returns:
        Optimised DNA sequence
    """
    sequence = clean_sequence(sequence)
    
    if len(sequence) % 3 != 0:
        raise ValueError(f"Sequence length ({len(sequence)}) must be divisible by 3")
    
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
    amino_acids = [CODON_TABLE.get(c, 'X') for c in codons]
    n = len(amino_acids)
    
    if n == 0:
        return ""
    
    # Get codon options for each position
    codon_options = [AA_TO_CODONS.get(aa, ['NNN']) for aa in amino_acids]
    
    # Epsilon weight for GpC - small enough to never override CpG differences
    # Max possible GpC in a sequence is roughly n*2, so 0.001 ensures
    # that even n*2*0.001 < 1 (one CpG) for sequences up to 500 codons
    GPC_EPSILON = 0.001 if reduce_gpc else 0
    
    # Frequency weight - even smaller, only for tiebreaking after CpG and GpC
    FREQ_EPSILON = 0.0001 if use_freq else 0
    
    # DP tables
    INF = float('inf')
    dp = [[INF] * len(codon_options[i]) for i in range(n)]
    parent = [[None] * len(codon_options[i]) for i in range(n)]
    
    # Base case: first codon (only internal CpG/GpC)
    for j, codon in enumerate(codon_options[0]):
        # Primary: CpG count (integer weight = 1)
        score = codon.count('CG')
        
        # Secondary: GpC count (small epsilon weight)
        score += codon.count('GC') * GPC_EPSILON
        
        # Tertiary: frequency penalty (even smaller epsilon)
        if use_freq:
            score += (1 - HUMAN_CODON_FREQ.get(codon, 10) / 40) * FREQ_EPSILON
        
        dp[0][j] = score
    
    # Fill DP table
    for i in range(1, n):
        for j, curr_codon in enumerate(codon_options[i]):
            # Internal CpG (weight = 1)
            internal_cpg = curr_codon.count('CG')
            
            # Internal GpC (weight = epsilon)
            internal_gpc = curr_codon.count('GC') * GPC_EPSILON
            
            for k, prev_codon in enumerate(codon_options[i-1]):
                # Cross-boundary CpG (weight = 1)
                boundary_cpg = 1 if (prev_codon[-1] == 'C' and curr_codon[0] == 'G') else 0
                
                # Cross-boundary GpC (weight = epsilon)
                boundary_gpc = GPC_EPSILON if (prev_codon[-1] == 'G' and curr_codon[0] == 'C') else 0
                
                # Total score
                score = dp[i-1][k] + internal_cpg + internal_gpc + boundary_cpg + boundary_gpc
                
                # Frequency penalty
                if use_freq:
                    score += (1 - HUMAN_CODON_FREQ.get(curr_codon, 10) / 40) * FREQ_EPSILON
                
                if score < dp[i][j]:
                    dp[i][j] = score
                    parent[i][j] = k
    
    # Backtrack
    best_end = min(range(len(codon_options[n-1])), key=lambda j: dp[n-1][j])
    
    result = []
    curr = best_end
    for i in range(n-1, -1, -1):
        result.append(codon_options[i][curr])
        if i > 0:
            curr = parent[i][curr]
    
    result.reverse()
    return ''.join(result)


# =============================================================================
# ANALYSIS AND REPORTING
# =============================================================================

def analyse_sequence(sequence: str) -> dict:
    """Calculate various sequence metrics."""
    seq = clean_sequence(sequence)
    length = len(seq)
    
    if length == 0:
        return {'length': 0, 'cpg_count': 0, 'gpc_count': 0, 
                'gc_content': 0, 'cpg_density': 0}
    
    gc_count = seq.count('G') + seq.count('C')
    
    return {
        'length': length,
        'cpg_count': count_cpg(seq),
        'gpc_count': count_gpc(seq),
        'gc_content': gc_count / length * 100,
        'cpg_density': count_cpg(seq) / (length / 100),
    }


def highlight_differences(seq1: str, seq2: str) -> str:
    """Create a string highlighting differences between two sequences."""
    return ''.join([' ' if a == b else '↑' for a, b in zip(seq1, seq2)])


def highlight_cpg(sequence: str) -> str:
    """Create a string marking CpG positions."""
    result = []
    seq = sequence.upper()
    i = 0
    while i < len(seq):
        if i < len(seq) - 1 and seq[i:i+2] == 'CG':
            result.append('^^')
            i += 2
        else:
            result.append(' ')
            i += 1
    return ''.join(result)


def format_codons(sequence: str, per_line: int = 20) -> str:
    """Format sequence with codon spacing."""
    seq = clean_sequence(sequence)
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    
    lines = []
    for i in range(0, len(codons), per_line):
        lines.append(' '.join(codons[i:i+per_line]))
    
    return '\n'.join(lines)


def print_report(original: str, optimised: str, name: str = "Sequence", 
                 reduce_gpc: bool = False, use_freq: bool = True):
    """Print a comprehensive optimisation report."""
    orig_analysis = analyse_sequence(original)
    opt_analysis = analyse_sequence(optimised)
    
    orig_protein = translate(original)
    opt_protein = translate(optimised)
    
    print(f"\n{'='*70}")
    print(f"  CpG OPTIMISATION REPORT: {name}")
    print(f"{'='*70}")
    
    # Mode indicator
    gpc_status = "CpG + GpC" if reduce_gpc else "CpG only"
    freq_status = "human codon frequencies" if use_freq else "no frequency weighting"
    print(f"\n  Mode: {gpc_status} reduction, {freq_status}")
    
    # Verification
    protein_match = orig_protein == opt_protein
    print(f"\n✓ Protein sequence preserved: {'YES ✓' if protein_match else 'NO ✗ ERROR!'}")
    
    if not protein_match:
        print(f"  Original protein:  {orig_protein}")
        print(f"  Optimised protein: {opt_protein}")
        return
    
    # Metrics comparison
    print(f"\n{'METRIC':<20} {'ORIGINAL':>12} {'OPTIMISED':>12} {'CHANGE':>12}")
    print(f"{'-'*56}")
    print(f"{'Length (bp)':<20} {orig_analysis['length']:>12} {opt_analysis['length']:>12} {'—':>12}")
    
    # CpG (always targeted)
    cpg_change = opt_analysis['cpg_count'] - orig_analysis['cpg_count']
    cpg_pct = (cpg_change / orig_analysis['cpg_count'] * 100) if orig_analysis['cpg_count'] > 0 else 0
    print(f"{'CpG count (target)':<20} {orig_analysis['cpg_count']:>12} {opt_analysis['cpg_count']:>12} {cpg_change:>+12} ({cpg_pct:+.1f}%)")
    
    # GpC (targeted only if reduce_gpc=True)
    gpc_change = opt_analysis['gpc_count'] - orig_analysis['gpc_count']
    if reduce_gpc:
        gpc_pct = (gpc_change / orig_analysis['gpc_count'] * 100) if orig_analysis['gpc_count'] > 0 else 0
        print(f"{'GpC count (target)':<20} {orig_analysis['gpc_count']:>12} {opt_analysis['gpc_count']:>12} {gpc_change:>+12} ({gpc_pct:+.1f}%)")
    else:
        print(f"{'GpC count (info)':<20} {orig_analysis['gpc_count']:>12} {opt_analysis['gpc_count']:>12} {gpc_change:>+12} (not targeted)")
    
    gc_change = opt_analysis['gc_content'] - orig_analysis['gc_content']
    print(f"{'GC content (%)':<20} {orig_analysis['gc_content']:>12.1f} {opt_analysis['gc_content']:>12.1f} {gc_change:>+12.1f}")
    
    print(f"{'CpG/100bp':<20} {orig_analysis['cpg_density']:>12.2f} {opt_analysis['cpg_density']:>12.2f}")
    
    # Codon changes
    orig_codons = [original[i:i+3] for i in range(0, len(original), 3)]
    opt_codons = [optimised[i:i+3] for i in range(0, len(optimised), 3)]
    
    changes = [(i, orig_codons[i], opt_codons[i]) 
               for i in range(len(orig_codons)) 
               if orig_codons[i] != opt_codons[i]]
    
    print(f"\n{'CODON SUBSTITUTIONS':<20} {len(changes)} of {len(orig_codons)} codons changed")
    
    if changes and len(changes) <= 30:
        print(f"\n{'Position':<10} {'Original':<10} {'Optimised':<10} {'Amino Acid':<10}")
        print(f"{'-'*40}")
        for pos, orig, opt in changes:
            aa = CODON_TABLE.get(orig, '?')
            print(f"{pos+1:<10} {orig:<10} {opt:<10} {aa:<10}")
    
    # Sequence alignment (for short sequences)
    if len(original) <= 120:
        print(f"\n{'SEQUENCE ALIGNMENT':}")
        print(f"Original:  {original}")
        print(f"Optimised: {optimised}")
        print(f"Changes:   {highlight_differences(original, optimised)}")
        print(f"CpG orig:  {highlight_cpg(original)}")
        print(f"CpG opt:   {highlight_cpg(optimised)}")
    
    # Warnings
    print(f"\n{'WARNINGS':}")
    
    # Check for rare Arg codons
    aga_count = sum(1 for i in range(0, len(optimised), 3) if optimised[i:i+3] in ('AGA', 'AGG'))
    if aga_count > 0:
        print(f"  ⚠ Contains {aga_count} AGA/AGG (Arg) codons - may affect translation in some systems")
    
    if opt_analysis['gc_content'] < 40:
        print(f"  ⚠ Low GC content ({opt_analysis['gc_content']:.1f}%) - may affect mRNA stability")
    
    if opt_analysis['cpg_count'] > 0:
        print(f"  ℹ {opt_analysis['cpg_count']} CpG sites remain (unavoidable for this amino acid sequence)")
    else:
        print(f"  ✓ All CpG sites eliminated")


# =============================================================================
# FILE I/O
# =============================================================================

def parse_fasta(filename: str) -> List[Tuple[str, str]]:
    """Parse a FASTA file."""
    sequences = []
    current_name = None
    current_seq = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name:
                    sequences.append((current_name, ''.join(current_seq)))
                current_name = line[1:].split()[0]
                current_seq = []
            elif line and not line.startswith('#'):
                current_seq.append(line)
        
        if current_name:
            sequences.append((current_name, ''.join(current_seq)))
    
    return sequences


def write_fasta(filename: str, sequences: List[Tuple[str, str]], line_width: int = 60):
    """Write sequences to a FASTA file."""
    with open(filename, 'w') as f:
        for name, seq in sequences:
            f.write(f">{name}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i+line_width] + '\n')


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='CpG-reduced codon optimiser for transgenic ORF methylation studies',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -s ATGCGTTCGACGCCGACGGCGATCACGTAA
  %(prog)s -f input.fasta -o output.fasta
  %(prog)s -s ATGCGTACG --method greedy 
  %(prog)s -s ATGCGTACG --reduce-gpc
  %(prog)s -s ATGCGTACG --no-use-freq
  %(prog)s -s ATGCGTACG --quiet
  
Notes:
  - Input sequences must be in-frame (length divisible by 3)
  - The 'dp' method (default) finds the globally optimal solution
  - The 'greedy' method is faster but may miss some optimisations
  - Human codon usage frequences are applied by default for tiebreaking (disable with --no-use-freq)
  - Use --reduce-gpc to also minimise GpC (non-canonical methylation target)
        """
    )
    
    input_group = parser.add_mutually_exclusive_group()
    input_group.add_argument('-s', '--sequence', type=str,
                             help='Input DNA sequence (coding strand, 5\' to 3\')')
    input_group.add_argument('-f', '--file', type=str,
                             help='Input FASTA file')
    
    parser.add_argument('-o', '--output', type=str,
                        help='Output FASTA file (if not specified, prints to stdout)')
    parser.add_argument('--method', choices=['greedy', 'dp'], default='dp',
                        help='Optimisation algorithm (default: dp)')
    
    # Frequency weighting: ON by default, can be disabled
    freq_group = parser.add_mutually_exclusive_group()
    freq_group.add_argument('--use-freq', action='store_true', dest='use_freq',
                            help='Use human codon frequency weighting (default: enabled)')
    freq_group.add_argument('--no-use-freq', action='store_false', dest='use_freq',
                            help='Disable codon frequency weighting')
    parser.set_defaults(use_freq=True)  # ← Default is now ON
    
    parser.add_argument('--reduce-gpc', action='store_true',
                        help='Also minimise GpC dinucleotides (non-canonical methylation)')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='Minimal output (just the optimised sequence)')
    
    args = parser.parse_args()
    
    # Get input sequences
    if args.sequence:
        sequences = [('input', args.sequence)]
    elif args.file:
        try:
            sequences = parse_fasta(args.file)
        except FileNotFoundError:
            print(f"Error: File '{args.file}' not found", file=sys.stderr)
            sys.exit(1)
    else:
        # Demo mode
        print("No input provided. Running demonstration...\n")
        demo_seq = "ATGCGTTCGACGCCGACGGCGATCACGTAA"
        sequences = [('demo_ORF', demo_seq)]
        args.quiet = False
    
    # Process sequences
    results = []
    
    for name, seq in sequences:
        seq = clean_sequence(seq)
        
        # Validate
        if not seq:
            print(f"Warning: Empty sequence '{name}', skipping", file=sys.stderr)
            continue
            
        if len(seq) % 3 != 0:
            print(f"Warning: '{name}' length ({len(seq)}) not divisible by 3, truncating", 
                  file=sys.stderr)
            seq = seq[:len(seq) - (len(seq) % 3)]
        
        # Validate nucleotides
        invalid = set(seq) - set('ATCG')
        if invalid:
            print(f"Warning: '{name}' contains invalid characters: {invalid}", file=sys.stderr)
            continue
        
        # Optimise
        if args.method == 'dp':
            optimised = optimise_dp(seq, args.use_freq, args.reduce_gpc)
        else:
            optimised = optimise_greedy(seq, args.use_freq, args.reduce_gpc)
        
        # Verify
        if translate(seq) != translate(optimised):
            print(f"ERROR: Protein sequence changed for '{name}'!", file=sys.stderr)
            continue
        
        # Generate output name
        suffix = "_CpG_reduced" if not args.reduce_gpc else "_CpG_GpC_reduced"
        results.append((f"{name}{suffix}", optimised))
        
        # Report
        if not args.quiet:
            print_report(seq, optimised, name, args.reduce_gpc, args.use_freq)
    
    # Output
    if args.output:
        write_fasta(args.output, results)
        print(f"\nOptimised sequences written to: {args.output}")
    elif args.quiet:
        for name, seq in results:
            print(f">{name}")
            print(seq)


if __name__ == '__main__':
    main()