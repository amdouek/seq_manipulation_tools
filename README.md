# Sequence Manipulation Tools

A collection of Python tools for DNA sequence manipulation, with a focus on 
codon optimisation for transgenic expression studies.

## Tools

### CpG Optimiser (`cpg_optimiser.py`)

Codon-optimises input ORFs to minimise the presence of CpG dinucleotides while preserving the original amino acid sequence. Designed for studying methylation dynamics of transgenic ORFs.

#### Features

- **CpG minimisation**: Removes canonical methylation targets (5'-CG-3')
- **Optional GpC reduction**: For non-canonical methylation studies (`--reduce-gpc`, targets 5'-GC-3' dinucleotides)
- **Human codon optimisation**: Uses human codon frequency tables by default. Parameterisation of species is planned but currently not implemented.
- **Multiple algorithms**: Dynamic programming (optimal) or greedy (fast)
- **Batch processing**: Accepts FASTA files with multiple sequences

#### Installation

```bash
# Clone the repository
git clone https://github.com/amdouek/seq_manipulation_tools.git
cd seq_manipulation_tools

# No additional dependencies required (uses Python standard library only)
```

#### Usage

```bash
# Basic usage with a sequence
python cpg_optimiser.py -s ATGCGTTCGACGCCGACGGCGATCACGTAA

# Process a FASTA file
python cpg_optimiser.py -f input.fasta -o output.fasta

# Include GpC reduction (for neuronal/stem cell contexts)
python cpg_optimiser.py -s ATGCTGCCGTAA --reduce-gpc

# Disable codon frequency weighting
python cpg_optimiser.py -s ATGCGTTCG --no-use-freq

# Quiet mode (output sequence only)
python cpg_optimiser.py -s ATGCGTTCG -q
```

#### Options

```bash
-s, --sequence  # Input DNA sequence
-f, --file      # Input .fasta file
-o, --output    # Output .fasta file (if not specified, prints to stdout)
--method        # Optimisation algorithm (dp [default] or greedy)
--reduce-gpc    # Also target GpC dinucleotides for optimisation
--no-use-freq   # Disable human codon frequency weighting
-q, --quiet     # Minimal printed output (no report)
```

#### Example output
```bash
======================================================================
  CpG OPTIMISATION REPORT: demo_ORF
======================================================================

  Mode: CpG only reduction, human codon frequencies

✓ Protein sequence preserved: YES ✓

METRIC                   ORIGINAL    OPTIMISED       CHANGE
--------------------------------------------------------
Length (bp)                    30           30            —
CpG count (target)              7            0     -7 (-100.0%)
GpC count (info)                4            2     -2 (not targeted)
GC content (%)               56.7         40.0        -16.7
```

#### Background

This tool was developed to facilitate the study of methylation of transgenic ORFs. By creating CpG-depleted versions of coding sequences that encode identical proteins, we can disambiguate between methylation as a *cause* of ORF expression change, and methylation as a *consequence* of expression state.

#### Citation

If you use this tool, please cite:

Douek, A.M. (2024). seq-manipulation-tools: CpG-reducing ORF codon optimiser.
GitHub repository: https://github.com/amdouek/seq_manipulation_tools


#### Requirements

- Python 3.8+
- No external dependencies

#### Licence

GPL-3.0

#### Author

Alon M Douek

#### Contributions

Contributions are always welcome! Please open an issue or submit a PR.