# Sequence Manipulation Tools

A collection of custom tools for DNA/RNA sequence manipulation.

## Table of Contents
1. [CpG Optimiser](#cpg-optimiser-cpg_optimiserpy)
2. [Some Assembly Required](#some-assembly-required-some_assembly_requiredr)

## Tools

### CpG Optimiser (`cpg_optimiser.py`)

Codon-optimises input ORFs to minimise the presence of CpG dinucleotides while preserving the original amino acid sequence. Designed for studying methylation dynamics of transgenic ORFs.

#### Features

- **CpG minimisation**: Removes canonical methylation targets (5'-CG-3')
- **Optional GpC reduction**: For non-canonical methylation studies (`--reduce-gpc` and `--reduce-gpc-only`, targets 5'-GC-3' dinucleotides).
- **Human codon optimisation**: Uses human codon frequency tables by default. Parameterisation of species is planned but currently not implemented.
- **Multiple algorithms**: Dynamic programming (optimal) or greedy (fast).
- **Batch processing**: Accepts FASTA files with multiple sequences.

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

# Reduce both CpG and GpC dinucleotides
python cpg_optimiser.py -s ATGCTGCCGTAA --reduce-gpc

# Reduce ONLY GpC dinucleotides (ignore CpG)
python cpg_optimiser.py -s ATGCTGCCGTAA --reduce-gpc-only

# Disable codon frequency weighting
python cpg_optimiser.py -s ATGCGTTCG --no-use-freq

# Quiet mode (output sequence only)
python cpg_optimiser.py -s ATGCGTTCG -q
```

#### Options

```bash
-s, --sequence      # Input DNA sequence
-f, --file          # Input .fasta file
-o, --output        # Output .fasta file (if not specified, prints to stdout)
--method            # Optimisation algorithm (dp [default] or greedy)
--reduce-gpc        # Target both CpG and GpC dinucleotides for optimisation
--reduce-gpc-only   # Only target GpC dinucleotides for optimisation
--no-use-freq       # Disable human codon frequency weighting
-q, --quiet         # Minimal printed output (no report)
--version           # Check version
```

#### Example output
```bash
======================================================================
  CpG/GpC OPTIMISATION REPORT: demo_ORF
======================================================================

  Mode: CpG + GpC reduction, with human codon frequencies

✓ Protein sequence preserved: YES ✓

METRIC                   ORIGINAL    OPTIMISED       CHANGE
--------------------------------------------------------
Length (bp)                    30           30            —
CpG count (target)              7            0           -7 (-100.0%)
GpC count (target)              3            1           -2 (-66.7%)
GC content (%)               60.0         56.7         -3.3
CpG/100bp                   23.33         0.00
GpC/100bp                   10.00         3.33

CODON SUBSTITUTIONS  8 of 10 codons changed

Position   Original   Optimised  Amino Acid
----------------------------------------
2          CGT        AGA        R
3          TCG        TCC        S
4          ACG        ACC        T
5          CCG        CCC        P
6          ACG        ACA        T
7          GCG        GCC        A
9          ACG        ACC        T
10         TAA        TGA        *

SEQUENCE ALIGNMENT
Original:  ATGCGTTCGACGCCGACGGCGATCACGTAA
Optimised: ATGAGATCCACCCCCACAGCCATCACCTGA
Changes:      ↑ ↑  ↑  ↑  ↑  ↑  ↑     ↑ ↑
CpG orig:     ^^  ^^ ^^ ^^ ^^ ^^    ^^
CpG opt:
GpC orig:    ^^       ^^     ^^
GpC opt:                     ^^

WARNINGS
  ⚠ Contains 1 AGA/AGG (Arg) codons - may affect translation in some systems
  ✓ All CpG sites eliminated
  ℹ 1 GpC sites remain (unavoidable for this amino acid sequence)
```

#### Background
As with most of the stuff I make, this tool's origin can be traced back to a brief chat over coffee. Due credit to the Australian Regenerative Medicine Institute, Monash University for the provision of a coffee machine, without which this tool wouldn't exist. BYO beans is still bullshit though.

This tool was developed to facilitate the study of methylation of transgenic ORFs. By creating CpG/GpC-depleted versions of coding sequences that encode identical proteins, we can disambiguate between methylation as a *cause* of ORF expression change, and methylation as a *consequence* of expression state.

#### Citation

If you use this tool, please cite:

Douek, A.M. (2026). seq_manipulation_tools: CpG/GpC-reducing ORF codon optimiser.
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

#### Changelog

**V1.1.0** Implemented `--reduce-gpc-only` flag to allow targeting of GpC while ignoring CpG (`--reduce-gpc` still targets both dinucleotides).

**V1.0.0** Initial release.

---

### Some Assembly Required (some_assembly_required.R)

End-to-end R script for forging a [BSgenome](https://bioconductor.org/packages/release/bioc/html/BSgenome.html) data package from any NCBI assembly accession (`GCF_*` or `GCA_*`). Designed to bypass some brittle steps in `BSgenomeForge::forgeBSgenomeDataPkgFromNCBI()` that fail on brand-new or unregistered assemblies.

#### Features
- **Robust against new assemblies**: Sidesteps the FASTA/assembly-report cross-check in `forgeBSgenomeDataPkgFromNCBI()` that fails when `GenomeInfoDb`'s cached chromosome info hasn't caught up with NCBI's live FASTA.
- **Smart chromosome naming**: Derives `chr*` names from `SequenceName` (always populated) rather than `AssignedMolecule` (often `NA` for new assemblies). Avoids a `chrNA` failure mode.
- **Automatic deduplication**: De-collides scaffold names that sanitisation might collapse (e.g. `Scaffold_1.5` and `Scaffold-1-5`) by appending `.1, .2, ...`.
- **Circular sequence detection**: Reads the `circular` column from chrom-info; works for mitochondrial, plastid, and bacterial assemblies. Falls back to name-based heuristics if needed.
- **Download integrity checks**: Tuneable download timeout (vs R's 60 s default) and MD5 verification against NCBI's published checksums. Catches silent truncation that would otherwise produce a partial BSgenome.
- **2bit pre-validation**: Smoke-tests the intermediate 2bit file before forging, so errors surface with useful messages for subsequent diagnostics as needed.
- **Generalised**: Configurable `chrom_prefix` (default `""`, set `"chr"` for assemblies whose `SequenceName` lacks the prefix) and `unplaced_prefix`. Works for any organism.
- **Verification helper**: Companion `verify_bsgenome()` function for sanity-checking the forged object against a user-supplied panel of (gene, expected chromosome, approximate coordinate) tuples.

#### Installation

```bash
# Clone the repository
git clone https://github.com/amdouek/seq_manipulation_tools.git
cd seq_manipulation_tools
```

Dependencies are installed automatically on first run via Bioconductor. Requires Rtools (Windows) or equivalent build tools for installing the forged package.

#### Usage

```r
# Edit the PARAMS block at the top of the script, then source:
source("some_assembly_required.R")

# Or call run_forge() directly with named arguments, for example:
run_forge(
    assembly_accession = "GCF_049306965.1",
    organism           = "Danio rerio",
    genome_name        = "GRCz12tu",
    pkg_maintainer     = "Your Name ",
    workdir            = "E:/custom_bsgenome/GRCz12tu"
)

# Verify the result:
panel <- data.frame(
    symbol           = c("sox2", "pax6a"),
    expected_chr     = c("chr22", "chr25"),
    approx_start_bp  = c(39480000, 8420000)
)
verify_bsgenome("BSgenome.Drerio.NCBI.GRCz12tu", panel)
```

#### Options

```r
assembly_accession   # NCBI accession, e.g. "GCF_049306965.1" or "GCA_052040795.1"
organism             # Full Latin binomial, e.g. "Danio rerio"
genome_name          # Short tag used in the package name, e.g. "GRCz12tu"
pkg_maintainer       # Maintainer string per R package conventions
workdir              # Output directory (created if absent)
install_after_forge  # If TRUE, build .tar.gz and install (default TRUE)
chrom_prefix         # Prefix for assembled molecules (default "chr"; "" if SequenceName already includes it)
unplaced_prefix      # Prefix for unplaced/unlocalised scaffolds (default "chrUn_")
keep_intermediate    # Keep cleaned FASTA + intermediate 2bit (default TRUE)
```

#### Background

NCBI's BSgenomeForge convenience function (`forgeBSgenomeDataPkgFromNCBI()`) sometimes fails ungracefully on brand-new assemblies whose live FASTA and cached chrom-info disagree on sequence counts. This script was built to account for this failure mode for newly-released NCBI-indexed genomes.

#### Citation

If you use this tool, please cite:

Douek, A.M. (2026). seq_manipulation_tools: NCBI-to-BSgenome forge script.
GitHub repository: https://github.com/amdouek/seq_manipulation_tools

#### Requirements

- R ≥ 4.3
- Bioconductor packages (auto-installed): `BSgenomeForge`, `Biostrings`, `rtracklayer`, `GenomeInfoDb`, `GenomicRanges`, `BSgenome`
- CRAN: `devtools` (for building the forged package)
- Build tools: Rtools (Windows), Xcode CLT (macOS), or `r-base-dev` (Linux)
- Memory: ~8 GB RAM for vertebrate-scale genomes
- Disk: 2–5× the compressed assembly size during forging

#### Licence

GPL-3.0

#### Author

Alon M Douek

#### Contributions

Contributions are always welcome! Please open an issue or submit a PR.

#### Changelog

**V1.0.0** Initial release. Handles new and unregistered NCBI assemblies, automatic deduplication of sanitised scaffold names, circular sequence detection, MD5-verified downloads with extended timeout, 2bit pre-validation, and a verification helper.