###############################################################################
##  **some_assembly_required.R**
##  Forge a BSgenome data package from any NCBI assembly
##  =====================================================
##
##  An end-to-end workflow that takes an NCBI assembly accession
##  (GCF_*** or GCA_***) and produces an installed BSgenome data package
##  with clean, conventional chromosome names.
##
##  This script bakes in workarounds for several issues encountered when
##  forging brand-new NCBI assemblies (where GenomeInfoDb's cached
##  chromosome-info can drift from what NCBI serves, and where the
##  AssignedMolecule column may be NA for newly-deposited assemblies):
##
##    - forgeBSgenomeDataPkgFromNCBI() is bypassed; we go FASTA -> 2bit ->
##      forgeBSgenomeDataPkgFromTwobitFile() to avoid its assembly-report
##      cross-check.
##    - Chromosome names are derived from the SequenceName column (which
##      is always populated) rather than AssignedMolecule (which may be
##      NA for the nuclear chromosomes of new assemblies).
##    - Sequence renaming is followed by a duplicate check; collisions
##      from sanitised scaffold names are de-duplicated by suffix.
##    - Circular sequences are read from the chrom-info `circular` column,
##      not hardcoded to "MT" (so this works for bacteria, plastids,
##      assemblies that use "MtDNA"/"mitochondrion", etc).
##    - The 2bit is validated before forging.
##    - Stale outputs are removed before each regeneration step.
##
##  Usage
##  -----
##  Edit the USER-CONFIGURABLE section below and source the script. Or
##  call run_forge() with named arguments.
##
##  Tested on Windows with Rtools. Memory: budget ~8 Gb RAM for vertebrate
##  genomes (the cleaned FASTA is loaded into a DNAStringSet at one point).
##
##  Author: Alon Douek <Alon.Douek@monash.edu>
###############################################################################

## ===== 0. USER-CONFIGURABLE PARAMETERS =======================================
## Edit these for each new assembly, or pass them to run_forge() directly.

PARAMS <- list(
  assembly_accession = "",        # "GCF_*" or "GCA_*"
  organism           = "",        # full Latin binomial, e.g. "Danio rerio"
  genome_name        = "",        # short tag for package name, e.g. "GRCz12tu"
  pkg_maintainer     = "Your Name <Your.Name@affiliation.org>",
  workdir            = "",        # Path to save to, e.g. "E:/custom_bsgenome/GRCz12tu"

  ## Optional advanced settings -- safe defaults
  install_after_forge = TRUE,     # build .tar.gz and install via R CMD?
  chrom_prefix        = "",       # default "", set if you want an additional prefix
  unplaced_prefix     = "chrUn_", # prefix for non-assembled-molecule seqs
  keep_intermediate   = TRUE      # keep cleaned FASTA + intermediate 2bit
)

## ===== 1. DEPENDENCIES =======================================================

ensure_pkgs <- function() {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  bioc_pkgs <- c("BSgenomeForge", "Biostrings", "rtracklayer",
                 "GenomeInfoDb", "GenomicRanges", "BSgenome")
  for (p in bioc_pkgs) {
    if (!requireNamespace(p, quietly = TRUE))
      BiocManager::install(p, update = FALSE, ask = FALSE)
  }
  if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")   # for build/install of the forged pkg
}

## ===== 2. HELPER FUNCTIONS ===================================================

## Assign a clean seqname from a chrom-info row.
## Uses SequenceName (always populated) rather than AssignedMolecule (sometimes NA for new genomes).
assign_chrname <- function(role, seqname, circular, chrom_prefix, unplaced_prefix) {
  role    <- as.character(role)
  seqname <- as.character(seqname)
  if (is.na(seqname)) return(NA_character_)
  sanitised <- gsub("[^A-Za-z0-9]", "_", seqname)
  if (!is.na(role) && role == "assembled-molecule") {
    ## Detect mitochondrial / plastid / circular molecules by NAME, not by
    ## AssignedMolecule. Common conventions: "MT", "MtDNA", "Mt",
    ## "mitochondrion", "chrM", "Pt" (plastid), "Chloroplast".
    if (grepl("^(MT|MtDNA|Mt|mitochondrion|chrM)$", seqname, ignore.case = TRUE))
      return(paste0(chrom_prefix, "MT"))
    if (grepl("^(Pt|Chloroplast|chrPt|plastid)$", seqname, ignore.case = TRUE))
      return(paste0(chrom_prefix, "Pt"))
    return(paste0(chrom_prefix, sanitised))
  }
  paste0(unplaced_prefix, sanitised)
}

## De-duplicate names by appending .1, .2, ... to repeated occurrences.
fix_dup_names <- function(x) {
  ave(x, x, FUN = function(v) {
    if (length(v) == 1L) v
    else c(v[1], paste0(v[-1], ".", seq_len(length(v) - 1L)))
  })
}

## Order names: chrN (numeric) first, then chrMT/chrPt, then everything else.
## Works for any chrom_prefix (including ""). Always returns an index vector
## the same length as `nms`.
preferred_order <- function(nms, chrom_prefix) {
  pfx <- if (nzchar(chrom_prefix)) paste0("^", chrom_prefix) else "^"
  chrN_pat <- paste0(pfx, "([0-9]+)$")
  mt_pat   <- paste0(pfx, "MT$")
  pt_pat   <- paste0(pfx, "Pt$")
  is_chrN  <- grepl(chrN_pat, nms)
  is_chrMT <- grepl(mt_pat,   nms)
  is_chrPt <- grepl(pt_pat,   nms)
  is_other <- !(is_chrN | is_chrMT | is_chrPt)
  chrN_nums <- suppressWarnings(
    as.integer(sub(chrN_pat, "\\1", nms[is_chrN]))
  )
  chrN_order <- order(chrN_nums)
  ord <- c(which(is_chrN)[chrN_order],
           which(is_chrMT),
           which(is_chrPt),
           which(is_other))
  if (length(ord) != length(nms))
    stop("preferred_order(): produced ", length(ord),
         " indices for ", length(nms), " names. Names: ",
         paste(head(nms, 5), collapse = ", "), " ...")
  ord
}

## ===== 3. MAIN ENTRY POINT ===================================================

run_forge <- function(assembly_accession,
                      organism,
                      genome_name,
                      pkg_maintainer,
                      workdir,
                      install_after_forge = TRUE,
                      chrom_prefix        = "chr",
                      unplaced_prefix     = "chrUn_",
                      keep_intermediate   = TRUE) {

  ensure_pkgs()
  suppressPackageStartupMessages({
    library(BSgenomeForge)
    library(Biostrings)
    library(rtracklayer)
    library(GenomeInfoDb)
    library(GenomicRanges)
    library(BSgenome)
  })

  dir.create(workdir,
             showWarnings = FALSE,
             recursive = TRUE)
  fastadir <- file.path(workdir, "fasta")
  dir.create(fastadir,
             showWarnings = FALSE,
             recursive = TRUE)
  clean_fasta <- file.path(workdir, paste0(genome_name, ".clean.fa.gz"))
  twobit_path <- file.path(workdir, paste0(genome_name, ".2bit"))

  cat("================================================================\n")
  cat(" BSgenome forge: ", genome_name, " (", assembly_accession, ")\n", sep = "")
  cat(" workdir: ", workdir, "\n", sep = "")
  cat("================================================================\n")

  ## ----- 3a. Download FASTA -----------------------------------------------
  ##
  ## NCBI's vertebrate assemblies are 400-1000+ MB; the default R timeout
  ## of 60 s is far too short. We raise it for the duration of the download
  ## and verify against NCBI's published MD5 to catch silent truncation
  ## (which previously slipped through and produced BSgenomes containing
  ## only a partial subset of chromosomes).

  DOWNLOAD_TIMEOUT_SEC <- 1800L   # This is deliberately generous; tune for network capability / genome size

  raw_fasta <- list.files(fastadir, pattern = "\\.(fa|fna|fasta)\\.gz$",
                          full.names = TRUE)

  ## Helper: fetch the NCBI md5checksums.txt for this assembly and look up
  ## the expected MD5 for the genomic FASTA filename.
  fetch_expected_md5 <- function(accn, fasta_filename) {
    ## Build the assembly's NCBI directory URL from the accession.
    ## e.g. GCA_052040795.1 -> /GCA/052/040/795/
    parts <- strsplit(sub("\\..*$", "", accn), "_")[[1]]
    prefix <- parts[1]              # GCA or GCF
    digits <- parts[2]              # e.g. 052040795
    dir1 <- substr(digits, 1, 3)
    dir2 <- substr(digits, 4, 6)
    dir3 <- substr(digits, 7, 9)
    ## We need the assembly-name suffix to find the actual subdir;
    ## NCBI's "all_assembly_versions/<accn>" symlink resolves it.
    base_url <- sprintf(
      "https://ftp.ncbi.nlm.nih.gov/genomes/all/%s/%s/%s/%s/",
      prefix, dir1, dir2, dir3
    )
    ## List directory to find the <accn>_<name> subfolder
    idx <- tryCatch(readLines(base_url, warn = FALSE),
                    error = function(e) character(0))
    subdir <- grep(paste0("\"", accn, "_"), idx, value = TRUE)
    subdir <- sub('.*href="([^"]+)/".*', "\\1", subdir)
    subdir <- subdir[nzchar(subdir)][1]
    if (is.na(subdir)) return(NA_character_)
    md5_url <- paste0(base_url, subdir, "/md5checksums.txt")
    md5_lines <- tryCatch(readLines(md5_url, warn = FALSE),
                          error = function(e) character(0))
    if (length(md5_lines) == 0) return(NA_character_)
    m <- grep(paste0("\\./", fasta_filename, "$"), md5_lines, value = TRUE)
    if (length(m) == 0)
      m <- grep(fasta_filename, md5_lines, value = TRUE, fixed = TRUE)
    if (length(m) == 0) return(NA_character_)
    sub("\\s+.*$", "", m[1])
  }

  ## Helper: compute a file's MD5 (uses base tools::md5sum)
  file_md5 <- function(path) unname(tools::md5sum(path))

  need_download <- length(raw_fasta) == 0
  if (!need_download) {
    ## We have a candidate file. Verify it before reusing.
    raw_fasta <- raw_fasta[1]
    expected <- fetch_expected_md5(assembly_accession, basename(raw_fasta))
    if (is.na(expected)) {
      warning("Could not fetch NCBI MD5 for ", basename(raw_fasta),
              "; reusing existing file without verification.")
    } else if (!identical(file_md5(raw_fasta), expected)) {
      message("Existing FASTA failed MD5 check (likely truncated from ",
              "a previous timeout). Re-downloading.")
      file.remove(raw_fasta)
      need_download <- TRUE
    } else {
      message("Reusing existing FASTA (MD5 verified): ", raw_fasta)
    }
  }

  if (need_download) {
    old_timeout <- getOption("timeout")
    on.exit(options(timeout = old_timeout), add = TRUE)
    options(timeout = DOWNLOAD_TIMEOUT_SEC)
    message("Downloading FASTA for ", assembly_accession,
            " (timeout = ", DOWNLOAD_TIMEOUT_SEC, " s) ...")
    raw_fasta <- downloadGenomicSequencesFromNCBI(
      assembly_accession = assembly_accession,
      destdir            = fastadir,
      method             = "auto",
      quiet              = FALSE
    )
    ## Post-download MD5 verification - fail loudly on truncation
    expected <- fetch_expected_md5(assembly_accession, basename(raw_fasta))
    if (is.na(expected)) {
      warning("Downloaded FASTA but could not verify against NCBI MD5. ",
              "Proceed with caution.")
    } else {
      got <- file_md5(raw_fasta)
      if (!identical(got, expected)) {
        stop("MD5 mismatch after download of ", basename(raw_fasta),
             ".\n  Expected: ", expected,
             "\n  Got:      ", got,
             "\nThe file is incomplete or corrupted. Delete it and re-run.")
      }
      message("Download MD5 verified.")
    }
  }

  ## ----- 3b. Read headers & strip to accession -----------------------------

  message("Reading FASTA ...")
  dna <- readDNAStringSet(raw_fasta)
  orig_names <- names(dna)
  accessions <- sub("\\s.*$", "", orig_names)   # first whitespace-token
  names(dna) <- accessions
  cat("Sequences in FASTA:", length(dna), "\n")

  ## ----- 3c. Fetch chrom-info; tolerate failure ----------------------------

  chrom_info <- tryCatch(
    {
      ci <- getChromInfoFromNCBI(assembly_accession)
      ci[!is.na(ci$RefSeqAccn) | !is.na(ci$GenBankAccn), , drop = FALSE]
    },
    error = function(e) {
      warning("getChromInfoFromNCBI() failed: ", conditionMessage(e),
              "\nProceeding with raw FASTA headers as seqnames.")
      NULL
    }
  )

  ## ----- 3d. Build the accession -> chrname map ----------------------------
  ## Try RefSeq accessions first (GCF_* assemblies), then GenBank.

  new_names <- accessions   # Default: keep raw accession
  circ_seqs <- character(0)
  chrom_info_used <- FALSE

  if (!is.null(chrom_info) && nrow(chrom_info) > 0) {
    chrom_info$chrname <- vapply(seq_len(nrow(chrom_info)), function(i) {
      assign_chrname(chrom_info$SequenceRole[i],
                     chrom_info$SequenceName[i],
                     if ("circular" %in% colnames(chrom_info))
                       chrom_info$circular[i] else NA,
                     chrom_prefix, unplaced_prefix)
    }, character(1))

    ## Try RefSeq first
    rs_map <- setNames(chrom_info$chrname, chrom_info$RefSeqAccn)
    rs_map <- rs_map[!is.na(names(rs_map))]
    gb_map <- setNames(chrom_info$chrname, chrom_info$GenBankAccn)
    gb_map <- gb_map[!is.na(names(gb_map))]

    n_rs <- sum(accessions %in% names(rs_map))
    n_gb <- sum(accessions %in% names(gb_map))
    which_map <- if (n_rs >= n_gb) rs_map else gb_map
    which_label <- if (n_rs >= n_gb) "RefSeq" else "GenBank"
    cat("Resolving seqnames via", which_label,
        "accessions (matched", max(n_rs, n_gb), "of", length(accessions), ")\n")

    new_names <- ifelse(accessions %in% names(which_map),
                        unname(which_map[accessions]),
                        accessions)   # fall back to raw accn
    chrom_info_used <- TRUE

    ## Circular sequences: trust the chrom-info `circular` column.
    if ("circular" %in% colnames(chrom_info)) {
      circ_rows <- chrom_info[!is.na(chrom_info$circular) &
                                chrom_info$circular, , drop = FALSE]
      circ_seqs <- unique(circ_rows$chrname)
      circ_seqs <- circ_seqs[!is.na(circ_seqs)]
    }
    if (length(circ_seqs) == 0) {
      ## Fallback: name-based detection (mitochondrial / plastid)
      mt_like <- grepl(paste0("^", chrom_prefix, "(MT|Pt)$"),
                       new_names, ignore.case = FALSE)
      circ_seqs <- unique(new_names[mt_like])
    }
  } else {
    message("No chrom-info available; using raw accessions as seqnames.")
    ## Heuristic circ detection in this fallback path
    circ_seqs <- character(0)
  }

  ## ----- 3e. De-duplicate names --------------------------------------------

  dup_count_before <- sum(duplicated(new_names))
  if (dup_count_before > 0) {
    cat("De-duplicating", dup_count_before, "name collisions ...\n")
    new_names <- fix_dup_names(new_names)
    stopifnot(!anyDuplicated(new_names))
  }
  names(dna) <- new_names
  cat("Sample of final seqnames (first 10): ",
      paste(head(new_names, 10), collapse = ", "), "\n", sep = "")

  ## ----- 3f. Reorder -------------------------------------------------------

  ord <- preferred_order(new_names, chrom_prefix)
  dna <- dna[ord]

  ## ----- 3g. Write clean FASTA --------------------------------------------

  if (file.exists(clean_fasta)) file.remove(clean_fasta)
  message("Writing cleaned FASTA: ", clean_fasta)
  writeXStringSet(dna, filepath = clean_fasta, compress = TRUE, format = "fasta")

  ## ----- 3h. Convert to 2bit (without assembly_accession!) ----------------

  if (file.exists(twobit_path)) file.remove(twobit_path)
  message("Converting FASTA -> 2bit ...")
  ## fastaTo2bit() signature is (origfile, destfile, assembly_accession=NA).
  ## We deliberately omit assembly_accession to avoid the cross-check that
  ## fails on brand-new / unregistered assemblies. The default BSgenomeForge
  ## commands should still work well with common/well-used genomes.
  fastaTo2bit(origfile = clean_fasta, destfile = twobit_path)

  ## ----- 3i. Validate 2bit -------------------------------------------------

  message("Validating 2bit ...")
  tb <- TwoBitFile(twobit_path)
  tb_seqs <- seqlevels(tb)
  stopifnot(length(tb_seqs) == length(dna))
  stopifnot(!anyDuplicated(tb_seqs))
  ## Smoke read of the first sequence
  sl <- seqlengths(tb)
  smoke_chr <- tb_seqs[1]
  smoke_len <- min(1000L, sl[[smoke_chr]])
  test_seq <- import(tb,
                     which = GRanges(smoke_chr, IRanges(1, smoke_len)))
  stopifnot(sum(width(test_seq)) == smoke_len)
  cat("2bit validated:", length(tb_seqs), "sequences,",
      round(sum(as.numeric(sl)) / 1e9, 3), "Gb total\n")

  ## Reconcile circular-sequence list with what's actually in the 2bit
  circ_seqs <- intersect(circ_seqs, tb_seqs)
  cat("Circular sequences:",
      if (length(circ_seqs)) paste(circ_seqs, collapse = ", ") else "(none)",
      "\n")

  ## ----- 3j. Move aside any old forged package directory -------------------

  old_pkgs <- list.dirs(workdir, recursive = FALSE)
  old_pkgs <- old_pkgs[grepl("^BSgenome\\.", basename(old_pkgs))]
  for (p in old_pkgs) {
    bak <- paste0(p, ".old_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    message("Moving previous package aside: ", basename(p),
            " -> ", basename(bak))
    file.rename(p, bak)
  }

  ## ----- 3k. Forge ---------------------------------------------------------

  message("Forging BSgenome data package ...")
  forgeBSgenomeDataPkgFromTwobitFile(
    filepath       = twobit_path,
    organism       = organism,
    provider       = "NCBI",
    genome         = genome_name,
    pkg_maintainer = pkg_maintainer,
    circ_seqs      = circ_seqs,
    destdir        = workdir
  )

  pkg_dirs <- list.dirs(workdir, recursive = FALSE)
  pkg_dirs <- pkg_dirs[grepl("^BSgenome\\.", basename(pkg_dirs))]
  pkg_dirs <- pkg_dirs[!grepl("\\.old_", basename(pkg_dirs))]
  if (length(pkg_dirs) == 0)
    stop("Forge produced no package directory; something went wrong.")
  pkg_dir <- pkg_dirs[which.max(file.info(pkg_dirs)$mtime)]
  pkg_name <- basename(pkg_dir)
  cat("Forged package source:", pkg_dir, "\n")

  ## ----- 3l. Optionally build + install ------------------------------------

  if (install_after_forge) {
    message("Building source tarball ...")
    tarball <- devtools::build(pkg_dir, path = workdir)
    message("Installing ", tarball, " ...")
    install.packages(tarball, repos = NULL, type = "source")

    library(pkg_name, character.only = TRUE)
    genome <- get(pkg_name)
    cat("\n--- Installed BSgenome object summary ---\n")
    print(genome)
    cat("Total length (Gb):",
        round(sum(as.numeric(seqlengths(genome))) / 1e9, 3), "\n")
  }

  ## ----- 3m. Optional cleanup ----------------------------------------------

  if (!keep_intermediate) {
    for (f in c(clean_fasta, twobit_path)) {
      if (file.exists(f)) file.remove(f)
    }
  }

  invisible(list(
    pkg_dir       = pkg_dir,
    pkg_name      = pkg_name,
    twobit_path   = twobit_path,
    clean_fasta   = clean_fasta,
    n_sequences   = length(tb_seqs),
    total_length  = sum(as.numeric(sl)),
    circ_seqs     = circ_seqs,
    chrom_info_used = chrom_info_used
  ))
}

## ===== 4. VERIFICATION HELPER ================================================
##
## A lightweight sanity check to be run after the forge. Pass it a
## data frame of (symbol, expected_chr, approx_start_bp) tuples, or just
## let it do object-level checks only.

verify_bsgenome <- function(pkg_name,
                            panel = NULL,
                            window = 2000L) {
  suppressPackageStartupMessages({
    library(BSgenome)
    library(Biostrings)
  })
  library(pkg_name, character.only = TRUE)
  genome <- get(pkg_name)

  cat("=== Object summary ===\n"); print(genome)
  sn <- seqnames(genome); sl <- seqlengths(genome)
  cat("seqs:", length(sn), " total Gb:",
      round(sum(as.numeric(sl)) / 1e9, 3),
      " circular:", paste(sn[isCircular(genome) %in% TRUE], collapse = ","),
      "\n")

  if (is.null(panel)) return(invisible(NULL))
  stopifnot(all(c("symbol", "expected_chr", "approx_start_bp")
                %in% colnames(panel)))
  out <- data.frame()
  for (i in seq_len(nrow(panel))) {
    sym <- panel$symbol[i]; chr <- panel$expected_chr[i]
    if (!(chr %in% sn)) {
      cat(sprintf("[%s] expected %s -- MISSING\n", sym, chr)); next
    }
    s <- max(1L, min(as.integer(panel$approx_start_bp[i]), sl[[chr]] - window))
    e <- s + window - 1L
    seq <- getSeq(genome, chr, start = s, end = e)
    lf  <- letterFrequency(seq, c("A","C","G","T","N"))
    gc  <- 100 * (lf["G"] + lf["C"]) / max(1, sum(lf[c("A","C","G","T")]))
    cat(sprintf("[%-8s] %s:%d-%d GC=%.1f%% N=%d\n",
                sym, chr, s, e, gc, lf["N"]))
    cat("    first60: ",
        as.character(subseq(seq, 1L, min(60L, length(seq)))), "\n", sep="")
    out <- rbind(out, data.frame(symbol=sym, chr=chr, start=s, end=e,
                                 gc=round(gc,2), n=unname(lf["N"]),
                                 first60=as.character(subseq(seq, 1L,
                                                             min(60L, length(seq)))),
                                 stringsAsFactors=FALSE))
  }
  invisible(out)
}

## ===== 5. RUN ================================================================
## Uncomment to execute with the PARAMS above. Or call run_forge() / verify_bsgenome()
## directly with your own arguments.

if (interactive() && exists("PARAMS")) {
  info <- do.call(run_forge, PARAMS)
  cat("\nForge complete. Installed package:", info$pkg_name, "\n")

  ## Optional: provide a verification panel and call verify_bsgenome(). E.g.
  ## panel <- data.frame(
  ##     symbol           = c("sox2",  "pax6a"),
  ##     expected_chr     = c("chr22", "chr25"),
  ##     approx_start_bp  = c(39480000, 8420000),
  ##     stringsAsFactors = FALSE
  ## )
  ## verify_bsgenome(info$pkg_name, panel)
}
