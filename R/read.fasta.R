# canprot/read.fasta.R
# Read FASTA sequence files to get amino acid compositions of proteins
# 20240227 Moved from CHNOSZ

read.fasta <- function(file, iseq = NULL, ret = "count", lines = NULL, ihead = NULL,
  start = NULL, stop = NULL, type = "protein", id = NULL) {
  # Read sequences from a fasta file
  # Some of the following code was adapted from 
  # read.fasta in package seqinR
  # value of 'iseq' is what sequences to read (default is all)
  # value of 'ret' determines format of return value:
  #   count: amino acid composition (output can be used by CHNOSZ::add.protein())
  #        or nucleic acid base composition (A-C-G-T)
  #   seq: amino acid sequence
  #   fas: fasta entry
  # value of 'id' is used for 'protein' in output table,
  #   otherwise ID is parsed from FASTA header (can take a while)
  
  # Check if the file is in an archive (https://github.com/jimhester/archive)
  if(inherits(file, "archive_read")) {
    is.archive <- TRUE
    filebase <- gsub("]", "", basename(summary(file)$description))
  } else {
    is.archive <- FALSE
    filebase <- basename(file)
  }
  if(is.null(lines)) {
    message("read.fasta: reading ", filebase, " ... ", appendLF = FALSE)
    is.nix <- Sys.info()[[1]] == "Linux"
    if(is.archive) {
      # We can't use scan here?
      lines <- readLines(file)
    } else if(is.nix) {
      # Retrieve contents using system command (seems slightly faster even than scan())
      # Figure out whether to use 'cat', 'zcat' or 'xzcat'
      suffix <- substr(file,nchar(file)-2,nchar(file))
      if(suffix == ".gz") mycat <- "zcat"
      else if(suffix == ".xz") mycat <- "xzcat"
      else mycat <- "cat"
      lines <- system(paste(mycat,' "',file,'"',sep = ""),intern = TRUE)
    } else lines <- scan(file, what = character(), sep = "\n", quiet = TRUE)
  }
  nlines <- length(lines)
  message(nlines, " lines ... ", appendLF = FALSE)
  if(is.null(ihead)) ihead <- which(substr(lines, 1, 1) == ">")
  message(length(ihead), " sequences")
  # Identify the lines that begin and end each sequence
  begin <- ihead + 1
  end <- ihead - 1
  end <- c(end[-1], nlines)
  # Use all or selected sequences
  if(is.null(iseq)) iseq <- seq_along(begin)
  # Just return the lines from the file
  if(ret == "fas") {
    iline <- numeric()
    for(i in iseq) iline <- c(iline, (begin[i]-1):end[i])
    return(lines[iline])
  }
  # Get each sequence from the begin to end lines
  seqfun <- function(i) paste(lines[begin[i]:end[i]], collapse = "")
  sequences <- lapply(iseq, seqfun)
  # Organism name is from file name
  # (basename minus extension)
  bnf <- strsplit(filebase, split = ".", fixed = TRUE)[[1]][1]
  organism <- bnf
  # Protein/gene name is from header line for entry
  # (strip the ">" and go to the first space)
  missid <- missing(id)
  if(is.null(id)) id <- as.character(lapply(iseq, function(j) {
    # Get the text of the line
    f1 <- lines[ihead[j]]
    # Stop if the first character is not ">"
    # or the first two charaters are "> "
    if(substr(f1, 1, 1) != ">" | length(grep("^> ", f1)>0))
      stop(paste("file", filebase, "line", j, "doesn't begin with FASTA header '>'."))
    # Discard the leading '>'
    f2 <- substr(f1,  2,  nchar(f1))
    # Keep everything before the first space
    return(strsplit(f2, " ")[[1]][1])
  } ))
  if(ret == "count") {
    counts <- count.aa(sequences, start, stop, type)
    ref <- abbrv <- NA
    chains <- 1
    if(type == "protein") {
      colnames(counts) <- c("Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu",
                            "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr")
      # 20090507 Made stringsAsFactors FALSE
      out <- cbind(protein = id, organism = organism,
        ref = ref, abbrv = abbrv, chains = chains, stringsAsFactors = FALSE, counts)
      # 20170117 Extra processing for files from UniProt
      isUniProt <- grepl("\\|......\\|.*_", out$protein[1])
      if(isUniProt & missid) {
        p1 <- sapply(strsplit(out$protein, "\\|"), "[", 1)
        p2 <- sapply(strsplit(out$protein, "\\|"), "[", 2)
        p3 <- sapply(strsplit(out$protein, "\\|"), "[", 3)
        out$abbrv <- sapply(strsplit(p3, "_"), "[", 1)
        out$organism <- sapply(strsplit(p3, "_"), "[", 2)
        out$protein <- paste0(p1, "|", p2)
      }
      out
    } else if(type %in% c("DNA", "RNA")) {
      cbind(gene = id, organism = organism,
        ref = ref, abbrv = abbrv, chains = chains, stringsAsFactors = FALSE, counts)
    }
  } else return(sequences)
}

# Count amino acids or DNA bases in one or more sequences
# 20090423 Use table(strsplit(toupper(seq), "")[[1]]) for counting (CHNOSZ 0.8)
# 20180517 Use C code for counting (CHNOSZ 1.2.0)
# 20240229 Use stringi::stri_count_fixed for counting (canprot 1.1.2-22)
count.aa <- function(sequence, start = NULL, stop = NULL, type = "protein") {
  if(type == "protein") alphabet <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y") else
  if(type == "DNA") alphabet <- c("A", "C", "G", "T") else
  if(type == "RNA") alphabet <- c("A", "C", "G", "U") else
  stop(paste("unknown sequence type", type))
  # Convert input to uppercase
  sequence <- toupper(sequence)
  # Loop over sequences
  counts <- lapply(sequence, function(seq) {
    # Get a substring between start and stop positions
    # If only one of start or stop is given, get a default value for the other
    if(!is.null(start)) {
      if(is.null(stop)) stop <- nchar(seq)
      seq <- substr(seq, start, stop)
    } else if(!is.null(stop)) {
      seq <- substr(seq, 1, stop)
    }
    # Count all letters
    stri_count_fixed(seq, LETTERS)
  })
  # Convert to data frame
  counts <- as.data.frame(do.call(rbind, counts))
  # Keep only letters for protein, DNA, or RNA
  iab <- match(alphabet, LETTERS)
  ina <- colSums(counts[, -iab, drop = FALSE]) > 0
  if(any(ina)) message(paste("count.aa: unrecognized letter(s) in", type, "sequence:", paste(LETTERS[-iab][ina], collapse = " ")))
  counts <- counts[, iab, drop = FALSE]
  # Add column and row names
  colnames(counts) <- alphabet
  if(is.null(names(sequence))) rownames(counts) <- 1:nrow(counts)
  counts
}

# Sum or average amino acid counts (weighted by abundance if given)
aasum <- function(AAcomp, abundance = 1, average = FALSE) {
  # Recycle abundance values into same length as number of proteins
  abundance <- rep(abundance, length.out = nrow(AAcomp))
  # Find columns with names for the amino acids and "chains"
  AA_names <- c(
    "Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu",
    "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr",
    "chains"
  )
  isAA <- tolower(colnames(AAcomp)) %in% tolower(AA_names)
  # Multiply amino acid counts by protein abundance
  AAcomp[, isAA] <- AAcomp[, isAA] * abundance
  # Sum amino acid counts
  out <- AAcomp[1, ]
  out[, isAA] <- colSums(AAcomp[, isAA], na.rm = TRUE)
  # Take the average
  if(average) out[, isAA] <- out[, isAA] / sum(abundance)
  out
}
