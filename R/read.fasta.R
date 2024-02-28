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
  linefun <- function(i1, i2) lines[i1:i2]
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
  seqfun <- function(i) paste(linefun(begin[i], end[i]), collapse = "")
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
    f1 <- linefun(ihead[j], ihead[j])
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
      out <- cbind(data.frame(protein = id, organism = organism,
        ref = ref, abbrv = abbrv, chains = chains, stringsAsFactors = FALSE), counts)
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
      cbind(data.frame(gene = id, organism = organism,
        ref = ref, abbrv = abbrv, chains = chains, stringsAsFactors = FALSE), counts)
    }
  } else return(sequences)
}

count.aa <- function(seq, start = NULL, stop = NULL, type = "protein") {
  # Count amino acids or DNA bases in one or more sequences given as elements of the list seq
  if(type == "protein") letts <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  else if(type == "DNA") letts <- c("A", "C", "G", "T")
  else if(type == "RNA") letts <- c("A", "C", "G", "U")
  else stop(paste("unknown sequence type", type))
  # The numerical positions of the letters in alphabetical order
  ilett <- match(letts, LETTERS)
  # The letters A-Z represented by raw values
  rawAZ <- charToRaw("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
  # To count the letters in each sequence
  countfun <- function(seq, start, stop) {
    # Get a substring if one or both of start or stop are given
    # If only one of start or stop is given, get a default value for the other
    if(!is.null(start)) {
      if(is.null(stop)) stop <- nchar(seq)
      seq <- substr(seq, start, stop)
    } else if(!is.null(stop)) {
      seq <- substr(seq, 1, stop)
    }
    ## The actual counting ...
    #nnn <- table(strsplit(toupper(seq), "")[[1]])
    # ... Replaced with C version 20180217
    counts <- .C(C_count_letters, seq, integer(26))[[2]]
    # which is equivalent to this R code:
    #rawseq <- charToRaw(toupper(seq))
    #counts <- sapply(rawAZ, function(x) sum(rawseq == x))
    return(counts)
  }
  # Counts for each sequence
  counts <- lapply(seq, countfun, start, stop)
  counts <- do.call(rbind, counts)
  # Check for letters that aren't in our alphabet
  ina <- colSums(counts[, -ilett, drop = FALSE]) > 0
  if(any(ina)) {
    message(paste("count.aa: unrecognized letter(s) in", type, "sequence:", paste(LETTERS[-ilett][ina], collapse = " ")))
  }
  counts <- counts[, ilett, drop = FALSE]
  # Clean up row/column names
  colnames(counts) <- letts
  rownames(counts) <- 1:nrow(counts)
  return(counts)
}

# Combine amino acid counts (sum, average, or weighted sum by abundance)
aasum <- function(aa, abundance = 1, average = FALSE, protein = NULL, organism = NULL) {
  # Returns the sum of the amino acid counts in aa,
  #   multiplied by the abundances of the proteins
  abundance <- rep(abundance, length.out = nrow(aa))
  # Drop any NA rows or abundances
  ina.aa <- is.na(aa$chains)
  ina.ab <- is.na(abundance)
  ina <- ina.aa | ina.ab
  if(any(ina)) {
    aa <- aa[!ina, ]
    abundance <- abundance[!ina]
    message("aasum: dropped ", sum(ina), " proteins with NA composition and/or abundance")
  }
  # Multiply
  aa[, 6:25] <- aa[, 6:25] * abundance
  # Sum
  out <- aa[1, ]
  out[, 5:25] <- colSums(aa[, 5:25])
  # Average if told to do so
  if(average) {
    # Polypeptide chains by number of proteins, residues by frequency
    out[, 5] <- out[, 5]/nrow(aa)
    out[, 6:25] <- out[, 6:25]/sum(abundance)
  }
  # Add protein and organism names if given
  if(!is.null(protein)) out$protein <- protein
  if(!is.null(organism)) out$organism <- organism
  return(out)
}

