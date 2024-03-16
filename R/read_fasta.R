# canprot/read_fasta.R
# Read FASTA sequence files to get amino acid compositions of proteins
# 20240227 Moved from CHNOSZ

read_fasta <- function(file, iseq = NULL, type = "count", lines = NULL, ihead = NULL,
  start = NULL, stop = NULL, molecule = "protein", id = NULL) {
  # Read sequences from a fasta file
  # Some of the following code was adapted from 
  # read.fasta in package seqinR
  # value of 'iseq' is what sequences to read (default is all)
  # value of 'type' determines format of return value:
  #   count: amino acid composition (output can be used by CHNOSZ::add.protein())
  #        or nucleic acid base composition (A-C-G-T)
  #   sequence: amino acid sequences
  #   lines: lines from the file (including headers)
  #   headers: header lines only
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
    message("read_fasta: reading ", filebase, " ... ", appendLF = FALSE)
    if(is.archive) {
      # We can't use scan here?
      lines <- readLines(file)
#    } else if(Sys.info()[[1]] == "Linux") {
#      # Retrieve contents using system command (seems slightly faster even than scan())
#      # Figure out whether to use 'cat', 'zcat' or 'xzcat'
#      suffix <- substr(file, nchar(file) - 2, nchar(file))
#      if(suffix == ".gz") mycat <- "zcat" else
#      if(suffix == ".xz") mycat <- "xzcat" else mycat <- "cat"
#      lines <- system(paste(mycat, ' "', file, '"', sep = ""), intern = TRUE)
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
  if(tolower(substr(type, 1, 3)) == "lin") {
    iline <- numeric()
    for(i in iseq) iline <- c(iline, (begin[i] - 1):end[i])
    return(lines[iline])
  }
  # Just return the headers 20240304
  if(tolower(substr(type, 1, 3)) == "hea") {
    return(lines[ihead][iseq])
  }
  # Get each sequence from the begin to end lines
  seqfun <- function(i) paste(lines[begin[i]:end[i]], collapse = "")
  sequences <- lapply(iseq, seqfun)
  # Organism name is from file name (basename minus extension)
  organism <- file_path_sans_ext(filebase)
  # Protein/gene name is from header line for entry
  # (strip the ">" and go to the first space)
  missid <- missing(id)
  if(is.null(id)) id <- as.character(lapply(iseq, function(j) {
    # Get the text of the line
    f1 <- lines[ihead[j]]
    # Stop if the first character is not ">"
    # or the first two charaters are "> "
    if(substr(f1, 1, 1) != ">" | length(grep("^> ", f1) > 0))
      stop(paste("file", filebase, "line", j, "doesn't begin with FASTA header '>'."))
    # Discard the leading '>'
    f2 <- substr(f1,  2,  nchar(f1))
    # Keep everything before the first space
    return(strsplit(f2, " ")[[1]][1])
  } ))
  if(tolower(substr(type, 1, 3)) == "cou") {
    
    if(molecule == "protein") {
      if(length(iseq) == 0) {
        # Deal with length-0 iseq 20240308
        out <- structure(list(protein = character(0), organism = character(0), 
          ref = character(0), abbrv = character(0), chains = integer(0), 
          Ala = numeric(0), Cys = numeric(0), Asp = numeric(0), Glu = numeric(0), 
          Phe = numeric(0), Gly = numeric(0), His = numeric(0), Ile = numeric(0), 
          Lys = numeric(0), Leu = numeric(0), Met = numeric(0), Asn = numeric(0), 
          Pro = numeric(0), Gln = numeric(0), Arg = numeric(0), Ser = numeric(0), 
          Thr = numeric(0), Val = numeric(0), Trp = numeric(0), Tyr = numeric(0)), row.names = integer(0), class = "data.frame")
      } else {
        # Count the number of occurrences of each amino acid
        counts <- count_aa(sequences, start, stop, molecule)
        colnames(counts) <- c("Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu",
                              "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr")
        out <- cbind(protein = id, organism = organism, ref = NA, abbrv = NA, chains = 1, counts)
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
      }
    } else if(molecule %in% c("DNA", "RNA")) {
      if(length(iseq) == 0) {
        out <- structure(list(gene = character(0), organism = character(0),
          ref = character(0), abbrv = character(0), chains = integer(0),
          A = numeric(0), C = numeric(0), G = numeric(0), T = numeric(0)),
          row.names = integer(0), class = "data.frame")
        if(molecule == "RNA") colnames(out)[9] <- "U"
      } else {
        # Count the number of occurrences of each base
        counts <- count_aa(sequences, start, stop, molecule)
        out <- cbind(gene = id, organism = organism, ref = NA, abbrv = NA, chains = 1, counts)
      }
    }
    out
  } else return(sequences)
}

# Count amino acids or DNA bases in one or more sequences
# 20090423 Use table(strsplit(toupper(seq), "")[[1]]) for counting (CHNOSZ 0.8)
# 20180517 Use C code for counting (CHNOSZ 1.2.0)
# 20240229 Use stringi::stri_count_fixed for counting (canprot 1.1.2-22)
count_aa <- function(sequence, start = NULL, stop = NULL, molecule = "protein") {
  if(molecule == "protein") alphabet <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y") else
  if(molecule == "DNA") alphabet <- c("A", "C", "G", "T") else
  if(molecule == "RNA") alphabet <- c("A", "C", "G", "U") else
  stop(paste0("unknown molecule '", molecule, "'"))
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
  if(any(ina)) message(paste("count_aa: unrecognized letter(s) in", molecule, "sequence:", paste(LETTERS[-iab][ina], collapse = " ")))
  counts <- counts[, iab, drop = FALSE]
  # Add column and row names
  colnames(counts) <- alphabet
  if(is.null(names(sequence))) rownames(counts) <- 1:nrow(counts)
  counts
}

# Sum or average amino acid counts (weighted by abundance if given)
sum_aa <- function(AAcomp, abundance = 1, average = FALSE) {
  # Recycle abundance values into same length as number of proteins
  abundance <- rep(abundance, length.out = nrow(AAcomp))
  # Find columns with names for the amino acids and 'chains'
  # (number of polypeptide chains)
  AA_names <- c(
    "Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu",
    "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr",
    "chains"
  )
  isAA <- tolower(colnames(AAcomp)) %in% tolower(AA_names)
  # Remove proteins with NA amino acid composition or abundance
  ina.aa <- is.na(rowSums(AAcomp[, isAA]))
  ina.ab <- is.na(abundance)
  ina <- ina.aa | ina.ab
  if (any(ina)) {
      AAcomp <- AAcomp[!ina, ]
      abundance <- abundance[!ina]
      message("sum_aa: dropped ", sum(ina), " proteins with NA amino acid composition and/or abundance")
  }
  # Multiply amino acid counts by protein abundance
  AAcomp[, isAA] <- AAcomp[, isAA] * abundance
  # Sum amino acid counts
  out <- AAcomp[1, ]
  out[, isAA] <- colSums(AAcomp[, isAA], na.rm = TRUE)
  # Take the average
  if(average) out[, isAA] <- out[, isAA] / sum(abundance)
  out
}
