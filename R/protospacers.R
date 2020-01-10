#' find_proto
#'
#' Function to identify all candidate guide (protospacer) sequences from a given genomic region.

#' @rdname find_proto
#' @param d_seq input DNA sequence for genomic region from which to identify protospacer sequences (in capital letters).
#' @param l protospacer length.
#' @param PAM sequence to match.
#' @return a dataframe with columns.
#' \enumerate{
#'     \item  start.
#'     \item  end.
#'     \item  protospacer sequence.
#'     \item  PAM sequence (e.g. CGG as the sequence that matched the PAM sequence NGG).
#'     \item  Strand (+ or -).
#'     }
#' @details For reverse strand reverseComplement PAM, comolement d_seq, match after PAM
#' @export
find_proto <- function(d_seq = "TGATCTACTAGAGACTACTAACGGGGATACATAG", l = 20, PAM = "NGG") {

  # arguments to uppercase
  argum <- lapply(list(PAM = PAM, d_seq = d_seq), toupper)
  list2env(argum, envir =  environment())

  # PAMs for forward and reverse strand ----
  # substitute N with [ACGT]
  # reverse complement of PAM
  PAM_r <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(PAM)))
  PAMs <- lapply(list(PAM, PAM_r), function(x) gsub("N", "[ACGT]", x))
  names(PAMs) <- c("forward", "reverse")

  # dynamically create the pattern to match ----
  pro_pattern <- list(
    # protospacer pattern forward strand
    # example: ([ACGT])(?=[ACGT]{19}[ACGT]GG)
    # looks for START of protospacer: a character that is followed by l -  1 characters and PAM.
    forward = paste0("([ACGT])(?=[ACGT]{", l - 1, "}", PAMs$forward, ")", collapse = ""),

    # protospacer pattern reverse strand
    # example: (?<=CC[ACGT][ACGT]{19})([ACGT])
    # looks for END of protospacer:  that is a character preceded  by l - 1  characters and reverse complement PAM.
    reverse = paste0("(?<=", PAMs$reverse, "[ACGT]", "{", l - 1, "}", ")", "([ACGT])", collapse = "")
  )

  # extract the protospacers ----
  # Note that we could actually create a function to apply to both strands, and thus
  # eliminate some duplication

  # forward strand
  # locate the start of the protospacers in forward
  pro_start_fw <- as.data.frame(stringr::str_locate_all(d_seq, pro_pattern$forward))

  # extract protospacers in forward
  protospacers_fw <- purrr::map_dfr(pro_start_fw$start, function(x) {
    start_p <- x
    end_p <- x + l - 1
    protospacer <- stringr::str_sub(d_seq, start = start_p, end = end_p)
    PAM <- stringr::str_sub(d_seq, start = end_p + 1, end = end_p + 3)
    strand <- "+"

    return(list(
      start_p = start_p,
      end_p = end_p,
      protospacer = protospacer,
      PAM = PAM,
      strand = strand
    ))
  })

  # reverse strand
  # locate the END of the protospacers in reverse strand

  # we need the complement of d_seq
  d_seq_rev <- as.character(Biostrings::complement(Biostrings::DNAString(d_seq)))
  pro_end_rev <- as.data.frame(stringr::str_locate_all(d_seq_rev, pro_pattern$reverse))

  protospacers_neg <- purrr::map_dfr(pro_end_rev$end, function(x) {
    start_p <- x - l + 1
    end_p <- x
    protospacer <- stringr::str_sub(d_seq_rev, start = start_p, end = end_p)
    PAM <- stringr::str_sub(d_seq_rev, start = start_p - nchar(PAM) + 1, end = start_p)
    strand <- "-"

    return(list(
      start = start_p,
      end = end_p,
      protospacer = protospacer,
      PAM = PAM,
      strand = strand
    ))
  })

  # bind them together
  rbind(protospacers_fw, protospacers_neg)
}

#' find_proto2
#'
#' Function to identify all candidate guide (protospacer) sequences from a given genomic region.
#' @details For reverse strand: complement PAM, complement d_seq, match before PAM
#' @export
#' @rdname find_proto
find_proto2 <- function(d_seq = "TGATCTACTAGAGACTACTAACGGGGATACATAG", l = 20, PAM = "NGG") {

  # arguments to uppercase
  argum <- lapply(list(PAM = PAM, d_seq = d_seq), toupper)
  list2env(argum, envir =  environment())

  # PAMs for forward and reverse strand ----
  # substitute N with [ACGT]
  # reverse complement of PAM
  PAM_r <- as.character(Biostrings::complement(Biostrings::DNAString(PAM)))
  PAMs <- lapply(list(PAM, PAM_r), function(x) gsub("N", "[ACGT]", x))
  names(PAMs) <- c("forward", "reverse")

  # dynamically create the pattern to match ----
  # protospacer pattern forward strand
  # example: ([ACGT])(?=[ACGT]{19}[ACGT]GG)
  # looks for START of protospacer: a character that is followed by l -  1 characters and PAM.
  pro_pattern <-     paste0("([ACGT])(?=[ACGT]{", l - 1, "}", PAMs$forward, ")", collapse = "")

  # extract the protospacers ----
  # Note that we could actually create a function to apply to both strands, and thus
  # eliminate some duplication

  # locate the start of the protospacers ---
  d_seq_rev <- as.character(Biostrings::complement(Biostrings::DNAString(d_seq)))

  pro_start <- purrr::map2_dfr(list(d_seq, d_seq_rev), list("+", "-"), function(x, y){
                      pro_start <- as.data.frame(stringr::str_locate_all(x, pro_pattern))
                      pro_start[,"strand"] <- y
                      pro_start
  })

  # extract protospacers ---
  # note: needs d_seq and d_seq_rev!
  protospacers <- purrr::map2_dfr(pro_start$start, pro_start$strand, function(x, y) {
    if(y == "+") {
      sequen <- d_seq
    } else {
      sequen <- d_seq_rev
    }
    start_p <- x
    end_p <- x + l - 1
    protospacer <- stringr::str_sub(sequence, start = start_p, end = end_p)
    PAM <- stringr::str_sub(d_seq, start = end_p + 1, end = end_p + 3)
    strand <- y

    return(list(
      start_p = start_p,
      end_p = end_p,
      protospacer = protospacer,
      PAM = PAM,
      strand = strand
    ))
  })

  protospacers

}




#' find_FASTA
#'
#' Function to identify all candidate guide (protospacer) sequences in a FASTA file. Note that the genomic
#' coordinates are \strong{1-indexed} and \strong{fully closed}.
#'
#' @param file_fasta  path to a file FASTA (either compressed or not).
#' @param chr chromosome.
#' @param start start of the DNA sequence to scan for protospacers.
#' @param end  end of the DNA sequence to scan for protospacers.
#' @param l protospacer length.
#' @param PAM PAM sequence to match.
#' @return a dataframe with columns:
#' \enumerate{
#'     \item  chr.
#'     \item  start.
#'     \item  end.
#'     \item  protospacer sequence.
#'     \item  PAM sequence (e.g. CGG as the sequence that matched the PAM sequence NGG).
#'     \item  Strand (+ or -).
#' }
#' @export

find_FASTA <- function(file_fasta, chr = 7, start = 117465784, end = 117466784, l = 20, PAM = "NGG") {
  browser()
  message(">>> Importing Data...\n")
  gen <- Biostrings::readDNAStringSet(file_fasta, format = "fasta")
  message(">>> Data Imported!\n")

  # example to match: "NC_000001.11 Homo sapiens chromosome 1, GRCh38.p13 Primary Assembly"
  # pattern = "chromosome 1,.*Primary Assembly$"
  chr_pattern <- paste0("chromosome ", chr, ",.*Primary Assembly$", collapse = "")

  # subset by chromosome and coordinates
  gen_chr <- gen[grepl(chr_pattern, gen@ranges@NAMES)]
  gen_chr_sub <- as.character(XVector::subseq(gen_chr, start = start, end = end))

  protospacers <- find_proto(d_seq = gen_chr_sub, l = l, PAM = PAM)

  solution <- cbind(rep(chr, nrow(protospacers)), protospacers)

  solution
}
