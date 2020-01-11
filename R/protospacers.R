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

  # if the sum of characters that is different than ACGTNacgtn is more than 0, then some wrong nucleotides
  if (sum(unlist(lapply(list(PAM, d_seq), function(x) stringr::str_count(x, "[^ACGTNacgtn]")))) > 0){
    stop("Error: unrecognized nucleotides")
  }

  if ( l == 0 || !is.numeric(l)) stop("Error: double check your numeric argument!")

  # arguments to uppercase
  argum <- lapply(list(PAM = PAM, d_seq = d_seq), toupper)
  list2env(argum, envir =  environment())

  # substitute N with [ACGT]
  PAM <- gsub("N", "[ACGT]", PAM)

  # dynamically create the pattern to match ----
  # protospacer pattern forward strand
  # example: ([ACGT])(?=[ACGT]{19}[ACGT]GG)
  # looks for START of protospacer: a character that is followed by l -  1 characters and PAM.
  pro_pattern <- paste0("([ACGT])(?=[ACGT]{", l - 1, "}", PAM, ")", collapse = "")

  # locate the start of the protospacers ---
  pro_start <- as.data.frame(stringr::str_locate_all(d_seq, pro_pattern))

  # extract protospacers ---
  # note: needs d_seq
  protospacers <- purrr::map_dfr(pro_start$start, function(x) {
    start_p <- x
    end_p <- x + l - 1
    protospacer <- stringr::str_sub(d_seq , start = start_p, end = end_p)
    PAM <- stringr::str_sub(d_seq , start = end_p + 1, end = end_p + 3)

    return(list(
      start_p = start_p,
      end_p = end_p,
      protospacer = protospacer,
      PAM = PAM
    ))
  })

  protospacers[, "strand"] <-  "+"
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

  # check arguments ---
  if (start == end || end - start <= l || !is.numeric(c(start, end))) stop("Error: double check your numeric parameters!")

  message(">>> Importing Data: this will take few seconds...\n")
  gen <- Biostrings::readDNAStringSet(file_fasta, format = "fasta")
  message(">>> Data Imported!\n")

  # example to match: "NC_000001.11 Homo sapiens chromosome 1, GRCh38.p13 Primary Assembly"
  # pattern = "chromosome 1,.*Primary Assembly$"
  chr_pattern <- paste0("chromosome ", chr, ",.*Primary Assembly$", collapse = "")

  # subset by chromosome and coordinates
  gen_chr <- gen[grepl(chr_pattern, gen@ranges@NAMES)]
  gen_chr_sub <- as.character(XVector::subseq(gen_chr, start = start, end = end))

  protospacers <- find_proto(d_seq = gen_chr_sub, l = l, PAM = PAM)

  protospacers
}
