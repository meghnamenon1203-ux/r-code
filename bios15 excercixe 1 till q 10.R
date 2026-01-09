library(seqinr)
mito_fasta <- read.fasta('NC_012920.fasta')
mito_seq <- mito_fasta$NC_012920.1
nucleo_content <- function(seq){
  seq <- tolower(seq)
  g_count <- sum(seq == 'g')
  c_count <- sum(seq == 'c')
  a_count <- sum(seq == 'a')
  t_count <- sum(seq == 't')
  total <- length(seq)
  g_percent <- (g_count / total) * 100
  c_percent <- (c_count / total) * 100
  a_percent <- (a_count / total) * 100
  t_percent <- (t_count / total) * 100
  
  cat(sprintf("G: %.3f%%\n", g_percent))
  cat(sprintf("C: %.3f%%\n", c_percent))
  cat(sprintf("A: %.3f%%\n", a_percent))
  cat(sprintf("T: %.3f%%\n", t_percent))
}
nucleo_content(mito_seq)
create_complement <- function(seq) {
  seq <- tolower(seq)
  comp_seq <- vector(mode = "character", length = length(seq))
  comp_seq[seq == "a"] <- "t"
  comp_seq[seq == "t"] <- "a"
  comp_seq[seq == "c"] <- "g"
  comp_seq[seq == "g"] <- "c"
  return(comp_seq)
}
#create_complement(mito_seq)
#nucleo_content(create_complement(mito_seq))
# here the nucleotide proportion remains same as gcat percentages got reversed
find_codon <- function(seq, codon, frame) {
seq_str <- paste(seq, collapse = "") #Convert vector of single nucleotides to one string
codon <- tolower(codon) # Convert codon to lowercase for case-insensitive matching
seq_str <- tolower(seq_str)# Convert sequence to lowercase for case-insensitive matching

starts <- seq(frame, nchar(seq_str) - 2, by = 3)# Generate starting positions in the reading frame(like 1,4,7 etc)

codons <- substring(seq_str, starts, starts + 2)#Extract triplets (codons) from sequence at each start position

return(starts[codons == codon])# Return positions where codon matches
}
find_codon( c('a','c','g','g','a','t'), 'gat', 1)

find_start_codon <- function(seq,frame) {
  start_codons <- c("atg", "ata", "att", "atc", "gtg")# vertebrate mitochondrial start codons
  positions <- c()#Initializes an empty vector to store the start positions of all found start codons.
  for (codon in start_codons) {
    positions <- c(positions, find_codon(seq, codon, frame))
  }#Loops through each codon in the start codons list.Calls the previously defined function find_codon that returns all positions of the current codon in the DNA sequence seq considering the reading frame frame.Combines  all found codon positions into one vector, accumulating results from all start codons.
  
  
  return(sort(positions))
}#Returns all found positions sorted in ascending order, so start codons are reported from sequence beginning to end.
find_start_codon(mito_seq[1:500], 1)

find_stop_codon <- function(seq, frame) {
  stop_codons <- c("taa", "tag", "aga", "agg")# vertebrate mitochondrial stop codons
  positions <- c()
  for (codon in stop_codons) {
    positions <- c(positions, find_codon(seq, codon, frame))
  }
  return(sort(positions))
}
find_stop_codon(mito_seq[1:500],1)
mito_seq[13:201]  
