knitr::opts_chunk$set(echo = TRUE, cache=TRUE, cache.path = "cache/CountKmers-")
library(tidyverse)
library(Biostrings)

## Define kmer-counting functions
make_kmer_list <- function(k=2L, alphabet = c("A","C","G","T")) {
  list(alphabet) %>%
    base::rep(times = k) %>%
    purrr::set_names(LETTERS[1:k]) %>%
    base::expand.grid() %>%
    tidyr::unite("kmer",LETTERS[1:k],sep="") %>%
    dplyr::pull(kmer)
}

count_kmer_one <- function(string,kmer_list,k=2L) {
  stopifnot( nchar(kmer_list[1]) == k )
  
  kCount <- rep(0L,length(kmer_list)) %>% purrr::set_names(kmer_list)
  kminus1 <- k - 1L
  
  for(i in 1:( nchar(string) - kminus1 ) ) {
    kmer <- stringr::str_sub(string, start=i, end = i + kminus1)
    if(kmer %in% kmer_list) {
      kCount[kmer] <- kCount[kmer] + 1L
    }
  }
  return(kCount)
}

count_kmer_tibble <- function(stringset,kmer_list,k=2L) {
  lapply(stringset,count_kmer_one,kmer_list=kmer_list,k=k) %>%
    lapply(as.list) %>%
    dplyr::bind_rows()
}


promoterseqs_file <- here::here("H99_allorfs_p500.fasta")
promoterseqs <- Biostrings::readDNAStringSet(promoterseqs_file)

## 3-mers
promoterseqids <- names(promoterseqs) %>%
  stringr::str_extract(pattern="\\w+")

all_3mers <- make_kmer_list(k=3L)

promoter_3mer_countsonly <- promoterseqs %>%
  as.character() %>%
  count_kmer_tibble(kmer_list = all_3mers, k=3L)


promoter_3mer_counts <- bind_cols(tibble(Gene=promoterseqids), promoter_3mer_countsonly)
write.table(promoter_3mer_counts, file= "promoter_3mer_counts.txt", sep=",", quote=FALSE, row.names = F)


## 4-mers
promoterseqids <- names(promoterseqs) %>%
  stringr::str_extract(pattern="\\w+")

all_4mers <- make_kmer_list(k=4L)

promoter_4mer_countsonly <- promoterseqs %>%
  as.character() %>%
  count_kmer_tibble(kmer_list = all_4mers, k=4L)


promoter_4mer_counts <- bind_cols(tibble(Gene=promoterseqids), promoter_4mer_countsonly)
write.table(promoter_4mer_counts, file= "promoter_4mer_counts.txt", sep=",", quote=FALSE, row.names = F)

## 6-mers
promoterseqids <- names(promoterseqs) %>%
  stringr::str_extract(pattern="\\w+")

all_6mers <- make_kmer_list(k=6L)

promoter_6mer_countsonly <- promoterseqs %>%
  as.character() %>%
  count_kmer_tibble(kmer_list = all_6mers, k=6L)


promoter_6mer_counts <- bind_cols(tibble(Gene=promoterseqids), promoter_6mer_countsonly)

write.table(promoter_6mer_counts, file= "promoter_6mer_counts.txt", sep=",", quote=FALSE, row.names = F)