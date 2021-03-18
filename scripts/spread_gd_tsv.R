library(tidyverse)

# 1 - gdtools TSV output
# 2 - output TSV

args = commandArgs(trailingOnly = T)
in_tbl = args[1]
out_tbl = args[2]

comparison_tbl = read_tsv(in_tbl)
samples = comparison_tbl %>% select(title) %>% unique()
if(TRUE %in% str_detect(colnames(comparison_tbl), 'repeat')){
  comparison_tbl = comparison_tbl %>% 
    mutate(presence = '+') %>% 
    select(aa_new_seq, aa_position, aa_ref_seq, codon_new_seq, codon_number, codon_ref_seq, gene_name, gene_position, gene_product, gene_strand, locus_tag, mutation_category, new_seq, position, seq_id, title, type, repeat_new_copies, repeat_ref_copies, repeat_seq, presence) %>% 
    spread(key = title, value = presence)
  comparison_tbl = comparison_tbl %>% mutate(aa_mutation = str_c(aa_ref_seq, aa_position, aa_new_seq)) %>% 
    select(-c(aa_ref_seq, aa_position, aa_new_seq))
  comparison_tbl = comparison_tbl %>% mutate(codon_mutation = str_c(codon_ref_seq, codon_number, codon_new_seq)) %>% 
    select(-c(codon_ref_seq, codon_number, codon_new_seq))
  comparison_tbl = comparison_tbl %>% mutate(repeat_mutation = str_c(repeat_seq, ' (', repeat_ref_copies, '->', repeat_new_copies, ')')) %>% 
      select(-c(repeat_seq, repeat_ref_copies, repeat_new_copies))
  comparison_tbl[is.na(comparison_tbl)] = '-'
  comparison_tbl %>%
    select(mutation_category, type, seq_id, position, aa_mutation, codon_mutation, repeat_mutation, locus_tag, gene_name, gene_position, new_seq, gene_product, gene_strand, samples$title) %>%
    write_tsv(out_tbl)
}else{
  comparison_tbl = comparison_tbl %>% 
    mutate(presence = '+') %>% 
    select(aa_new_seq, aa_position, aa_ref_seq, codon_new_seq, codon_number, codon_ref_seq, gene_name, gene_position, gene_product, gene_strand, locus_tag, mutation_category, new_seq, position, seq_id, title, type, presence) %>% 
    spread(key = title, value = presence)
  comparison_tbl = comparison_tbl %>% mutate(aa_mutation = str_c(aa_ref_seq, aa_position, aa_new_seq)) %>% 
    select(-c(aa_ref_seq, aa_position, aa_new_seq))
  comparison_tbl = comparison_tbl %>% mutate(codon_mutation = str_c(codon_ref_seq, codon_number, codon_new_seq)) %>% 
    select(-c(codon_ref_seq, codon_number, codon_new_seq))
  comparison_tbl[is.na(comparison_tbl)] = '-'
  comparison_tbl %>%
    select(mutation_category, type, seq_id, position, aa_mutation, codon_mutation, locus_tag, gene_name, gene_position, new_seq, gene_product, gene_strand, samples$title) %>%
    write_tsv(out_tbl)
}
