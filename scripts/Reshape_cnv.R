library(tidyverse)
# 1 - Directory with unformatted CNOGpro output
# 2 - Output directory
# 3 - Threshold to keed observation
args = commandArgs(trailingOnly = T)

cnog_dir = args[1]
out_dir = args[2]
thrsh = args[3]

cnog_files = Sys.glob(file.path(cnog_dir, 'CNOGpro_*.txt'))
cnv = cnog_files %>% map(~ read_tsv(.)) %>% reduce(rbind)
# cnv = cnv %>% mutate(Effect = as.numeric(str_match(Effect, 'CNV  (\\d+\\.\\d+)')[,2]))
cnv = cnv %>% select(Gene_ID, Genome, DNA_Mutation, Clone, Effect) %>% spread(key = Clone, value = Effect)
cnv = cnv %>% mutate_at(vars(-Gene_ID, -Genome, -DNA_Mutation), ~ replace(., is.na(.), 1))

cnv$max_deivation = apply(cnv %>% select(c(-Gene_ID, -Genome, -DNA_Mutation)), 1, function(x) max(abs(x - 1)))
cnv %>% filter(max_deivation >= thrsh) %>% write_tsv(file.path(out_dir, 'CNV_comparison.txt'), quote_escape = F)
