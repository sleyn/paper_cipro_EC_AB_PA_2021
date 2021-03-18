library(CNOGpro)
# 1 - reference GBK
# 2 - Chromosome name
# 3 - IS elements table
# 4 - Sample name
args = commandArgs(trailingOnly = T)

print(paste('CNV analysis for sample', args[4], 'chromosome', args[2]))

check_overlap <- function(region_list, reference_list){
  print('Checking overlap with IS elements')
  max = max(reference_list[,4], reference_list[,3])
  ch = rep(0,max)	# chromosome with positions occuied by IS elements marked
  for (row in 1:nrow(reference_list)){     # mark IS positions
    if( reference_list[row, 3] <= reference_list[row, 4] ){    # if left position of IS element less than right, skip it
      ch[reference_list[row, 3]:reference_list[row, 4]] = 1
    }
  }
  
  check_if_is = function(x){
    if( as.integer(x[5]) > length(ch) ){    # gene is definitly out of the IS element
      return(0)
    }
	    l = as.integer(x[4])
	    r = as.integer(x[5])
	    s = sum(ch[l:r])
	    if( s / (r - l) > 0.3 ){
	      return(1)
	    }else{
	      return(0)
	    }
  }
  
  is_checked = apply(region_list, 1, check_if_is)
  return(is_checked)
}

print('Running CNOGpro CNV identification')
cnv = CNOGpro('out.hits', args[1])
cnv_norm = normalizeGC(cnv)
cnvB = runBootstrap(cnv_norm)
significant = cnvB$genes[cnvB$genes$Type == "CDS",]
write.table(cnvB$genes, file = paste0('CNOGpro_bootstrap', args[4], '_', args[2], '.txt'), sep = '\t', row.names=F)
    
# Read IS table
if( args[3] != 'NA' ){			# if IS table is provided
  print(paste('Reading IS table', args[3]))
  IS.table = read.table(args[3], header = F, sep = '\t')
  IS.table = IS.table[IS.table[,2] == args[2],]	# filter by chromosome
  if( nrow(IS.table) > 0){
    is_rows = check_overlap(significant, IS.table)
    significant = significant[which(is_rows == 0),]
  }
}

print(paste0('Generating output: ', 'CNOGpro_', args[4], '_chr.', args[2], '.txt'))
output = data.frame(
			Clone = args[4],
			Reactor = "",
			Time = "",
			Gene_ID = significant$Locus,
			PA = "",
			AA_Mutation = "",
			In_unevolved = "",
			Annotation = "",
			MIC = "",
			Comment = "",
			X = "",
			Genome = args[2],
			DNA_Mutation = paste(significant$Left, '-', significant$Right),
			Effect = paste(significant$CN_boot)
	    	   )

write.table(output, paste0("CNOGpro_", args[4], "_chr.", args[2], ".txt"), sep = '\t', row.names = F ,quote = F)
