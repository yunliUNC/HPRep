## run example:
##  Rscript normalize_bins.R /proj/yunligrp/users/minghu/armen/HP_reprod/humanBrain_YinLab_hg38/cell_specific/IJ050_merged_trimmedto100_fastp_current/ 
##  reg_raw IJ050_merged_trimmedto100_fastp.5k 5000 22 None pospoisson
##
## arguments:
## INFDIR - dir with reg files
## PREFIX
## SUFFIX
## RESOLUTION - resolution (for example 5000 or 10000)
## chroms - number of shromosomes (19 for mouse, 22 for human)
## FILTER - file containing bins that need to be filtered out. Format: two columns "chrom", "bin". "chrom" contains 'chr1','chr2',.. "bin" is bin label
## regresison_type - pospoisson for positive poisson regression, negbinom for negative binomial. default is pospoisson

library(VGAM)
library(MASS)
library(data.table)
options(warn=-1)

### constants
chroms = NULL
runs = c(1)
RESOLUTION = NULL
###

args <- commandArgs(trailingOnly=TRUE)
fltr = data.frame(chr='chrNONE',bin=-1)

if (length(args) != 8) {
    print('Wrong number of arguments. Stopping.')
    print('Arguments needed (in this order): INFDIR, PREFIX, SUFFIX, RESOLUTION, chroms, FILTER, regression_type, OUTDIR.')
    print(paste('Number of arguments entered:',length(args)))
    print('Arguments entered:')
    print(args)
    quit()
} else {
    print(args)
    INFDIR = args[1]
    PREFIX = args[2]
    SUFFIX = args[3]
    RESOLUTION = as.integer(args[4])
    chroms = paste0('chr', seq(1, as.numeric(args[5]), 1))
    if (args[6] != 'None') {
        FILTER = args[6]
        fltr = read.table(FILTER,header=T)
    }
    REG_TYPE = args[7]
    OUTDIR = args[8]
}

print('filter used (if any):')
if (args[6] == 'None') { print('None')} else { print(fltr) }

#COUNT_CUTOFF = 12
#RATIO_CUTOFF = 2.0
#GAP = 15000

## loading data
mm_combined_and = data.frame()
mm_combined_xor = data.frame()
outf_names = c()
for (i in chroms) {
    for (j in c('.and','.xor')) {
        print(paste0('loading chromosome ',i,' ',j))
        inf_name = paste0(INFDIR, PREFIX, '.', i, '.', SUFFIX, j)
        #outf_names = c(outf_names, paste('regression.cell_specific.',i,'.',SET,j,'.MAPS2_',REG_TYPE,sep = ''))
	#mm = read.table(inf_name, header=T)
        mm = fread(inf_name, data.table = F)
        mm$chr = i
	mm = subset(mm, dist > 1) # removing adjacent bins
        mm = subset(mm, count != 0) # rmoving bins with zero count - only factors in downsampling
	mm = subset(mm, !(mm$chr %in% fltr$chr & (mm$bin1_mid %in% fltr$bin | mm$bin2_mid %in% fltr$bin ))) ## filtering out bad bins
        if (j == '.and') {
            mm_combined_and = rbind(mm_combined_and, mm)
        } else if (j == '.xor') {
            mm_combined_xor = rbind(mm_combined_xor, mm)
        }
    }
}

dataset_length_and = length(mm_combined_and$bin1_mid)
dataset_length_xor = length(mm_combined_xor$bin1_mid)
dataset_length = dataset_length_and + dataset_length_xor

## doing statistics and resampling

pospoisson_regression <- function(mm) {
    #fit <- vglm(count ~ logl + loggc + logm + logdist + logShortCount, family = pospoisson(), data = mm)
    fit <- vglm(count ~ logl + loggc + logm + logShortCount, family = pospoisson(), data = mm) # logdist removed for HP reproducibility 
    mm$expected = fitted(fit)
    mm$p_val = ppois(mm$count, mm$expected, lower.tail = FALSE, log.p = FALSE) / ppois(0, mm$expected, lower.tail = FALSE, log.p = FALSE)
    return(mm)
}

negbinom_regression <- function(mm) {
    #fit <- glm.nb(count ~ logl + loggc + logm + logdist + logShortCount, data = mm)
    fit <- glm.nb(count ~ logl + loggc + logm + logShortCount, data = mm) # logdist removed for HP reproducibility
    mm$expected = fitted(fit)
    sze = fit$theta ##size parameter
    mm$p_val = pnbinom(mm$count, mu = mm$expected, size = sze, lower.tail = FALSE)
    return(mm)
}


mx_combined_and = data.frame()
mx_combined_xor = data.frame()

for (r in runs) {
## in case you want to do resampling
## you'd put resampling code here
    #name_counter = 1
    for (i in chroms) {
        ## regression
        print(paste('run',r,': regression on chromosome',i))
        #print(outf_names[name_counter])
        mm = subset(mm_combined_and, chr == i)
        if (REG_TYPE == 'pospoisson') {
            mm = pospoisson_regression(mm)
        } else if (REG_TYPE == 'negbinom') {
            mm = negbinom_regression(mm)
        }
        #write.table(mm,outf_names[name_counter],row.names = TRUE,col.names = TRUE,quote=FALSE)
        #name_counter = name_counter + 1
        #print(outf_names[name_counter])
        mx_combined_and = rbind(mx_combined_and, mm)
        
        mm = subset(mm_combined_xor, chr == i)
        if (REG_TYPE == 'pospoisson') {
            mm = pospoisson_regression(mm)
        } else if (REG_TYPE == 'negbinom') {
            mm = negbinom_regression(mm)
        }
        #write.table(mm,outf_names[name_counter],row.names = TRUE,col.names = TRUE,quote=FALSE)
        #name_counter = name_counter + 1
        mx_combined_xor = rbind(mx_combined_xor, mm)
    }
}

mx_combined_all = rbind(mx_combined_and, mx_combined_xor)
mx_combined_all$norm = log(1 + mx_combined_all$count / mx_combined_all$expected, 2)

res = mx_combined_all[, c("chr", "bin1_mid", "bin2_mid", "count", "expected", "norm")]
colnames(res) = c("Chr", "Bin1", "Bin2", "Obs_count", "Exp_count", "Normalized_count") 

options(scipen = 999)
write.table(res, file = paste0(OUTDIR, SUFFIX, '.normalized.txt'), row.names = F, col.names = T, quote = F)



