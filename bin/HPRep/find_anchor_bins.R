## run example:
##  Rscript find_anchor_bins.R /proj/yunligrp/users/minghu/armen/HP_reprod/mouse_brain/anchors/P56_H3K4me3_peaks.bed 
##
## arguments:
## INFILE
## Resolution
## chroms
## OUTFILE


args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
    print('Wrong number of arguments. Stopping.')
    print('Arguments needed (in this order): INFILE, Resolution, chroms, OUTFILE')
    print(paste('Number of arguments entered:',length(args)))
    print('Arguments entered:')
    print(args)
    quit()
} else {
    print(args)
    INFILE = args[1]
    Resolution = as.numeric(args[2])
    chroms = as.numeric(args[3])
    OUTFILE = args[4]
}

x = read.table(INFILE)[, 1:3]

# BinSize
Bin.Size = Resolution

# find unique ChIP-Seq peaks
x = unique(x)
x = as.data.frame(x)

colnames(x) = c('chr','start','end')
#dim(x)

chrnum = chroms
x = x[ x[,1] %in% paste0('chr', seq(1, chrnum)), ]
#dim(x)

y = x
#dim(y)
#summary(y[,3]-y[,2])

# find bin contains ChIP-Seq peak
y$start.bin = floor(y$start/Bin.Size)*Bin.Size
y$end.bin = floor(y$end/Bin.Size)*Bin.Size

u = y[y$start.bin == y$end.bin, ]
v = y[y$start.bin != y$end.bin, ]
dim(u) # ChIP-Seq peak within one bin
dim(v) # ChIP-Seq peak within multiple bins

# break down multiple bins
w = NULL
for(chrid in 1:chrnum) {
	vsub = v[ v[,1] == paste0('chr', chrid) ,]
	tmp = NULL
	for(i in 1:nrow(vsub)) {
		tmp = c(tmp, seq(vsub$start.bin[i], vsub$end.bin[i], by=Bin.Size) )
	}
 	out = cbind(rep(chrid, length(tmp)), tmp)
	out = as.data.frame(out)
	colnames(out) = c('chr', 'Bin')
	out[,1] = paste0('chr', out[,1])
	w = rbind(w, out)
}
dim(w)

# ChIP-Seq peak within one bin, keep chr id, bin start position
u2 = u[, c(1,4)]
colnames(u2) = c('chr', 'Bin')

# combine ChIP-Seq peak within one 5Kb bin and ChIP-Seq peak within multiple 5Kb bins
# only keep unique bins
x = unique( rbind(u2, w) )
dim(x)

y = NULL
for(chrid in 1:chrnum) {
	xsub = x[ x[,1] == paste0('chr', chrid), ]
 	y = rbind(y, xsub[ order(xsub$Bin), ])
}
dim(y)

colnames(y) = c('chr', 'start')
y$end = y$start + Bin.Size


# output bed format
options(scipen=999)
write.table(y, file=paste0(OUTFILE, '.txt'), row.names = F, col.names = T, sep='\t', quote=F)




