# arguments.
args<-commandArgs(TRUE)
# load into R and plot pathways using PathView.
library(pathview)
ko <- read.table(args[1])
rownames(ko) <- ko[,1]
ko <- ko[,2:7]
ko[,2] = log2(ko[,2])	# meth ratio -> logMethRatio
colnames(ko) <- c("logFC","logMethRatio","DiffExpUpDown","PresentInGenome","Meth","DiffMethUpDown")
map <- read.table(args[2], colClasses=c("character"))
# paint colors.
for (i in 1:nrow(map))
{
	pathview(gene.data=ko[,3:6], pathway.id=map[i,1], species="ko", out.suffix=args[3], low=list(gene="blue"), mid=list(gene="transparent"), high=list(gene="red"))
}

save.image(args[4])
