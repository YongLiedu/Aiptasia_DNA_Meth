# arguments.
args<-commandArgs(TRUE)
# load into R and plot pathways using PathView.
library(pathview)
ko <- read.table(args[1])
rownames(ko) <- ko[,1]
ko <- ko[,2:3]
colnames(ko) <- c("logFCup","logFCdn")
map <- read.table(args[2], colClasses=c("character"))

# paint colors.
for (i in 1:nrow(map))
{
	pathview(gene.data=ko, pathway.id=map[i,1], species="ko", out.suffix=args[3], low=list(gene="blue"), mid=list(gene="transparent"), high=list(gene="red"))
}
save.image(args[4])
