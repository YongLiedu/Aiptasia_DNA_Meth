# arguments.
args<-commandArgs(TRUE)
# fisher exact test.
fisher_test <- read.table(args[1], sep="\t", quote="")
colnames(fisher_test) <- c("PathwayID", "Pathway_name", "mapped_allMappedGenes", "mapped_allMappedDEGs", "allMappedDEGs", "allMappedGenes")
fisher_p <- c()
for (i in 1:nrow(fisher_test)){
	fisher_p <- c(fisher_p, fisher.test(matrix(c(fisher_test[i,"mapped_allMappedDEGs"],fisher_test[i,"allMappedDEGs"]-fisher_test[i,"mapped_allMappedDEGs"],fisher_test[i,"mapped_allMappedGenes"]-fisher_test[i,"mapped_allMappedDEGs"],fisher_test[i,"allMappedGenes"]-fisher_test[i,"mapped_allMappedGenes"]-fisher_test[i,"allMappedDEGs"]+fisher_test[i,"mapped_allMappedDEGs"]), ncol=2))$p.value)
}
fisher_test$fisher_p <- fisher_p
fisher_test$fdr <- p.adjust(fisher_test$fisher_p, method="BH")
fisher_test <- fisher_test[order(fisher_test$fdr,fisher_test$fisher_p),]
write.table(fisher_test, file=args[2], sep="\t", quote=F, row.names=F, col.names=T)
save.image(args[3])
