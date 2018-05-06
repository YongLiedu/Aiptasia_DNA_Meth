#' topGO enrichment analysis.
#' input format.
#' universal: one gene per line with GO numbers.
#'     gene_01 GO:0003674,GO:0003824,GO:0004129
#'     gene_02 GO:0005575,GO:0006810,GO:0008150,GO:0009987
#' interest genes: one gene name per line.
#'     gene_01
#'     gene_02

library(topGO)

working_folder <- "."
mult_files = list.files(working_folder, pattern = "Meth.*?.txt")
annot_filename = paste0(working_folder, "/", "aip_go_annots.all.tsv")

for (go_category in c("bp", "cc", "mf")) {
    gene_id_to_go = readMappings(file = annot_filename)
    gene_id_to_go = gene_id_to_go[gene_id_to_go != "no_hit"]
    gene_names = names(gene_id_to_go)

    for (m in mult_files) {
        print(paste("Current file:", m))
        genes_of_interest_filename = paste0(working_folder, "/", m)
        genes_of_interest = scan(genes_of_interest_filename, character(0), sep = "\n")

        # hack: spis/smic genes have spaces in them, but gene_id_to_go has had spaces
        # stripped from them - so we've got to strip spaces from genes_of_interest as
        # well.
        genes_of_interest = gsub(" ", "", genes_of_interest)

        genelist = factor(as.integer(gene_names %in% genes_of_interest), levels = c(0,
            1))
        names(genelist) = gene_names

        GOdata = try(new("topGOdata", ontology = toupper(go_category), allGenes = genelist,
            gene2GO = gene_id_to_go, annotationFun = annFUN.gene2GO))

        # handle error
        if (class(GOdata) == "try-error") {
            print(paste0("Error for file", m, "!"))
            next
        }

        # weight01 is the default algorithm used in Alexa et al. (2006)
        weight01.fisher <- runTest(GOdata, statistic = "fisher")

        # generate a results table (for only the top 1000 GO terms) topNodes: highest
        # 1000 GO terms shown numChar: truncates GO term descriptions at 1000 chars
        # (basically, disables truncation)
        results_table = GenTable(GOdata, P_value = weight01.fisher, orderBy = "P_value",
            topNodes = 100, numChar = 100)

        # write it out into a file for python post-processing
        output_filename = paste0(working_folder, "/", go_category, "_", m)
        write.table(results_table, file = output_filename, quote = FALSE, sep = "\t")
    }
}
