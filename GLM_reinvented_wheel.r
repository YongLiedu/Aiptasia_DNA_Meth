unlogit <- function(x)
{
    exp(x) / (1 + exp(x))
}

load_meth_table <- function(path)
{
    
    d          <- read.delim(path, header=FALSE)
    names(d)   <- c("scaffold", "position", "meth", "unmeth", "gene", "treatment")
    
    # Sanity checks
    #MINMETHN   <- 2
    #MINMETHP   <- 0.1     at least one (or all?!) treatment should have proportion > 0.1
    #MINCOV     <- 4
    #MAXCOV     <- 100
    #d          <- d[(d$qm >= MINMETHN) | (d$wm >= MINMETHN), ]
    #d          <- d[(d$qm + d$qu > MINCOV) & (d$wm + d$wu > MINCOV), ]
    #d          <- d[(d$qm + d$qu < MAXCOV) & (d$wm + d$wu < MAXCOV), ]
    #d          <- d[(d$qm / (d$qm + d$qu) > MINMETHP) | (d$wm / (d$wm + d$wu) > MINMETHP), ]
    #d$position <- d$position + 1
    
    d
}

do_GLM_gene <- function(meth_table, min_cpg=2)
{
    # define constants
    MINMETHP   <- 0.1
    MAXCOV     <- 100
    TREATMENTS <- unique(meth_table$treatment)
    REPLICATES <- 6
    
    p.value.treatment <- c()
    p.value.inter     <- c()
    n_cpgs_orig       <- c()
    n_cpgs            <- c()
    treatment.ratio   <- c()
    converged         <- c()
    unique_genes      <- unique(meth_table$gene)
    
    counter <- 0
    for (g in unique_genes)
    {
        # progress bar
        counter = counter + 1
        cat (counter, '/', length(unique_genes), '; ', g, '\n')
        
        # get rows corresponding to gene g
        cpgs        <- meth_table[meth_table$gene == g, ]
        n_cpgs_orig <- c(n_cpgs_orig, length(unique(cpgs$position)))
        
        # sanity checks:
        # 1. retain positions where at least one treatment has all six 
        #    replicates with methylation ratio of > MINMETHP.
        # 2. retain positions where ALL replicates has coverage less than
        #    or equals MAXCOV.
        sane_rows <- c()
        for (p in unique(cpgs$position))
        {
            sliced_cpgs <- cpgs[cpgs$position == p, ]
            if (nrow(sliced_cpgs[sliced_cpgs$meth + sliced_cpgs$unmeth <= MAXCOV, ]) < length(TREATMENTS) * REPLICATES)
            {
                # criteria 2: skip position if not all replicates are under
                # coverage cutoff
                next
            }
            
            for (q in TREATMENTS)
            {
                if (nrow(sliced_cpgs[sliced_cpgs$treatment == q & sliced_cpgs$meth / (sliced_cpgs$meth + sliced_cpgs$unmeth) > MINMETHP, ]) >= REPLICATES)
                {
                    # at least a treatment has all replicates with MINMETHP
                    sane_rows <- c(sane_rows, row.names(cpgs[cpgs$position == p, ]))
                    break
                }
            }
        }
        cpgs <- cpgs[sane_rows, ]
        n_cpgs <- c(n_cpgs, length(unique(cpgs$position)))
        
        if (length(unique(cpgs$position)) < min_cpg)
        {
            p.value.treatment  <- c(p.value.treatment, 1)
            p.value.inter      <- c(p.value.inter,     1)
            treatment.ratio    <- c(treatment.ratio,   1)
            converged          <- c(converged,         0)
            next
        }
        
        # temp arrays that will be fed into GLM
        treat <- c()
        pos   <- c()
        freq  <- c()
        me    <- c()
        un    <- c()
        for (j in 1:nrow(cpgs))
        {
            rr    <- cpgs[j, ]
            pos   <- c(pos, rr$position)
            me    <- c(me, rr$meth)
            un    <- c(un, rr$unmeth)
            treat <- c(treat, rr$treatment)
        }
        treat <- factor(treat)
        pos   <- factor(pos)
        contrasts(treat) <- "contr.sum"
        contrasts(pos)   <- "contr.sum"
        d  <- data.frame(meth=me,
                         unmeth=un,
                         treatment=treat,
                         position=pos)
        l <- glm(matrix(c(meth, unmeth), ncol=2) ~ treatment * position, 
                 data=d, family=binomial)
        s <- step(l, trace=0)
        a <- anova(s, test="Chisq")
        t <- labels(terms(s))
        if (s$converged)
        {
            if ("treatment" %in% t)
            {
                p.value.treatment <- c(p.value.treatment, a["treatment", "Pr(>Chi)"])
                tmp.coef          <- coef(s)
                tmp.ratio         <- unlogit(tmp.coef[['(Intercept)']] + tmp.coef[['treatment1']]) / unlogit(tmp.coef[['(Intercept)']] - tmp.coef[['treatment1']])
                treatment.ratio   <- c(treatment.ratio, tmp.ratio)
            }
            else
            {
                p.value.treatment <- c(p.value.treatment, 1)
                treatment.ratio   <- c(treatment.ratio,   1)
            }
            
            if ("treatment:position" %in% t)
            {
                p.value.inter <- c(p.value.inter, a["treatment:position", "Pr(>Chi)"])
            }
            else
            {
                p.value.inter <- c(p.value.inter, 1)
            }
            converged <- c(converged, 1)
        }
        else
        {
            treatment.ratio   <- c(treatment.ratio,   1)
            converged         <- c(converged,         0)
            p.value.treatment <- c(p.value.treatment, 1)
            p.value.inter     <- c(p.value.inter,     1)
        }
    }
    genes <- data.frame(gene=as.vector(unique_genes),
                        ncpgs_orig=n_cpgs_orig,
                        ncpgs=n_cpgs,
                        treatment.ratio=treatment.ratio,
                        converged=converged,
                        p.value.treatment=p.value.treatment,
                        p.value.inter=p.value.inter)
    genes <- genes[genes$ncpgs >= min_cpg, ]
    genes <- cbind(genes,
                   adj.p.value.treatment=p.adjust(genes$p.value.treatment, method='BH'),
                   adj.p.value.inter=p.adjust(genes$p.value.inter, method='BH'))
    genes
}

meth_table  <- load_meth_table('compiled_meth_unmeth_gene_treatment.genic.tsv')

sub_meth_table <- meth_table[meth_table$treatment == 'apo' | meth_table$treatment == 'symb', ]
genes <- do_GLM_gene(sub_meth_table, min_cpg=5)
write.table(genes, 'reinv.apovsymb.5cpg.tsv', sep='\t', row.names=FALSE)
