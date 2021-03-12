#For iNS vs CD34
#restriction of RMA to only iNS(Sample2 only) and CD34 (both same hybridization batch)
gse44532_rma_1_cov1.1 <- oligo::rma(gse44532_celfiles[,1:4])

design_core1_cov1.1 <- model.matrix( ~ 0 + gse44532_rma_1_cov1.1[['cell type:ch1']])
colnames(design_core1_cov1.1) <- levels(as.factor(gse44532_rma_1_cov1.1[['cell type:ch1']]))
kable(design_core1_cov1.1)

contrast_matrix_core1_cov1.1 <- makeContrasts(iNS - CD34,
                                              levels=design_core1_cov1.1)
contrast_matrix_core1_cov1.1

fit1_cov1.1 <- lmFit(gse44532_rma_1_cov1.1,design_core1_cov1.1)
fit1.2_cov1.1 <- contrasts.fit(fit1_cov1.1,contrasts=contrast_matrix_core1_cov1.1)
fit1.2_cov1.1 <- eBayes(fit1.2_cov1.1)
summary(decideTests(fit1.2_cov1.1,lfc=1))

ps1_cov1.1 <- topTable(fit1.2_cov1.1, number=Inf,p.value = 0.05,lfc=1)
ps1_cov1.1_up <- rownames(ps1_cov1.1[ps1_cov1.1$logFC > 0,])
df1_core_cov1.1 <- AnnotationDbi::select(hugene20sttranscriptcluster.db,ps1_cov1.1_up,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
dplyr::mutate(df1_core_cov1.1,GENENAME=stringr::str_trunc(GENENAME,30))

#visualization
#volcano plots CORE
volcanoplot(fit1.2_cov1.1, coef=1)
#this one had a lot of features pass the cutoff
interesting_genes_core1_cov1.1 <- topTable(fit1.2_cov1.1, number=Inf,p.value = 0.05,lfc=2)
volcanoplot(fit1.2_cov1.1, coef=1, main=sprintf("%d features pass our cutoffs",nrow(interesting_genes_core1_cov1.1)))
points(interesting_genes_core1_cov1.1[['logFC']],-log10(interesting_genes_core1_cov1.1[['P.Value']]),col='red')

eset_of_interest_core1_cov1.1 <- gse44532_rma_1_cov1.1[rownames(interesting_genes_core1_cov1.1),]
heatmap(exprs(eset_of_interest_core1_cov1.1))
heatmap(exprs(eset_of_interest_core1_cov1.1),
        labCol=gse44532_rma_1_cov1.1[['cell type:ch1']] ,labRow=NA,
        col       = rev(brewer.pal(10, "RdBu")),
        distfun   = function(x) as.dist(1-cor(t(x))))

#restriction of RMA to only iNS(Sample3) and CD34 (different hybridization batch)
gse44532_rma_1_cov1.2 <- oligo::rma(gse44532_celfiles[,c(1,2,5)])

design_core1_cov1.2 <- model.matrix( ~ 0 + gse44532_rma_1_cov1.2[['cell type:ch1']])
colnames(design_core1_cov1.2) <- levels(as.factor(gse44532_rma_1_cov1.2[['cell type:ch1']]))
kable(design_core1_cov1.2)

contrast_matrix_core1_cov1.2 <- makeContrasts(iNS - CD34,
                                              levels=design_core1_cov1.2)
contrast_matrix_core1_cov1.2

fit1_cov1.2 <- lmFit(gse44532_rma_1_cov1.2,design_core1_cov1.2)
fit1.2_cov1.2 <- contrasts.fit(fit1_cov1.2,contrasts=contrast_matrix_core1_cov1.2)
fit1.2_cov1.2 <- eBayes(fit1.2_cov1.2)
summary(decideTests(fit1.2_cov1.2,lfc=1))

ps1_cov1.2 <- topTable(fit1.2_cov1.2, number=Inf,p.value = 0.05,lfc=1)
ps1_cov1.2_up <- rownames(ps1_cov1.2[ps1_cov1.2$logFC > 0,])
df1_core_cov1.2 <- AnnotationDbi::select(hugene20sttranscriptcluster.db,ps1_cov1.2_up,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
dplyr::mutate(df1_core_cov1.2,GENENAME=stringr::str_trunc(GENENAME,30))

#visualization
#volcano plots CORE
volcanoplot(fit1.2_cov1.2, coef=1)
#this one had a lot of features pass the cutoff
interesting_genes_core1_cov1.2 <- topTable(fit1.2_cov1.2, number=Inf,p.value = 0.05,lfc=2)
volcanoplot(fit1.2_cov1.2, coef=1, main=sprintf("%d features pass our cutoffs",nrow(interesting_genes_core1_cov1.2)))
points(interesting_genes_core1_cov1.2[['logFC']],-log10(interesting_genes_core1_cov1.2[['P.Value']]),col='red')

eset_of_interest_core1_cov1.2 <- gse44532_rma_1_cov1.2[rownames(interesting_genes_core1_cov1.2),]
heatmap(exprs(eset_of_interest_core1_cov1.2))
heatmap(exprs(eset_of_interest_core1_cov1.2),
        labCol=gse44532_rma_1_cov1.2[['cell type:ch1']] ,labRow=NA,
        col       = rev(brewer.pal(10, "RdBu")),
        distfun   = function(x) as.dist(1-cor(t(x))))


#compare with original
summary(decideTests(fit1.2,lfc=1))
dplyr::mutate(df1_core,GENENAME=stringr::str_trunc(GENENAME,30))
volcanoplot(fit1.2, coef=1, main=sprintf("%d features pass our cutoffs",nrow(interesting_genes_core1)))
points(interesting_genes_core1[['logFC']],-log10(interesting_genes_core1[['P.Value']]),col='red')

heatmap(exprs(eset_of_interest_core1))
heatmap(exprs(eset_of_interest_core1),
        labCol=gse44532_rma_1[['cell type:ch1']] ,labRow=NA,
        col       = rev(brewer.pal(10, "RdBu")),
        distfun   = function(x) as.dist(1-cor(t(x))))
