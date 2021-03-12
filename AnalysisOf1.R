#restriction of RMA to only iNS and CD34
gse44532_rma_1 <- oligo::rma(gse44532_celfiles[,1:5])


design_core1 <- model.matrix( ~ 0 + gse44532_rma_1[['cell type:ch1']])
colnames(design_core1) <- levels(as.factor(gse44532_rma_1[['cell type:ch1']]))
kable(design_core1)

contrast_matrix_core1 <- makeContrasts(iNS - CD34,
                                       levels=design_core1)
contrast_matrix_core1

fit1 <- lmFit(gse44532_rma_1,design_core1)
fit1.2 <- contrasts.fit(fit1,contrasts=contrast_matrix_core1)
fit1.2 <- eBayes(fit1.2)
summary(decideTests(fit1.2,lfc=1))

ps1 <- topTable(fit1.2, number=Inf,p.value = 0.05,lfc=2, adjust.method = 'BH')
ps1.1 <- topTable(fit1.2, number=Inf)
ps1_up <- rownames(ps1[ps1$logFC > 0,])
df1_core <- AnnotationDbi::select(hugene20sttranscriptcluster.db,ps1_up,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
dplyr::mutate(df1_core,GENENAME=stringr::str_trunc(GENENAME,30))

#visualization
#volcano plots CORE
volcanoplot(fit1.2, coef=1)
#this one had a lot of features pass the cutoff
interesting_genes_core1 <- topTable(fit1.2, number=Inf,p.value = 0.05,lfc=2)
volcanoplot(fit1.2, coef=1, main=sprintf("%d features pass our cutoffs",nrow(interesting_genes_core1)))
points(interesting_genes_core1[['logFC']],-log10(interesting_genes_core1[['P.Value']]),col='red')

eset_of_interest_core1 <- gse44532_rma_1[rownames(interesting_genes_core1),]
heatmap(exprs(eset_of_interest_core1))
heatmap(exprs(eset_of_interest_core1),
        labCol=gse44532_rma_1[['cell type:ch1']] ,labRow=NA,
        col       = rev(brewer.pal(10, "RdBu")),
        distfun   = function(x) as.dist(1-cor(t(x))))


#determining ontology
write.table(df1_core$ENTREZID, "iNS-CD34.txt")