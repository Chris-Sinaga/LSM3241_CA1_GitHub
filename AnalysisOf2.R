#restriction of RMA to only iNS derived from iPSC and CD34
gse44532_rma_2 <- oligo::rma(gse44532_celfiles[, c(1,2,6,7)])


design_core2 <- model.matrix( ~ 0 + gse44532_rma_2[['cell type:ch1']])
colnames(design_core2) <- levels(as.factor(gse44532_rma_2[['cell type:ch1']]))

# Change colname of iNS derived from iPSC to iNS_derived_from_iPSC for downstream analysis
colnames(design_core2)[colnames(design_core2) == "iNS derived from iPSC"] <- "iNS_derived_from_iPSC"

kable(design_core2)

contrast_matrix_core2 <- makeContrasts(iNS_derived_from_iPSC - CD34,
                                       levels=design_core2)
contrast_matrix_core2

fit2 <- lmFit(gse44532_rma_2,design_core2)
fit2.2 <- contrasts.fit(fit2,contrasts=contrast_matrix_core2)
fit2.2 <- eBayes(fit2.2)
summary(decideTests(fit2.2,lfc=1))

ps2 <- topTable(fit2.2, number=Inf,p.value = 0.05,lfc=2, adjust.method = 'BH')
ps2.1 <- topTable(fit2.2, number=Inf)
ps2_up <- rownames(ps2[ps2$logFC > 0,])
df2_core <- AnnotationDbi::select(hugene20sttranscriptcluster.db,ps2_up,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
dplyr::mutate(df2_core,GENENAME=stringr::str_trunc(GENENAME,30))

#visualization
volcanoplot(fit2.2, coef=1)

interesting_genes_core2 <- topTable(fit2.2, number=Inf,p.value = 0.05,lfc=2)
volcanoplot(fit2.2, coef=1, main=sprintf("%d features pass our cutoffs",nrow(interesting_genes_core2)))
points(interesting_genes_core2[['logFC']],-log10(interesting_genes_core2[['P.Value']]),col='red')

eset_of_interest_core2 <- gse44532_rma_2[rownames(interesting_genes_core2),]
heatmap(exprs(eset_of_interest_core2))
heatmap(exprs(eset_of_interest_core2),
        labCol=gse44532_rma_2[['cell type:ch1']] ,labRow=NA,
        col       = rev(brewer.pal(10, "RdBu")),
        distfun   = function(x) as.dist(1-cor(t(x))))

#determining ontology
write.table(df2_core$ENTREZID, "iPSCiNS-CD34.txt")
