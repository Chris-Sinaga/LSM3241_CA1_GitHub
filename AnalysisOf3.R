#restriction of RMA to only NPC and CD34
gse44532_rma_3 <- oligo::rma(gse44532_celfiles[, c(1,2,8:11)])


design_core3 <- model.matrix( ~ 0 + gse44532_rma_3[['cell type:ch1']])
colnames(design_core3) <- levels(as.factor(gse44532_rma_3[['cell type:ch1']]))
kable(design_core3)

contrast_matrix_core3 <- makeContrasts(NPC - CD34,
                                       levels=design_core3)
contrast_matrix_core3

fit3 <- lmFit(gse44532_rma_3,design_core3)
fit3.2 <- contrasts.fit(fit3,contrasts=contrast_matrix_core3)
fit3.2 <- eBayes(fit3.2)
summary(decideTests(fit3.2,lfc=1))

ps3 <- topTable(fit3.2, number=Inf,p.value = 0.05,lfc=2, adjust.method = 'BH')
ps3.1 <- topTable(fit3.2,number = Inf)
ps3_up <- rownames(ps3[ps3$logFC > 0,])
df3_core <- AnnotationDbi::select(hugene20sttranscriptcluster.db,ps3_up,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
dplyr::mutate(df3_core,GENENAME=stringr::str_trunc(GENENAME,30))

#visualization
volcanoplot(fit3.2, coef=1)

interesting_genes_core3 <- topTable(fit3.2, number=Inf,p.value = 0.05,lfc=2)
volcanoplot(fit3.2, coef=1, main=sprintf("%d features pass our cutoffs",nrow(interesting_genes_core3)))
points(interesting_genes_core3[['logFC']],-log10(interesting_genes_core3[['P.Value']]),col='red')

eset_of_interest_core3 <- gse44532_rma_3[rownames(interesting_genes_core3),]
heatmap(exprs(eset_of_interest_core3))
heatmap(exprs(eset_of_interest_core3),
        labCol=gse44532_rma_3[['cell type:ch1']] ,labRow=NA,
        col       = rev(brewer.pal(10, "RdBu")),
        distfun   = function(x) as.dist(1-cor(t(x))))

#determining ontology
write.table(df3_core$ENTREZID, "NPC-CD34.txt")
