#For Core-level output, NPC vs CD34 [for diff hybridization batch]
#restriction of RMA to only NPC(Hybridization Batch 1) and CD34
gse44532_rma_3_cov4 <- oligo::rma(gse44532_celfiles[,c(1,2,8,9)])

design_core3_cov4 <- model.matrix( ~ 0 + gse44532_rma_3_cov4[['cell type:ch1']])
colnames(design_core3_cov4) <- levels(as.factor(gse44532_rma_3_cov4[['cell type:ch1']]))
kable(design_core3_cov4)

contrast_matrix_core3_cov4 <- makeContrasts(NPC - CD34,
                                            levels=design_core3_cov4)
contrast_matrix_core3_cov4

fit3_cov4 <- lmFit(gse44532_rma_3_cov4,design_core3_cov4)
fit3.2_cov4 <- contrasts.fit(fit3_cov4,contrasts=contrast_matrix_core3_cov4)
fit3.2_cov4 <- eBayes(fit3.2_cov4)
summary(decideTests(fit3.2_cov4,lfc=1))

ps3_cov4 <- topTable(fit3.2_cov4, number=Inf,p.value = 0.05,lfc=1)
ps3_cov4_up <- rownames(ps3_cov4[ps3_cov4$logFC > 0,])
df3_core_cov4 <- AnnotationDbi::select(hugene20sttranscriptcluster.db,ps3_cov4_up,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
dplyr::mutate(df3_core_cov4,GENENAME=stringr::str_trunc(GENENAME,30))

#visualization
#volcano plots CORE
volcanoplot(fit3.2_cov4, coef=1)
#this one had a lot of features pass the cutoff
interesting_genes_core3_cov4 <- topTable(fit3.2_cov4, number=Inf,p.value = 0.05,lfc=2)
volcanoplot(fit3.2_cov4, coef=1, main=sprintf("%d features pass our cutoffs",nrow(interesting_genes_core3_cov4)))
points(interesting_genes_core3_cov4[['logFC']],-log10(interesting_genes_core3_cov4[['P.Value']]),col='red')

eset_of_interest_core3_cov4 <- gse44532_rma_3_cov4[rownames(interesting_genes_core3_cov4),]
heatmap(exprs(eset_of_interest_core3_cov4))
heatmap(exprs(eset_of_interest_core3_cov4),
        labCol=gse44532_rma_3_cov4[['cell type:ch1']] ,labRow=NA,
        col       = rev(brewer.pal(10, "RdBu")),
        distfun   = function(x) as.dist(1-cor(t(x))))

#restriction of RMA to only NPC(Hybridization Batch 2) and CD34
gse44532_rma_3_cov5 <- oligo::rma(gse44532_celfiles[,c(1,2,10,11)])

design_core3_cov5 <- model.matrix( ~ 0 + gse44532_rma_3_cov5[['cell type:ch1']])
colnames(design_core3_cov5) <- levels(as.factor(gse44532_rma_3_cov5[['cell type:ch1']]))
kable(design_core3_cov5)

contrast_matrix_core3_cov5 <- makeContrasts(NPC - CD34,
                                            levels=design_core3_cov5)
contrast_matrix_core3_cov5

fit3_cov5 <- lmFit(gse44532_rma_3_cov5,design_core3_cov5)
fit3.2_cov5 <- contrasts.fit(fit3_cov5,contrasts=contrast_matrix_core3_cov5)
fit3.2_cov5 <- eBayes(fit3.2_cov5)
summary(decideTests(fit3.2_cov5,lfc=1))

ps3_cov5 <- topTable(fit3.2_cov5, number=Inf,p.value = 0.05,lfc=1)
ps3_cov5_up <- rownames(ps3_cov5[ps3_cov5$logFC > 0,])
df3_core_cov5 <- AnnotationDbi::select(hugene20sttranscriptcluster.db,ps3_cov5_up,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
dplyr::mutate(df3_core_cov5,GENENAME=stringr::str_trunc(GENENAME,30))

#visualization
#volcano plots CORE
volcanoplot(fit3.2_cov5, coef=1)
#this one had a lot of features pass the cutoff
interesting_genes_core3_cov5 <- topTable(fit3.2_cov5, number=Inf,p.value = 0.05,lfc=2)
volcanoplot(fit3.2_cov5, coef=1, main=sprintf("%d features pass our cutoffs",nrow(interesting_genes_core3_cov5)))
points(interesting_genes_core3_cov5[['logFC']],-log10(interesting_genes_core3_cov5[['P.Value']]),col='red')

eset_of_interest_core3_cov5 <- gse44532_rma_3_cov5[rownames(interesting_genes_core3_cov5),]
heatmap(exprs(eset_of_interest_core3_cov5))
heatmap(exprs(eset_of_interest_core3_cov5),
        labCol=gse44532_rma_3_cov5[['cell type:ch1']] ,labRow=NA,
        col       = rev(brewer.pal(10, "RdBu")),
        distfun   = function(x) as.dist(1-cor(t(x))))


#compare with original
summary(decideTests(fit3.2,lfc=1))
dplyr::mutate(df3_core,GENENAME=stringr::str_trunc(GENENAME,30))
interesting_genes_core3 <- topTable(fit3.2, number=Inf,p.value = 0.05,lfc=2)
volcanoplot(fit3.2, coef=1, main=sprintf("%d features pass our cutoffs",nrow(interesting_genes_core3)))
points(interesting_genes_core3[['logFC']],-log10(interesting_genes_core3[['P.Value']]),col='red')

eset_of_interest_core3 <- gse44532_rma_3[rownames(interesting_genes_core3),]
heatmap(exprs(eset_of_interest_core3))

heatmap(exprs(eset_of_interest_core3),
        labCol=gse44532_rma_3[['cell type:ch1']] ,labRow=NA,
        col       = rev(brewer.pal(10, "RdBu")),
        distfun   = function(x) as.dist(1-cor(t(x))))

