#For Core-level output, NPC vs CD34 [for diff sample source]
#restriction of RMA to only NPC(Sample 5) and CD34
gse44532_rma_3_cov1 <- oligo::rma(gse44532_celfiles[,c(1,2,8,9)])

design_core3_cov1 <- model.matrix( ~ 0 + gse44532_rma_3_cov1[['cell type:ch1']])
colnames(design_core3_cov1) <- levels(as.factor(gse44532_rma_3_cov1[['cell type:ch1']]))
kable(design_core3_cov1)

contrast_matrix_core3_cov1 <- makeContrasts(NPC - CD34,
                                            levels=design_core3_cov1)
contrast_matrix_core3_cov1

fit3_cov1 <- lmFit(gse44532_rma_3_cov1,design_core3_cov1)
fit3.2_cov1 <- contrasts.fit(fit3_cov1,contrasts=contrast_matrix_core3_cov1)
fit3.2_cov1 <- eBayes(fit3.2_cov1)
summary(decideTests(fit3.2_cov1,lfc=1))

ps3_cov1 <- topTable(fit3.2_cov1, number=Inf,p.value = 0.05,lfc=1)
ps3_cov1_up <- rownames(ps3_cov1[ps3_cov1$logFC > 0,])
df3_core_cov1 <- AnnotationDbi::select(hugene20sttranscriptcluster.db,ps3_cov1_up,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
dplyr::mutate(df3_core_cov1,GENENAME=stringr::str_trunc(GENENAME,30))

#visualization
#volcano plots CORE
volcanoplot(fit3.2_cov1, coef=1)
#this one had a lot of features pass the cutoff
interesting_genes_core3_cov1 <- topTable(fit3.2_cov1, number=Inf,p.value = 0.05,lfc=2)
volcanoplot(fit3.2_cov1, coef=1, main=sprintf("%d features pass our cutoffs",nrow(interesting_genes_core3_cov1)))
points(interesting_genes_core3_cov1[['logFC']],-log10(interesting_genes_core3_cov1[['P.Value']]),col='red')

eset_of_interest_core3_cov1 <- gse44532_rma_3_cov1[rownames(interesting_genes_core3_cov1),]
heatmap(exprs(eset_of_interest_core3_cov1))
heatmap(exprs(eset_of_interest_core3_cov1),
        labCol=gse44532_rma_3_cov1[['cell type:ch1']] ,labRow=NA,
        col       = rev(brewer.pal(10, "RdBu")),
        distfun   = function(x) as.dist(1-cor(t(x))))

#restriction of RMA to only NPC(Sample 6 only) and CD34
gse44532_rma_3_cov2 <- oligo::rma(gse44532_celfiles[,c(1,2,10)])

design_core3_cov2 <- model.matrix( ~ 0 + gse44532_rma_3_cov2[['cell type:ch1']])
colnames(design_core3_cov2) <- levels(as.factor(gse44532_rma_3_cov2[['cell type:ch1']]))
kable(design_core3_cov2)

contrast_matrix_core3_cov2 <- makeContrasts(NPC - CD34,
                                            levels=design_core3_cov2)
contrast_matrix_core3_cov2

fit3_cov2 <- lmFit(gse44532_rma_3_cov2,design_core3_cov2)
fit3.2_cov2 <- contrasts.fit(fit3_cov2,contrasts=contrast_matrix_core3_cov2)
fit3.2_cov2 <- eBayes(fit3.2_cov2)
summary(decideTests(fit3.2_cov2,lfc=1))

ps3_cov2 <- topTable(fit3.2_cov2, number=Inf,p.value = 0.05,lfc=1)
ps3_cov2_up <- rownames(ps3_cov2[ps3_cov2$logFC > 0,])
df3_core_cov2 <- AnnotationDbi::select(hugene20sttranscriptcluster.db,ps3_cov2_up,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
dplyr::mutate(df3_core_cov2,GENENAME=stringr::str_trunc(GENENAME,30))

#visualization
#volcano plots CORE
volcanoplot(fit3.2_cov2, coef=1)
#this one had a lot of features pass the cutoff
interesting_genes_core3_cov2 <- topTable(fit3.2_cov2, number=Inf,p.value = 0.05,lfc=2)
volcanoplot(fit3.2_cov2, coef=1, main=sprintf("%d features pass our cutoffs",nrow(interesting_genes_core3_cov2)))
points(interesting_genes_core3_cov2[['logFC']],-log10(interesting_genes_core3_cov2[['P.Value']]),col='red')

eset_of_interest_core3_cov2 <- gse44532_rma_3_cov2[rownames(interesting_genes_core3_cov2),]
heatmap(exprs(eset_of_interest_core3_cov2))
heatmap(exprs(eset_of_interest_core3_cov2),
        labCol=gse44532_rma_3_cov2[['cell type:ch1']] ,labRow=NA,
        col       = rev(brewer.pal(10, "RdBu")),
        distfun   = function(x) as.dist(1-cor(t(x))))


#restriction of RMA to only NPC(Sample 7 only) and CD34
gse44532_rma_3_cov3 <- oligo::rma(gse44532_celfiles[,c(1,2,11)])

design_core3_cov3 <- model.matrix( ~ 0 + gse44532_rma_3_cov3[['cell type:ch1']])
colnames(design_core3_cov3) <- levels(as.factor(gse44532_rma_3_cov3[['cell type:ch1']]))
kable(design_core3_cov3)

contrast_matrix_core3_cov3 <- makeContrasts(NPC - CD34,
                                            levels=design_core3_cov3)
contrast_matrix_core3_cov3

fit3_cov3 <- lmFit(gse44532_rma_3_cov3,design_core3_cov3)
fit3.2_cov3 <- contrasts.fit(fit3_cov3,contrasts=contrast_matrix_core3_cov3)
fit3.2_cov3 <- eBayes(fit3.2_cov3)
summary(decideTests(fit3.2_cov3,lfc=1))

ps3_cov3 <- topTable(fit3.2_cov3, number=Inf,p.value = 0.05,lfc=1)
ps3_cov3_up <- rownames(ps3_cov3[ps3_cov3$logFC > 0,])
df3_core_cov3 <- AnnotationDbi::select(hugene20sttranscriptcluster.db,ps3_cov3_up,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
dplyr::mutate(df3_core_cov3,GENENAME=stringr::str_trunc(GENENAME,30))

#visualization
#volcano plots CORE
volcanoplot(fit3.2_cov3, coef=1)
#this one had a lot of features pass the cutoff
interesting_genes_core3_cov3 <- topTable(fit3.2_cov3, number=Inf,p.value = 0.05,lfc=2)
volcanoplot(fit3.2_cov3, coef=1, main=sprintf("%d features pass our cutoffs",nrow(interesting_genes_core3_cov3)))
points(interesting_genes_core3_cov3[['logFC']],-log10(interesting_genes_core3_cov3[['P.Value']]),col='red')

eset_of_interest_core3_cov3 <- gse44532_rma_3_cov3[rownames(interesting_genes_core3_cov3),]
heatmap(exprs(eset_of_interest_core3_cov3))
heatmap(exprs(eset_of_interest_core3_cov3),
        labCol=gse44532_rma_3_cov3[['cell type:ch1']] ,labRow=NA,
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
