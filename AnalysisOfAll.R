gse44532_rma_all <- oligo::rma(gse44532_celfiles)

oligo::hist(gse44532_rma_all, target="core", main = 'Density Plot of Normalised Data')
oligo::boxplot(gse44532_rma_all, target="core", main = 'Boxplot of Normaliseda')


pd_all <- pData(gse44532_rma_all)
pd_all <- rename(pd_all,cell_type="cell type:ch1")
pd_all$cell_type <- as.factor(pd_all$cell_type)
levels(pd_all$cell_type) <- c("CD34","iNS","iNS.derived.from.iPSC","NPC")
design_all <- model.matrix(~ 0 + pd_all$cell_type)
colnames(design_all) <- levels(pd_all$cell_type)
design_all
contrasts_matrix_all <- makeContrasts(iNS_CD34=iNS - CD34,
                                      iNS.derived.from.iPSC_CD34=iNS.derived.from.iPSC - CD34,
                                      NPC_CD34=NPC - CD34,
                                      levels=design)
kable(contrasts_matrix_all)


fit.all <- lmFit(gse44532_rma_all,design)
fit.all2 <- contrasts.fit(fit.all,contrasts=contrasts_matrix)
fit.all2 <- eBayes(fit.all2)
fisummary(decideTests(fit.all2,lfc=1))

hypo_test <- decideTests(fit.all2,lfc = 1)
summary(hypo_test)
vennDiagram(hypo_test, names = c('iNS/CD34', 'iNS.derived.from.iPSC/CD34', 'NPC/CD34'),
            cex = c(0.8, 0.8, 0.8, 0.8))

#All genes
colnames(ps1.1) <- c('logFC_1', 'AveExpr_1', 't_1', 'P.value_1', 'adj.P.Val_1', 'B_1')
ps1.1$PROBEID <- rownames(ps1.1)
colnames(ps2.1) <- c('logFC_2', 'AveExpr_2', 't_2', 'P.value_2', 'adj.P.Val_2', 'B_2')
ps2.1$PROBEID <- rownames(ps2.1)
colnames(ps3.1) <- c('logFC_3', 'AveExpr_3', 't_3', 'P.value_3', 'adj.P.Val_3', 'B_3')
ps3.1$PROBEID <- rownames(ps3.1)

overlap_1.2 <- inner_join(ps1.1, ps2.1)
overlap_all <- inner_join(overlap_1.2, ps3.1)
rownames(overlap_all) <- overlap_all$PROBEID 

ps_all <- rownames(overlap_all)
df_all <- AnnotationDbi::select(hugene20sttranscriptcluster.db,ps_all,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
dplyr::mutate(df_all,GENENAME=stringr::str_trunc(GENENAME,30))

#Significant genes Only
colnames(ps1) <- c('logFC_1', 'AveExpr_1', 't_1', 'P.value_1', 'adj.P.Val_1', 'B_1')
ps1$PROBEID <- rownames(ps1)
colnames(ps2) <- c('logFC_2', 'AveExpr_2', 't_2', 'P.value_2', 'adj.P.Val_2', 'B_2')
ps2$PROBEID <- rownames(ps2)
colnames(ps3) <- c('logFC_3', 'AveExpr_3', 't_3', 'P.value_3', 'adj.P.Val_3', 'B_3')
ps3$PROBEID <- rownames(ps3)

overlap_1.2_sig <- inner_join(ps1, ps2)
overlap_all_sig <- inner_join(overlap_1.2_sig, ps3)
rownames(overlap_all_sig) <- overlap_all_sig$PROBEID 

ps_all_sig <- rownames(overlap_all_sig)
df_all_sig <- AnnotationDbi::select(hugene20sttranscriptcluster.db,ps_all_sig,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
dplyr::mutate(df_all_sig,GENENAME=stringr::str_trunc(GENENAME,30))

#finding if similar genes that are upregulated in iNS is found in the other two
ps_all_sig_up <- rownames(overlap_all_sig[overlap_all_sig$logFC_1 > 0,])
df_all_sig_up <- AnnotationDbi::select(hugene20sttranscriptcluster.db,ps_all_sig_up,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
dplyr::mutate(df_all_sig_up,GENENAME=stringr::str_trunc(GENENAME,30))

#determining ontology
write.table(df_all_sig_up$ENTREZID, "UpregulatedInALL.txt")

#visualization
eset_of_interest_all <- gse44532_rma_all[rownames(overlap_all_sig),]
heatmap(exprs(eset_of_interest_all))
heatmap(exprs(eset_of_interest_all),
        labCol=gse44532_rma_all[['cell type:ch1']] ,labRow=NA,
        col       = rev(brewer.pal(10, "RdBu")),
        distfun   = function(x) as.dist(1-cor(t(x))))

