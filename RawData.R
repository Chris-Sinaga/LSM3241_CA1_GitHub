gse44532 <- getGEO('GSE44532')
gse44532 <- gse44532[[1]]

gse44532_SuppFiles <- getGEOSuppFiles('GSE44532')
untar("GSE44532/GSE44532_RAW.tar", exdir = "GSE44532")

gse44532_celfiles <- read.celfiles(list.celfiles('GSE44532',
                                                 full.names=TRUE,
                                                 listGzipped=TRUE))

oligo::hist(gse44532_celfiles, target="core", main = 'Density Plot of Raw Data')
oligo::boxplot(gse44532_celfiles, target="core", main = 'Boxplot of Raw Data')

pd <- pData(gse44532)
pd['cel_file'] <- str_split(pd$supplementary_file,"/") %>% map_chr(tail,1)

gse44532_celfiles <- read.celfiles(paste0('GSE44532/',pd$cel_file),
                                   phenoData=phenoData(gse44532))

#to identify which columns contain the information we need
pData(gse44532_celfiles)

#initial data visualization
image(gse44532_celfiles[,1])

#information of interest are in the following columns
pData(gse44532)[,c("geo_accession","cell type:ch1","source sample:ch1",
                   "hybridization batch:ch1")]

#exporting table (u can change the directory for your own if you wanna run this again)
write.table(pData(gse44532)[,c("geo_accession","cell type:ch1","source sample:ch1",
                               "hybridization batch:ch1")], 
            "C:/Users/Christopher Sinaga/Desktop/LSM3241_CA1_FinalDraft/ges44532.txt", 
            sep="\t")

#moving forward, 
#1 <- analysis of iNS vs CD34
#2 <- analysis of iNS derived from iPSC vs CD34
#3 <- analysis of NPC vs CD34