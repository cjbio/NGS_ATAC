source("http://bioconductor.org./biocLite.R")
biocLite("ChIPseeker")
biocLite("clusterProfiler")
biocLite("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")

library("ChIPseeker")
library("clusterProfiler")
library("GenomicRanges")
library("org.Hs.eg.db")
x <- as.data.frame(chaoDBA.DB)
#seqinfo(chaoDBA.DB)
#x[-x$seqnames=="GL000205.2",]



y <- x[x$seqnames=="1"|x$seqnames=="2"|x$seqnames=="3"|x$seqnames=="4"|x$seqnames=="5"|x$seqnames=="6"|x$seqnames=="7"|x$seqnames=="8"
  |x$seqnames=="9"|x$seqnames=="10"|x$seqnames=="11"|x$seqnames=="12"|x$seqnames=="13"|x$seqnames=="14"|x$seqnames=="15"
  |x$seqnames=="4"|x$seqnames=="4"|x$seqnames=="4"|x$seqnames=="4"|x$seqnames=="4"|x$seqnames=="4"|x$seqnames=="4"
  |x$seqnames=="16"|x$seqnames=="17"|x$seqnames=="18"|x$seqnames=="19"|x$seqnames=="20"|x$seqnames=="21"|x$seqnames=="22"
  |x$seqnames=="X"|x$seqnames=="Y",]

y$seqnames <- sub("^", "chr", y$seqnames)
#View(y)
gr <- makeGRangesFromDataFrame(y, ignore.strand = T, seqnames.field = "seqnames",
                               start.field = "start", end.field = "end")
#View(gr)
anno_chaoDBA.DB <- annotatePeak(gr,tssRegion = c(-3000,3000), TxDb =TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb = "org.Hs.eg.db", level = "gene")
#anno_chaoDBA.DB <- annotatePeak(gr)
?annotatePeak()
#View(anno_chaoDBA.DB)
anno_go <- cbind(as.data.frame(anno_chaoDBA.DB),y$Fold)
#View(anno_go)

anno_up <- as.data.frame(anno_go)[anno_go$`y$Fold`>0,]
anno_down <- as.data.frame(anno_go)[anno_go$`y$Fold`<0,]
go_up <- enrichGO(as.data.frame(anno_go)[anno_go$`y$Fold`>0,]$geneId, OrgDb = "org.Hs.eg.db", ont ="BP")
go_down <- enrichGO(as.data.frame(anno_go)[anno_go$`y$Fold`<0,]$geneId, OrgDb = "org.Hs.eg.db", ont ="BP")
write.csv(anno_go, "anno_consensus_peak.csv")
write.csv(go_up, "enrichGO_OpenWithAge_sig.csv")
write.csv(go_down, "enrichGO_CloseWithAge_sig.csv")
View(go_up)
View(go_down)

up_skel_genes <- unlist(strsplit(go_up[1,]$geneID, split="/"))
down_skel_genes <- unlist(strsplit(go_down[1,]$geneID, split="/"))



library(annotate)
unlist(lookUp(up_skel_genes, "org.Hs.eg.db", "SYMBOL"))
down_skel_genes <- unlist(lookUp(down_skel_genes, "org.Hs.eg.db", "SYMBOL"))
write.csv(down_skel_genes, "enrichGO_SigDownSkeletalGenes.csv")
View(anno_go)
##########################################################

anno_up_2 <- as.data.frame(anno_up)[anno_up$`y$Fold`>1,]
anno_down_2 <- anno_down[anno_down$`y$Fold`< -1,]

genes_up_2 <- anno_up$SYMBOL
genes_down_2 <- anno_down_2$SYMBOL
genes_up_2 <- genes_up_2[order(genes_up_2)]
genes_down_2 <- genes_down_2[order(genes_down_2)]
genes_up_2 <- unique(genes_up_2)
genes_down_2 <- unique(genes_down_2)

write.csv(anno_up_2,"genes_up_log2fold>1_MSC0_only.csv")
write.csv(anno_down_2,"genes_down_log2fold>1_MSC0_only.csv")

go_up_2 <- enrichGO(anno_up_2$geneId, OrgDb = "org.Hs.eg.db", ont ="BP")
go_down_2 <- enrichGO(anno_down_2$geneId, OrgDb = "org.Hs.eg.db", ont ="BP")


View(go_up_2)
View(go_down_2)
up_skel_genes_2 <- unlist(strsplit(go_up_2[1,]$geneID, split="/"))
down_skel_genes_2 <- unlist(strsplit(go_down_2[1,]$geneID, split="/"))

#library(annotate)
unlist(lookUp(up_skel_genes_2, "org.Hs.eg.db", "SYMBOL"))
unlist(lookUp(down_skel_genes_2, "org.Hs.eg.db", "SYMBOL"))

save.image("complete_analysis_Diffbind_chIPseeker_GO.RData")
load("complete_analysis_Diffbind_chIPseeker_GO.RData")

