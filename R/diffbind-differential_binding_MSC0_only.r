library(DiffBind)

samples <- read.csv("design.csv")
samples



samples <- samples[samples$Condition=="day0",]

chaoDBA <- dba(sampleSheet=samples)

plot(chaoDBA)

#Go to QC-consensus-peaks.R first, 
#chaoDBA <- dba.count(chaoDBA, summits=250)

#chaoDBA


###################################################
### code chunk number 12: tamox_aff_corhm
###################################################
plot(chaoDBA)


###################################################
### code chunk number 13: DiffBind.Rnw:243-244
###################################################
#chaoDBA <- dba.contrast(chaoDBA, categories=)
chaoDBA <- dba.contrast(chaoDBA, categories = DBA_TREATMENT)


###################################################
### code chunk number 14: DiffBind.Rnw:254-256
###################################################
chaoDBA <- dba.analyze(chaoDBA)
chaoDBA

plot(chaoDBA, contrast=1)

chaoDBA.DB <- dba.report(chaoDBA, bCounts = TRUE)

###################################################
### code chunk number 17: DiffBind.Rnw:294-295
###################################################
chaoDBA.DB


###################################################
### code chunk number 18: tamox_aff_pca
###################################################
#data(tamoxifen_analysis)
dba.plotPCA(chaoDBA,DBA_TREATMENT,label=DBA_TREATMENT)


###################################################
### code chunk number 19: tamox_sdb_pca
###################################################
dba.plotPCA(chaoDBA, contrast=1,label=DBA_TREATMENT)


###################################################
### code chunk number 20: tamox_sdb_ma
###################################################
dba.plotMA(chaoDBA)


###################################################
### code chunk number 21: tamox_sdb_volcano
###################################################
dba.plotVolcano(chaoDBA)


###################################################
### code chunk number 22: DiffBind.Rnw:399-401
###################################################
sum(chaoDBA.DB$Fold<0)
sum(chaoDBA.DB$Fold>0)


###################################################
### code chunk number 23: tamox_sdb_box
###################################################
pvals <- dba.plotBox(chaoDBA)


###################################################
### code chunk number 24: DiffBind.Rnw:427-428
###################################################
pvals


###################################################
### code chunk number 25: DiffBind.Rnw:440-441
###################################################
corvals <- dba.plotHeatmap(chaoDBA)


###################################################
### code chunk number 26: tamox_sdb_hm
###################################################
corvals <- dba.plotHeatmap(chaoDBA, contrast=1, correlations=FALSE)


###################################################
### code chunk number 27: DiffBind.Rnw:492-495
###################################################
#data(tamoxifen_counts)
chaoDBA_2 <- dba.contrast(chaoDBA,categories=DBA_CONDITION,
                          block=chaoDBA$masks$day0)


###################################################
### code chunk number 28: DiffBind.Rnw:500-502
###################################################
tamoxifen <- dba.analyze(tamoxifen)
tamoxifen


###################################################
### code chunk number 29: tamox_block_ma
###################################################
dba.plotMA(chaoDBA,method=DBA_DESEQ2_BLOCK)


###################################################
### code chunk number 30: DiffBind.Rnw:543-545
###################################################
tamoxifen <- dba.analyze(tamoxifen,method=DBA_ALL_METHODS)
dba.show(tamoxifen,bContrasts=T)[9:12]


###################################################
### code chunk number 31: tamox_block_venn
###################################################
tam.block <- dba.report(tamoxifen,method=DBA_ALL_METHODS_BLOCK,bDB=TRUE,bAll=TRUE)
tam.block
dba.plotVenn(tam.block,1:4,label1="edgeR",label2="DESeq2",
             label3="edgeR Blocked", label4="DESeq2 Blocked")


###################################################
### code chunk number 32: DiffBind.Rnw:573-574
###################################################
data(tamoxifen_peaks)


###################################################
### code chunk number 33: DiffBind.Rnw:595-597
###################################################
olap.rate <- dba.overlap(chaoDBA,mode=DBA_OLAP_RATE)
olap.rate


###################################################
### code chunk number 34: tamox_rate
###################################################
plot(olap.rate,type='b',ylab='# peaks', xlab='Overlap at least this many peaksets') 


###################################################
### code chunk number 35: DiffBind.Rnw:617-618
###################################################
names(chaoDBA$masks)


###################################################
### code chunk number 36: DiffBind.Rnw:623-625
###################################################
dba.overlap(tamoxifen,tamoxifen$masks$MCF7 & tamoxifen$masks$Responsive,
            mode=DBA_OLAP_RATE)


###################################################
### code chunk number 37: tamox_mcf7_venn
###################################################
dba.plotVenn(tamoxifen, tamoxifen$masks$MCF7 & tamoxifen$masks$Responsive)


###################################################
### code chunk number 38: DiffBind.Rnw:647-649
###################################################
tamoxifen_consensus <- dba.peakset(tamoxifen, consensus=c(DBA_TISSUE,DBA_CONDITION),
                                   minOverlap=0.66)


###################################################
### code chunk number 39: DiffBind.Rnw:657-660
###################################################
tamoxifen_consensus <- dba(tamoxifen_consensus, mask=tamoxifen_consensus$masks$Consensus,
                           minOverlap=1)
tamoxifen_consensus


###################################################
### code chunk number 40: DiffBind.Rnw:665-666
###################################################
consensus_peaks <- dba.peakset(tamoxifen_consensus, bRetrieve=TRUE)


###################################################
### code chunk number 41: tamox_lines_venn
###################################################
data(tamoxifen_peaks)
tamoxifen <- dba.peakset(tamoxifen, consensus=DBA_TISSUE, minOverlap=0.66)
dba.plotVenn(tamoxifen, tamoxifen$masks$Consensus)


###################################################
### code chunk number 42: DiffBind.Rnw:690-691
###################################################
data(tamoxifen_peaks)


###################################################
### code chunk number 43: DiffBind.Rnw:696-698
###################################################
dba.overlap(tamoxifen,tamoxifen$masks$Resistant,mode=DBA_OLAP_RATE)
dba.overlap(tamoxifen,tamoxifen$masks$Responsive,mode=DBA_OLAP_RATE)


###################################################
### code chunk number 44: tamox_cons_venn
###################################################
tamoxifen <- dba.peakset(tamoxifen, consensus=DBA_CONDITION, minOverlap=0.33)
dba.plotVenn(tamoxifen,tamoxifen$masks$Consensus)


###################################################
### code chunk number 45: DiffBind.Rnw:727-728
###################################################
tamoxifen.OL <- dba.overlap(tamoxifen, tamoxifen$masks$Consensus)


###################################################
### code chunk number 46: DiffBind.Rnw:733-735
###################################################
tamoxifen.OL$onlyA
tamoxifen.OL$onlyB


###################################################
### code chunk number 47: tamox_compare_venn
###################################################
tamoxifen <- dba.peakset(tamoxifen,tamoxifen$masks$Consensus, 
                         minOverlap=1,sampID="OL Consensus")
tamoxifen <- dba.peakset(tamoxifen,!tamoxifen$masks$Consensus,
                         minOverlap=3,sampID="Consensus_3") 
dba.plotVenn(tamoxifen,14:15)


###################################################
### code chunk number 48: DiffBind.Rnw:762-763
###################################################
data(tamoxifen_analysis)


###################################################
### code chunk number 49: DiffBind.Rnw:768-769
###################################################
tamoxifen.rep <- dba.report(tamoxifen,bCalled=TRUE,th=1)


###################################################
### code chunk number 50: DiffBind.Rnw:774-780
###################################################
onlyResistant <- tamoxifen.rep$Called1>=2 & tamoxifen.rep$Called2<3
sum(onlyResistant )
onlyResponsive <- tamoxifen.rep$Called2>=3 &  tamoxifen.rep$Called1<2
sum(onlyResponsive)
bothGroups <- tamoxifen.rep$Called1>= 2 & tamoxifen.rep$Called2>=3
sum(bothGroups)


###################################################
### code chunk number 51: DiffBind.Rnw:792-799
###################################################
tamoxifen.DB <- dba.report(tamoxifen,bCalled=T)
onlyResistant.DB <- tamoxifen.DB$Called1>=2 & tamoxifen.DB$Called2<3
sum(onlyResistant.DB)
onlyResponsive.DB <- tamoxifen.DB$Called2>=3 & tamoxifen.DB$Called1<2
sum(onlyResponsive.DB)
bothGroups.DB <- tamoxifen.DB$Called1>=2 & tamoxifen.DB$Called2>=3
sum(bothGroups.DB)tamoxifen.DB <- dba.report(tamoxifen)


###################################################
### code chunk number 17: DiffBind.Rnw:294-295
###################################################
tamoxifen.DB


###################################################
### code chunk number 18: tamox_aff_pca
###################################################
data(tamoxifen_analysis)
dba.plotPCA(tamoxifen,DBA_TISSUE,label=DBA_CONDITION)


###################################################
### code chunk number 19: tamox_sdb_pca
###################################################
dba.plotPCA(tamoxifen, contrast=1,label=DBA_TISSUE)


###################################################
### code chunk number 20: tamox_sdb_ma
###################################################
dba.plotMA(tamoxifen)


###################################################
### code chunk number 21: tamox_sdb_volcano
###################################################
dba.plotVolcano(tamoxifen)


###################################################
### code chunk number 22: DiffBind.Rnw:399-401
###################################################
sum(tamoxifen.DB$Fold<0)
sum(tamoxifen.DB$Fold>0)


###################################################
### code chunk number 23: tamox_sdb_box
###################################################
pvals <- dba.plotBox(tamoxifen)


###################################################
### code chunk number 24: DiffBind.Rnw:427-428
###################################################
pvals


###################################################
### code chunk number 25: DiffBind.Rnw:440-441
###################################################
corvals <- dba.plotHeatmap(tamoxifen)


###################################################
### code chunk number 26: tamox_sdb_hm
###################################################
corvals <- dba.plotHeatmap(tamoxifen, contrast=1, correlations=FALSE)


###################################################
### code chunk number 27: DiffBind.Rnw:492-495
###################################################
data(tamoxifen_counts)
tamoxifen <- dba.contrast(tamoxifen,categories=DBA_CONDITION,
                          block=tamoxifen$masks$MCF7)


###################################################
### code chunk number 28: DiffBind.Rnw:500-502
###################################################
tamoxifen <- dba.analyze(tamoxifen)
tamoxifen


###################################################
### code chunk number 29: tamox_block_ma
###################################################
dba.plotMA(tamoxifen,method=DBA_DESEQ2_BLOCK)


###################################################
### code chunk number 30: DiffBind.Rnw:543-545
###################################################
tamoxifen <- dba.analyze(tamoxifen,method=DBA_ALL_METHODS)
dba.show(tamoxifen,bContrasts=T)[9:12]


###################################################
### code chunk number 31: tamox_block_venn
###################################################
tam.block <- dba.report(tamoxifen,method=DBA_ALL_METHODS_BLOCK,bDB=TRUE,bAll=TRUE)
tam.block
dba.plotVenn(tam.block,1:4,label1="edgeR",label2="DESeq2",
             label3="edgeR Blocked", label4="DESeq2 Blocked")


###################################################
### code chunk number 32: DiffBind.Rnw:573-574
###################################################
data(tamoxifen_peaks)


###################################################
### code chunk number 33: DiffBind.Rnw:595-597
###################################################
olap.rate <- dba.overlap(tamoxifen,mode=DBA_OLAP_RATE)
olap.rate


###################################################
### code chunk number 34: tamox_rate
###################################################
plot(olap.rate,type='b',ylab='# peaks', xlab='Overlap at least this many peaksets') 


###################################################
### code chunk number 35: DiffBind.Rnw:617-618
###################################################
names(tamoxifen$masks)


###################################################
### code chunk number 36: DiffBind.Rnw:623-625
###################################################
dba.overlap(tamoxifen,tamoxifen$masks$MCF7 & tamoxifen$masks$Responsive,
            mode=DBA_OLAP_RATE)


###################################################
### code chunk number 37: tamox_mcf7_venn
###################################################
dba.plotVenn(tamoxifen, tamoxifen$masks$MCF7 & tamoxifen$masks$Responsive)


###################################################
### code chunk number 38: DiffBind.Rnw:647-649
###################################################
tamoxifen_consensus <- dba.peakset(tamoxifen, consensus=c(DBA_TISSUE,DBA_CONDITION),
                                   minOverlap=0.66)


###################################################
### code chunk number 39: DiffBind.Rnw:657-660
###################################################
tamoxifen_consensus <- dba(tamoxifen_consensus, mask=tamoxifen_consensus$masks$Consensus,
                           minOverlap=1)
tamoxifen_consensus


###################################################
### code chunk number 40: DiffBind.Rnw:665-666
###################################################
consensus_peaks <- dba.peakset(tamoxifen_consensus, bRetrieve=TRUE)


###################################################
### code chunk number 41: tamox_lines_venn
###################################################
data(tamoxifen_peaks)
tamoxifen <- dba.peakset(tamoxifen, consensus=DBA_TISSUE, minOverlap=0.66)
dba.plotVenn(tamoxifen, tamoxifen$masks$Consensus)


###################################################
### code chunk number 42: DiffBind.Rnw:690-691
###################################################
data(tamoxifen_peaks)


###################################################
### code chunk number 43: DiffBind.Rnw:696-698
###################################################
dba.overlap(tamoxifen,tamoxifen$masks$Resistant,mode=DBA_OLAP_RATE)
dba.overlap(tamoxifen,tamoxifen$masks$Responsive,mode=DBA_OLAP_RATE)


###################################################
### code chunk number 44: tamox_cons_venn
###################################################
tamoxifen <- dba.peakset(tamoxifen, consensus=DBA_CONDITION, minOverlap=0.33)
dba.plotVenn(tamoxifen,tamoxifen$masks$Consensus)


###################################################
### code chunk number 45: DiffBind.Rnw:727-728
###################################################
tamoxifen.OL <- dba.overlap(tamoxifen, tamoxifen$masks$Consensus)


###################################################
### code chunk number 46: DiffBind.Rnw:733-735
###################################################
tamoxifen.OL$onlyA
tamoxifen.OL$onlyB


###################################################
### code chunk number 47: tamox_compare_venn
###################################################
tamoxifen <- dba.peakset(tamoxifen,tamoxifen$masks$Consensus, 
                         minOverlap=1,sampID="OL Consensus")
tamoxifen <- dba.peakset(tamoxifen,!tamoxifen$masks$Consensus,
                         minOverlap=3,sampID="Consensus_3") 
dba.plotVenn(tamoxifen,14:15)


###################################################
### code chunk number 48: DiffBind.Rnw:762-763
###################################################
data(tamoxifen_analysis)


###################################################
### code chunk number 49: DiffBind.Rnw:768-769
###################################################
tamoxifen.rep <- dba.report(tamoxifen,bCalled=TRUE,th=1)


###################################################
### code chunk number 50: DiffBind.Rnw:774-780
###################################################
onlyResistant <- tamoxifen.rep$Called1>=2 & tamoxifen.rep$Called2<3
sum(onlyResistant )
onlyResponsive <- tamoxifen.rep$Called2>=3 &  tamoxifen.rep$Called1<2
sum(onlyResponsive)
bothGroups <- tamoxifen.rep$Called1>= 2 & tamoxifen.rep$Called2>=3
sum(bothGroups)


###################################################
### code chunk number 51: DiffBind.Rnw:792-799
###################################################
tamoxifen.DB <- dba.report(tamoxifen,bCalled=T)
onlyResistant.DB <- tamoxifen.DB$Called1>=2 & tamoxifen.DB$Called2<3
sum(onlyResistant.DB)
onlyResponsive.DB <- tamoxifen.DB$Called2>=3 & tamoxifen.DB$Called1<2
sum(onlyResponsive.DB)
bothGroups.DB <- tamoxifen.DB$Called1>=2 & tamoxifen.DB$Called2>=3
sum(bothGroups.DB)
