#rm(list = ls())
olap.rate<- dba.overlap(chaoDBA, mode = DBA_OLAP_RATE)
olap.rate

plot(olap.rate,type='b',ylab='# peaks', xlab='Overlap at least this many peaksets')

names((chaoDBA$masks))

dba.overlap(chaoDBA, chaoDBA$masks$day0 & chaoDBA$masks$old, mode = DBA_OLAP_RATE)

pdf("peak overlap between biological replicates.pdf")
dba.plotVenn(chaoDBA, chaoDBA$masks$day0 & chaoDBA$masks$old)
dba.plotVenn(chaoDBA, chaoDBA$masks$day0 & chaoDBA$masks$young)
dba.plotVenn(chaoDBA, chaoDBA$masks$day3 & chaoDBA$masks$old)
dba.plotVenn(chaoDBA, chaoDBA$masks$day3 & chaoDBA$masks$young)
dba.plotVenn(chaoDBA, chaoDBA$masks$day7 & chaoDBA$masks$old)
dba.plotVenn(chaoDBA, chaoDBA$masks$day7 & chaoDBA$masks$young)
dba.plotVenn(chaoDBA, chaoDBA$masks$day21 & chaoDBA$masks$old)
dba.plotVenn(chaoDBA, chaoDBA$masks$day21 & chaoDBA$masks$young)
dev.off()

chaoDBA_consensus <- dba.peakset(chaoDBA, consensus = c(DBA_CONDITION, DBA_TREATMENT), minOverlap = 2)

chaoDBA_consensus
chaoDBA_consensus <- dba(chaoDBA_consensus, mask = chaoDBA_consensus$masks$Consensus, minOverlap = 1)
consensus_peaks <- dba.peakset(chaoDBA_consensus, bRetrieve = TRUE) 


chaoDBA_co <- dba.count(chaoDBA, peaks = consensus_peaks)
chaoDBA <- chaoDBA_co

#check the sum of counted reads for all samples, change [[1]] to [[6]]
sum(chaoDBA$peaks[[1]][,5])
chaoDBA$samples
#similar total read counts