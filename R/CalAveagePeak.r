
sigPeaks <- read.csv("sigPeaks.csv", header = TRUE)
sigPeaks <- as.data.frame(sigPeaks)


sigPeakAvg <- aggregate(sigPeaks$y.Fold, list(sigPeaks$SYMBOL), mean)
View(sigPeakAvg)
length(rownames(sigPeakAvg))
length(unique(sigPeakAvg$Group.1))

write.csv(sigPeakAvg, "sigPeakAvg.csv")

# name base on the age group and timepoint/condition
tp_age <- paste(as.character(coldata$condition), as.character(coldata$agegroup))
colnames(named_nc_sig) <- tp_age[-c(67,69,70,71)]
test1 <- named_nc_sig


# calculate the mean of the columns with the same column name
test2 <- as.data.frame(
  sapply(unique(names(test1)), 
         function(col) rowMeans(test1[names(test1) == col]) # calculate row means
  )
)

# check if the cal is correct
asd <- test1[1,]
asd[colnames(asd)=="0 young"]
359.601 + 277.0566 + 433.0526
3.565700e+02*3


colnames(test2)
old_sig <- test2[,c("0 old","0.5 old","1 old","2 old","4 old","6 old","8 old","12 old","16 old","24 old","48 old","72 old","96 old","120 old","168 old","336 old","504 old")]
young_sig <- test2[, c("0 young","0.5 young","1 young","2 young","4 young","6 young","8 young","12 young","16 young","24 young","48 young","72 young","96 young","120 young","168 young","336 young","504 young")]
#"0 old","0.5 old","1 old","2 old","4 old","6 old","8 old","12 old","16 old","24 old","48 old","72 old","96 old","120 old","168 old","336 old","504 old"
#"0 young","0.5 young","1 young","2 young","4 young","6 young","8 young","12 young","16 young","24 young","48 young","72 young","96 young","120 young","168 young","336 young","504 young"

