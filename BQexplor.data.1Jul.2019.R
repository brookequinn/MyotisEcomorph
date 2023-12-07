## Script to create boxplots of exploratory data for Myotis project ##

# Set working directory
setwd("C:/Users/brook/Documents/R/amnh_home")


# Read .txt file of data
mydata<-read.table(file="SpecimensBQ_6.26.2019.txt", header = TRUE, sep = "\t", skipNul = TRUE, check.names = F)

##### Create a basic set of boxplots for each trait #####
pdf(file="boxplot_6.19.2019.pdf")
par(mfrow=c(4,5))
traitCols=12:25
for(trait in traitCols){
  boxplot(mydata[,trait], main=colnames(mydata)[trait], col = "tomato")
}
dev.off()

##### Create set of boxplots by ecomorph #####
pdf(file="boxplot$ecomorph_6.19.2019.pdf")
par(mfrow=c(4,5))

for(trait in traitCols){
  boxplot(mydata[,trait]~mydata$Ecomorph, main=colnames(mydata)[trait], col = c("darkgoldenrod1", "deepskyblue1", "firebrick1"))
} 

dev.off()
  
##### Create a set of boxplots by biogeographic region #####
pdf(file="boxplot$biogeo_6.19.2019.pdf", width=10, height=8)
par(mfrow=c(4,5))

for(trait in traitCols){
  boxplot(mydata[,trait]~mydata$BioReg, main=colnames(mydata)[trait], mar=c(3,2,2,0.4), xlab="Biogeographic Region", ylab="Length (mm)", col = c(mypalette), las=3, labels=c("E", "EP", "Nr", "Nt", "Oc", "Or", "WP"))
} 
# Colors () correspond to ("E" EP" "Nr" "Nt" "Oc" "Or" "WP")
dev.off()

###########################################################################
##### Create tables of means and standard deviations for each species #####

#Aggregate data by species tree label and find means
traitCols=12:25
for(trait in traitCols){
  newmeans<-aggregate(mydata[,traitCols], list(mydata$SpeciesTreeLabel),FUN=mean)
}
	colnames(newmeans)[1]<-"SpeciesTreeLabel"
write.table(newmeans,file="SpeciesMeanTraitValues.txt")

#Aggregate data by species tree label and find standard deviation. If number of species is less than 2, results in NA
for(trait in traitCols){
  newsds<-aggregate(mydata[,traitCols], list(mydata$SpeciesTreeLabel),FUN=sd)
}
	colnames(newsds)[1]<-"SpeciesTreeLabel"
write.table(newsds,file="SpeciesSDTraitValues.txt")

#Aggregate data by species tree label and find min values
for(trait in traitCols){
  newmins<-aggregate(mydata[,traitCols], list(mydata$SpeciesTreeLabel),FUN=min)
}
# colnames(newmins)[1]<-c()
write.table(newmins,file="SpeciesMinTraitValues.txt", sep="\t")

#Aggregate data by species tree label and find max values
for(trait in traitCols){
  newmaxs<-aggregate(mydata[,traitCols], list(mydata$SpeciesTreeLabel),FUN=max)
}
colnames(newmaxs)[1]<-"SpeciesTreeLabel"
write.table(newmaxs,file="SpeciesMaxTraitValues.txt",sep="\t")


##### Create histograms for each trait #####
pdf(file = "myhist.pdf", width = 10, height = 8)
par(mfrow=c(4,5))
    for(trait in traitCols){
    hist = hist(mydata[,trait],freq=TRUE,main=colnames(mydata)[trait],xlab="Trait (mm)",col="gold",plot=TRUE)
}
dev.off()

#######################################################
##### Shapiro test for normality of all variables #####
traitCols=12:25
resultsNormTest<-c()
for(trait in traitCols){
  normTest<-shapiro.test(mydata[,trait])
  resultsNormTest<-rbind(resultsNormTest, c("w"=round(normTest$statistic,9), "p.value"=round(normTest$p.value, 9)))
}
rownames(resultsNormTest)<-colnames(mydata[,traitCols])
write.table(resultsNormTest, file="resultsNormTest.3ecomorphs.txt", quote=FALSE)
  
#Testing normality for all traits within Leuconoe ecomorph
resultsNormTest.L<-c()
L.eco = dplyr::filter(mydata,Ecomorph == "L")
for(trait in traitCols){
  normTest.L<-shapiro.test(L.eco[,trait])
  resultsNormTest.L<-rbind(resultsNormTest.L, c("w"=round(normTest.L$statistic, 9), "p-value"=round(normTest.L$p.value,9)))
}
rownames(resultsNormTest.L)=colnames(mydata[,traitCols])
write.table(resultsNormTest.L,file="resultsNormTest.Lecomorph.AM.txt", quote=FALSE)

#Testing normality for all traits within Myotis ecomorph
resultsNormTest.M<-c()
M.eco = dplyr::filter(mydata, Ecomorph == "M")
for(trait in traitCols){
  normTest.M<-shapiro.test(M.eco[,trait])
  resultsNormTest.M<-rbind(resultsNormTest.M, c("w"=round(normTest.M$statistic, 9), "p-value"=round(normTest.M$p.value,9)))
}
rownames(resultsNormTest.M)=colnames(mydata[,traitCols])
write.table(resultsNormTest.M,file="resultsNormTest.Mecomorph.AM.txt", quote=FALSE)

#Testing normality for all traits within Selysius ecomorph
resultsNormTest.S=c()
S.eco = dplyr::filter(mydata, Ecomorph == "S")
for(trait in traitCols){
  normTest.S<-shapiro.test(S.eco[,trait])
  resultsNormTest.S=rbind(resultsNormTest.S,c("w"=round(normTest.S$statistic,9),"p-value"=round(normTest.S$p.value,9)))
}
rownames(resultsNormTest.S)=colnames(mydata[,traitCols])
write.table(resultsNormTest.S, file="resultsNormTest.Secomorph.AM.txt", quote=FALSE)


##########################################
##### T tests between ecomorph pairs #####
#Leuconoe and myotis
traitCols=12:25
t.LM.summary<-c()
for(trait in traitCols){
  ttest.LM = t.test(L.eco[,trait], M.eco[,trait])
  results.ttest =c(ttest.LM$statistic,ttest.LM$p.value,ttest.LM$conf.int, ttest.LM$parameter,ttest.LM$estimate)
    t.LM.summary<-rbind(t.LM.summary, results.ttest)
  }
  colnames(t.LM.summary)=c("t stat","p value","lower conf int","upper conf int", "df","mean of x", "mean of y")
  rownames(t.LM.summary)<-colnames(L.eco[,12:25])
write.table(t.LM.summary,file="tTestLMecomorphs.AM.txt",quote=FALSE)



#Leuconoe and selysius
traitCols=12:25
t.LS.summary<-c()
for(trait in traitCols){
  ttest.LS = t.test(L.eco[,trait], S.eco[,trait])
  results.ttest =c(ttest.LS$statistic,ttest.LS$p.value,ttest.LS$conf.int, ttest.LS$parameter,ttest.LS$estimate)
    t.LS.summary<-rbind(t.LS.summary, results.ttest)
  }
  colnames(t.LS.summary)=c("t stat","p value","lower conf int","upper conf int", "df","mean of x", "mean of y")
  rownames(t.LS.summary)<-colnames(L.eco[,12:25])
write.table(t.LS.summary,file="tTestLSecomorphs.AM.txt",quote=FALSE)


#Selysius and myotis
traitCols=12:25
t.SM.summary<-c()
for(trait in traitCols){
  ttest.SM = t.test(S.eco[,trait], M.eco[,trait])
  results.ttest =c(ttest.SM$statistic,ttest.SM$p.value,ttest.SM$conf.int, ttest.SM$parameter,ttest.SM$estimate)
    t.SM.summary<-rbind(t.SM.summary, results.ttest)
}
  colnames(t.SM.summary)=c("t stat","p value","lower conf int","upper conf int", "df","mean of x", "mean of y")
  rownames(t.SM.summary)<-colnames(L.eco[,12:25])
write.table(t.SM.summary,file="tTestSMecomorphs.AM.txt",quote=FALSE)



############################
### T test between sexes ###

male = dplyr::filter(mydata,Sex == "M")
female = dplyr::filter(mydata,Sex == "F")
t.sex.summary<-c()
for(trait in traitCols){
  ttest.sexes = t.test(female[,trait], male[,trait])
  results.ttest = c(ttest.sexes$statistic,ttest.sexes$p.value,ttest.sexes$conf.int, ttest.sexes$parameter,ttest.sexes$estimate)
     t.sex.summary<-rbind(t.sex.summary, results.ttest)
  #row.names(results.ttest)=c("t stat","p value","lower conf int","upper conf int")
}
	colnames(t.sex.summary)=c("t stat","p value","lower conf int","upper conf int", "df","mean of x", "mean of y")
  	rownames(t.sex.summary)<-colnames(L.eco[,12:25])
write.table(results.ttest,file="tTest.sexes.AM.txt",quote=FALSE)



##### Supplementary table of ecomorph by biogeography #####
datadist = with(mydata, table(BioReg,Ecomorph))
write.table(datadist,file="Ecomorph.Bioreg.SuppTable.txt",quote=FALSE, sep="\t")

##### ANOVA test #####
# ANOVA for ecomorphs
traitCols=12:25
eco=mydata$Ecomorph
myanova.eco<-c()
for(trait in traitCols){
  myanova = aov(mydata[,trait]~eco, data=mydata)
  anova=summary(myanova)[[1]][1,]
  anova2=summary(myanova)[[1]][2,]
  anovafull=rbind(anova, anova2)
  myanova.eco=rbind(myanova.eco,anovafull)
}
rows =c(rep(colnames(mydata[,traitCols]),1,each=2))
rows2 = make.names(rows, unique=TRUE)
rownames(myanova.eco)<-rows2
colnames(myanova.eco)=c("Df","Sum sq","Mean sq","F","p")
write.table(myanova.eco,file="anovatest.txt", sep="\t")

# ANOVA for bioregions
bioreg=mydata$BioReg
myanova.bio<-c()
for(x in traitCols){
  anova.output = aov(mydata[,x]~bioreg, data=mydata)
  anova=summary(myanova)[[1]][1,]
  myanova.bio=rbind(myanova.bio,anova)
}
colnames(myanova.bio)=c("Df","Sum sq","Mean sq","F","p")
rownames(myanova.bio)<-colnames(mydata[,traitCols])
write.table(myanova.bio,file="ANOVAoutput.bio.txt")


##### Random forest to test variable importance #####
#Load library. version 4.6.14
library(randomForest)

#Rank variables in order of importance. Makes 1000 trees, permutes OOB data 100 times per tree to assess importance, and importance = true. 
RandFor.tmp <- randomForest(mydata$Ecomorph~.,data=mydata[,setdiff(9:24, seq(10))],importance=TRUE, ntree=1000, nperm=100)
#Save as data frame
#Order variables by importance and list names in order
imp <- data.frame(importance(RF.tmp))
INC <- imp[order(imp$MeanDecreaseAccuracy,decreasing = T),]

#Save output as pdf file graph
#List variable names by importance
impvar <- rownames(imp)[order(imp[, 1], decreasing=TRUE)]
pdf(file="RF.varImp.pdf", width=10, height=5)
ImpPlot<-varImpPlot(RF.tmp, sort=TRUE, main="Variables in Order of Decreasing Importance")
dev.off()
  




