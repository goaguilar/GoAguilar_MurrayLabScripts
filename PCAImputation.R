readdata = read.csv("probeStructureExpressionMasked9861.csv")
library(pcaMethods)

df<-readdata
df<-df[,-1]
df1<-df[complete.cases(df),]

## Optional to remove only rows with all na
##df <- df1[rowSums(is.na(df1))!=ncol(df1), ]

resPPCA  <- pca(df, method="ppca", center=FALSE, nPcs=5)
resBPCA  <- pca(df, method="bpca", center=FALSE, nPcs=5)
resSVDI  <- pca(df, method="svdImpute", center=FALSE, nPcs=5)
resNipals  <- pca(df, method="nipals", center=FALSE, nPcs=5)
resNLPCA <- pca(df, method="nlpca", center=FALSE, nPcs=5, maxSteps=300)

imputed <- completeObs(resPPCA)

## Generating performance vs fraction missing graph
corlist<-vector("list", 10)
try_values<-matrix(NA, nrow=10, ncol=10)
for(p in 1:10){
  performance_values<-vector("list", 10)
  for(n in 1:10){
    gene_probe_true<-sample(1:19268,1)
    
    x <- as.matrix(df1[gene_probe_true,])
    x1 <- as.vector(x)
    df2 <- na.omit(df)
    samples_to_remove <-sample(1:215,(20*n),replace=F)
    
    for(i in 1:(20*n)){
      df2[gene_probe_true,samples_to_remove[i]] <- NA
    }
    resPPCA <- pca(df2, method="svdImpute", center=FALSE, nPcs=5)
    imputed <- completeObs(resPPCA)
    y1 <- as.vector(imputed[gene_probe_true,])
    performance <- cor(x1,y1)
    performance_values[n] <- performance
  }
  try_values[p,]<-unlist(performance_values)
}
corlist<-matrix(NA, nrow=10, ncol=11)
for (l in 1:10){
  mean_performance<- mean(unlist(nipalsResults[,l]))
  corlist[l,1:10]<-nipalsResults[,l]
  corlist[l,11]<-mean_performance
}

meanlist<-matrix(NA, nrow=10, ncol=4)
#1st is ppca,2nd is bpca,3rd is svdImpute,4th is nipals
meanlist[,4]<-corlist[,11]

#plots for one pca method and its individual methods
#setup
cl<-rainbow(11)
for(z in 1:10){
  cl[z]<-"#FFBF00"
}
cl[11]<-"#000000"
fraction <- seq(20, 200, 20)
ylabeltick<-seq(1,10,1)
fraction <- (fraction/215) *100
fraction <-round(fraction)

matplot(corlist, type = "b",pch=1, col = 1:4,xaxt="n",xlab="Fraction Missing", ylab="Performance", main="Method Comparison")
axis(side=1,at=ylabeltick,labels=fraction)

#plots for comparing the mean plots of each method
matplot(meanlist, type = "b",pch=1, col = 1:4,xaxt="n",xlab="Fraction Missing", ylab="Performance", main="Method Comparison")
axis(side=1,at=ylabeltick,labels=fraction)
legend("bottomleft",legend =c("PPCA","BPCA","svdImpute","Nipals"),col=1:4,pch=1)
