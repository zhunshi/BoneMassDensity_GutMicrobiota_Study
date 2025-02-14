library(vegan)
library(cluster)
library(clusterSim)
library(tidyverse)
library(ggpubr)
library(metastat)
library(reshape2)



#-------------------------------------------------------------------------------
# alpha diversity (Shannon)
dat.shannon <- vegan::diversity(profile)                   

# beta diversity (PCoA)
bray_dis <- vegdist(profile, method = 'bray')  
pcoa <- cmdscale(bray_dis,k=2,eig=T)



#-------------------------------------------------------------------------------
# Enterotyping

# JSD dissimilarity
dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}


## Best clusters
kBest <- function(data, dist , method = "kmeans"){
  
  nclusters=NULL
  sil = NULL
  out <- list()
  res <- matrix(NA, 7, ncol(data))
  
  for (k in 2:8) {
    print(k)
    switch (method,
            kmeans = { data.cluster_temp <- kmeans(dist, k)$cluster},
            pam = { data.cluster_temp <- pam(dist, k)$clustering}
    )
    res[k-1,] <- data.cluster_temp
    nclusters[k-1] <- index.G1(t(data) , data.cluster_temp,  d = dist,
                               centrotypes = "medoids")
    sil[k-1] <- mean(silhouette(data.cluster_temp, dist = dist)[,3])
  }
  
  best <- which.max(nclusters)+1
  kCluster <- c(2:8)
  CH_index <- nclusters
  Silhouette <- sil
  cluster <- data.frame(kCluster,  CH_index, Silhouette)
  cluster <- melt(cluster, id = "kCluster")
  colnames(cluster) <- c("kCluster", "Index", "value")
  
  out <- list(res, best, cluster)
  
  return(out)
  
}




#-------------------------------------------------------------------------------
# Partial Spearman's associations
pcc_between <- function (dat1, dat2, dat_adj){
  inte <- intersect(rownames(dat1), rownames(dat2))
  dat1 <- dat1[inte, , drop = F]
  dat2 <- dat2[inte, , drop = F]
  dat_adj <- dat_adj[inte, , drop = F]
  out <- as.data.frame(matrix(NA, ncol(dat1), ncol(dat2) * 
                                3))
  rownames(out) <- colnames(dat1)
  colnames(out) <- paste(rep(colnames(dat2), rep(3, ncol(dat2))), 
                         rep(c("rho", "pvalue", "FDR"), ncol(dat2)))
  for (i in 1:ncol(dat2)) {
    x <- dat2[, i, drop = F]
    for (j in 1:ncol(dat1)) {
      y <- dat1[, j, drop = F]
      a <- cbind(x, y, dat_adj)
      a <- na.omit(a)
      res <- pcor.test(a[, 1, drop = F], a[, 2, drop = F], 
                       dat_adj[rownames(a), , drop = F], method = "spearman")
      out[j, i * 3 - 2] <- res$estimate
      out[j, i * 3 - 1] <- res$p.value
    }
    out[, i * 3] <- p.adjust(out[, i * 3 - 1], method = "BH")
  }
  return(out)
}





                   
#-------------------------------------------------------------------------------








                   
