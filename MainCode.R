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


# prodominant genera within ETB and ETP
ContributorsEnterotypes <- function(et,pr){
  a <- apply(pr,2,mean)
  a <- as.data.frame(a)
  a <- a[order(a[,1],decreasing = T),,drop=F]
  # only focus on top10 genera
  l <- rownames(a)[1:10]
  a <- cbind(et,pr[rownames(et),l])
  a <- a %>%
    group_by(ET) %>%
    summarise_all(mean) %>%
    gather(taxa,mean,-ET)
  a$taxa <- factor(a$taxa,levels = rev(l))

  # color selection
  col <- colorRampPalette(RColorBrewer::brewer.pal(8,"Set3"))(10)
  col <- rev(col)
  #col <- RColorBrewer::brewer.pal(9,"Set3")
  #col <- rev(col)

  
  ggplot(a,aes(x=ET,y=mean,fill=taxa,alluvium=taxa))+
    ggalluvial::geom_flow(stat = "alluvium",lode.guidance = "frontback", color = "darkgray") + 
    ggalluvial::geom_stratum(stat = "alluvium") +
    scale_fill_manual(values = col)+
    #scale_fill_brewer(palette = "Set3")+
    theme(legend.position = "right")+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    scale_x_discrete(expand = expansion(mult = c(0.15,0.15)))+
    xlab("")+ylab("Mean relative abundance") +
    theme_pubr() +
    theme(
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 13),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13),
      legend.position = "right"
    )
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
# barplot for phenotypes (Fig2e,d)
BarplotPercentage <- function (dat, x, group, yaxies = "p") {
    dat$group1 <- dat[, x]
    dat$group2 <- dat[, group]
    dat <- dat %>% 
    dplyr::select(group1, group2) %>% 
    group_by(group1) %>% 
    count(TargetGroup = group2, .drop = F) %>% 
    mutate(pct = prop.table(n), 
        lab.p = paste0(round(pct * 100, 1), "%"), 
        y.p = pct/2 + c(rev(cumsum(rev(pct))[-length(levels(TargetGroup))]), 0), 
        lab.n = n, 
        y.n = n/2 + c(rev(cumsum(rev(n))[-length(levels(TargetGroup))]), 0))
    y = ifelse(yaxies == "n", "n", "pct")
    ylab = ifelse(yaxies == "n", "Number", "Percentage")
    ggplot(dat, aes_string(x = "group1", y = y, fill = "TargetGroup")) + 
        geom_bar(stat = "identity") + 
        geom_text(aes_string(y = ifelse(ylab == "Number", "y.n", "y.p"), label = "lab.n")) + 
        scale_y_continuous(expand = expansion(mult = c(0,  0.1))) + 
        xlab("") + 
        ylab(ylab) + 
        theme_bw()
}




#-------------------------------------------------------------------------------
# scaltterplot for association between Age and BMD within ETB and ETP (Fig4)
scatterplot <- function (dat, x, y, group = NULL,adj=NULL) {
  dat <- dat[!is.na(dat[, x]), , drop = F]
  if(is.null(adj)){
    s0 <- cor.test(dat[, x], dat[, y], method = "s")
  }else{
    adj <- adj[rownames(dat),,drop=F]
    s0 <- pcor.test(dat[,x],dat[,y],adj,method = "s")
  }
  lst_levels <- levels(as.factor(dat[, group]))
  lab <- c()
  for (l in lst_levels) {
    id <- dat[, group] == l
    if(is.null(adj)){
      a <- cor.test(dat[id, x], dat[id, y], method = "s")
    }else{
      a <- pcor.test(dat[id,x],dat[id,y],adj[id,,drop=F],method = "s")
    }
    lab <- c(lab, paste0(l, " rho=", round(a$estimate, 3), "; p=", formatC(a$p.value, digits = 2)))
  }
  lab <- paste(lab, collapse = "\n")
  p <- ggplot(dat, aes_string(x, y, color = group)) + 
    geom_point(size = 2, alpha = 0.5,aes_string(fill=group),color="white",shape=21) + 
    geom_smooth(method = "lm", se = F, size = 1) + 
    annotate("text", x = -Inf, y = Inf, vjust = 1.2, hjust = 0, label = lab, size = 3) + 
    theme_bw()
  p
}
                   



                   

#-------------------------------------------------------------------------------
# Odds ratio
library(questionr)
dat <- cbind(et[,1,drop=F],dat.adj[rownames(et),],phe[rownames(et),"BMD_group",drop=F])

OR <- function(dat,lab.x,lab.y){
  formula <- as.formula(paste0(lab.y,"~",paste(lab.x,collapse = "+")))
  res <- glm(formula, data = dat, family = binomial)
  
  or <- odds.ratio(res)
  or <- as.data.frame(or)
  or <- or[2,,drop=F]
  
  return(or)
}

MultipleVariableOR <- function(dat,labs.x,labs.y,labs.adj){
  out <- list()
  for(i in labs.x){
    for(j in labs.y){
      name <- paste(j,i,sep = "-")
      out[[name]] <- OR(dat,c(i,labs.adj),j)
    }
  }
  out2 <- do.call(rbind,out)
  out2 <- out2 %>%
    rownames_to_column(var="Factors") %>%
    separate(Factors,c("Disease","Factors"),sep="-") %>%
    mutate(
      Factors = factor(Factors,levels = labs.x),
      Disease = factor(Disease,levels = labs.y)
    )
  return(out2)
}
or <- function(dat1,dat.adj){
  res1 <- MultipleVariableOR(dat1,"Age","y",colnames(dat.adj)[-2])
  res3 <- MultipleVariableOR(dat1[dat1$ET=="ET_B",],"Age","y",colnames(dat.adj)[-2])
  res5 <- MultipleVariableOR(dat1[dat1$ET=="ET_P",],"Age","y",colnames(dat.adj)[-2])

  a <- rbind(res1,res3,res5)
  a$Disease <- "ON/OP vs. NC"
  a$Group <- c("Overall","ET-B",'ET-P')
  a$Group <- factor(a$Group,levels = c("Overall","ET-B",'ET-P'))
  a
}



# forest plot for odds ratio
HR_Forest <- function(dat){
  dat <- dat[-1,]
  dat <- as.data.frame(dat)
  dat$` ` <- paste(rep(" ", 20), collapse = " ")
  dat[,"OR (95% CI)"] <- paste0(round(dat$OR,2)," (",round(dat$`2.5 %`,2)," to ",round(dat$`97.5 %`,2),")")
  dat$p <- formatC(dat$p,digits = 2)
  forestploter::forest(
    data = dat[,c(7,8,9,6),drop=F],
    est=dat$OR,
    lower = dat$`2.5 %`,
    upper = dat$`97.5 %`,
    ci_column = 2,
    ref_line = 1,
    sizes = 0.6,
    xlim = c(0.5,9),
    ticks_at = c(1,2,4,6,8)
  )
}                   

                   
