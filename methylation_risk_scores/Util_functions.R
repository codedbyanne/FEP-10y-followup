AUC <- function (actuals, predictedScores){
  fitted <- data.frame (Actuals=actuals, PredictedScores=predictedScores)
  colnames(fitted) <- c('Actuals','PredictedScores')
  ones <- fitted[fitted$Actuals==1, ] # Subset ones
  zeros <- fitted[fitted$Actuals==0, ] # Subsetzeros
  totalPairs <- nrow (ones) * nrow (zeros) # calculate total number of pairs to check
  conc <- sum (c(vapply(ones$PredictedScores, function(x) {((x > zeros$PredictedScores))}, FUN.VALUE=logical(nrow(zeros)))), na.rm=T)
  disc <- sum(c(vapply(ones$PredictedScores, function(x) {((x < zeros$PredictedScores))}, FUN.VALUE = logical(nrow(zeros)))), na.rm = T)
  concordance <- conc/totalPairs
  discordance <- disc/totalPairs
  tiesPercent <- (1-concordance-discordance)
  AUC = concordance + 0.5*tiesPercent
  Gini = 2*AUC - 1
  return(list("Concordance"=concordance, "Discordance"=discordance,
              "Tied"=tiesPercent, "AUC"=AUC, "Gini or Somers D"=Gini))
}


prs.pc <- function(dat,x){
  xo <- scale(as.matrix(dat[,-1]))  ## scale cols of matrix of only PRSs (remove ID)
  g <- prcomp(xo)   ## perform principal components
  pca.r2 <- g$sdev^2/sum(g$sdev^2)    ## calculate variance explained by each PC
  pc1.loadings <- g$rotation[,1];     ## loadings for PC1
  pc2.loadings <- g$rotation[,2]      ## loadings for PC2
  ## flip direction of PCs to keep direction of association
  ## (sign of loadings for PC1 is arbitrary so we want to keep same direction)
  if (mean(pc1.loadings>0)==0){     
    pc1.loadings <- pc1.loadings*(-1) 
    pc2.loadings <- pc2.loadings*(-1)
  }
  ## calculate PRS-PCA (outputs PC1 and PC2 even though PC1 sufficient)
  pc1 <- xo %*% pc1.loadings
  pc2 <- xo %*% pc2.loadings
  dat[,paste0(x,".prs.pc")] <- scale(pc1)   ## rescales PRS-PCA1 
  dat[,paste0(x,".prs.pc2")] <- scale(pc2)  ## rescales PRS-PCA2
  return(list(data=dat,r2=pca.r2,loadings=pc1.loadings))
}

