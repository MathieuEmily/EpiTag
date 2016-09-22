###
## getProba.Paire: estimate the marginal joint probabilities for the 2 pairs in a 4x4 joint probability table
###

getProba.Paire <- function(proba.jointe){
	p.AB <- matrix(NA,ncol=2,nrow=2)
	p.AB[1,1] <- sum(proba.jointe[,1])
	p.AB[1,2] <- sum(proba.jointe[,2])
	p.AB[2,1] <- sum(proba.jointe[,3])
	p.AB[2,2] <- sum(proba.jointe[,4])
	p.CD <- matrix(NA,ncol=2,nrow=2)
	p.CD[1,1] <- sum(proba.jointe[1,])
	p.CD[1,2] <- sum(proba.jointe[2,])
	p.CD[2,1] <- sum(proba.jointe[3,])
	p.CD[2,2] <- sum(proba.jointe[4,])
	return(list(p.AB=p.AB,p.CD=p.CD))
}

###
## getI.allelic: estimate the Information between two pairs of SNPs in allelic format
###

getI.allelic <- function(proba.jointe=NULL,SNP1=NA,SNP2=NA,tSNP1=NA,tSNP2=NA){
	pj.test <- proba.jointe
	if (is.null(pj.test)){
		n.obs <- length(SNP1)
		SNP1 <- factor(SNP1,levels=c(0,1))
		SNP2 <- factor(SNP2,levels=c(0,1))
		tSNP1 <- factor(tSNP1,levels=c(0,1))
		tSNP2 <- factor(tSNP2,levels=c(0,1))
		
		proba.jointe <- table(SNP1,SNP2,tSNP1,tSNP2)/n.obs
		proba.SNP <- table(SNP1,SNP2)/n.obs
		proba.tSNP <- table(tSNP1,tSNP2)/n.obs
  	} else{
		proba.marg <- getProba.Paire(proba.jointe)
		proba.SNP <- proba.marg$p.AB
		proba.tSNP <- proba.marg$p.CD
	}
  
  I <- 0
  for (s1 in 1:2){
    for (s2 in 1:2){
      for (ts1 in 1:2){
        for (ts2 in 1:2){
        	if (is.null(pj.test)){
	          if (proba.jointe[s1,s2,ts1,ts2] != 0){
    	        I <- I+proba.jointe[s1,s2,ts1,ts2]*log(proba.jointe[s1,s2,ts1,ts2]/(proba.SNP[s1,s2]*proba.tSNP[ts1,ts2]))
        	  }
        	}else{
	          if (proba.jointe[2*(ts1-1)+ts2,2*(s1-1)+s2] != 0){
				I <- I+proba.jointe[2*(ts1-1)+ts2,2*(s1-1)+s2]*log(proba.jointe[2*(ts1-1)+ts2,2*(s1-1)+s2]/(proba.SNP[s1,s2]*proba.tSNP[ts1,ts2]))
				}
			}
        }
      }
    }
  }
  return(I)
}

###
## getH.allelic: estimate the Information between two pairs of SNPs in genotypic format
###

getH.allelic <- function(proba.Paire=NULL,i=NA,j=NA){
	if (is.null(proba.Paire)){
		H <- getI.allelic(SNP1=i,SNP2=j,tSNP1=i,tSNP2=j)
	}
	else{
		H <- 0
		for (s1 in 1:2){
			for (s2 in 1:2){
				if (proba.Paire[s1,s2] != 0){
					H <- H-proba.Paire[s1,s2]*log(proba.Paire[s1,s2])
				}
			}
		}
	}
	return(H)
}

getMI.allelic <- function(proba.jointe=NULL,SNP1=NA,SNP2=NA,tSNP1=NA,tSNP2=NA,Hi=NULL,Hj=NULL){
	if (is.null(proba.jointe)){
		I <- getI.allelic(SNP1=SNP1,SNP2=SNP2,tSNP1=tSNP1,tSNP2=tSNP2)
		if (is.null(Hi)){
			Hi <- getH.allelic(i=SNP1,j=SNP2)
		}
		if (is.null(Hj)){
			Hj <- getH.allelic(i=tSNP1,j=tSNP2)
		}
	}
	else{
		I <- getI.allelic(proba.jointe=proba.jointe)
		pP <- getProba.Paire(proba.jointe)
		if (is.null(Hi)){
			Hi <- getH.allelic(proba.Paire=pP$p.AB)
		}
		if (is.null(Hj)){
			Hj <- getH.allelic(proba.Paire=pP$p.CD)
		}
	}
	return(I/sqrt(Hi*Hj))
}


###
## getI.geno: estimate the Information between two pairs of SNPs in genotypic format
###

getI.geno <- function(SNP1,SNP2,tSNP1,tSNP2){
  n.obs <- length(SNP1)
  SNP1 <- factor(SNP1,levels=c(0,1,2))
  SNP2 <- factor(SNP2,levels=c(0,1,2))
  tSNP1 <- factor(tSNP1,levels=c(0,1,2))
  tSNP2 <- factor(tSNP2,levels=c(0,1,2))
  
  # /n.obs et pas /n 
  proba.jointe <- table(SNP1,SNP2,tSNP1,tSNP2)/n.obs
  proba.SNP <- table(SNP1,SNP2)/n.obs
  proba.tSNP <- table(tSNP1,tSNP2)/n.obs
  
  I <- 0
  for (s1 in 1:3){
    for (s2 in 1:3){
      for (ts1 in 1:3){
        for (ts2 in 1:3){
          if (proba.jointe[s1,s2,ts1,ts2] != 0){
            I <- I+proba.jointe[s1,s2,ts1,ts2]*log(proba.jointe[s1,s2,ts1,ts2]/(proba.SNP[s1,s2]*proba.tSNP[ts1,ts2]))
          }
        }
      }
    }
  }
  return(I)
}

###
## getH.geno: estimate the Entropy between two pairs of SNPs in genotypic format
###
getH.geno <- function(i,j){
  return(getI.geno(i,j,i,j))
}

###
## getMI.geno: estimate the Mutual Information between two pairs of SNPs in genotypic format
###
getMI.geno <- function(SNP1,SNP2,tSNP1,tSNP2,Hi=NULL,Hj=NULL){
  I <- getI.geno(SNP1=SNP1,SNP2=SNP2,tSNP1=tSNP1,tSNP2=tSNP2)
  if (is.null(Hi)){
  	Hi <- getH.geno(SNP1,SNP2)}
  if (is.null(Hj)){
  	Hj <- getH.geno(tSNP1,tSNP2)}
  return(I/sqrt(Hi*Hj))
}

