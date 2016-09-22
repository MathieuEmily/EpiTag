getMatMI.1Region <- function(data,allelic=FALSE){
  n.snp <- ncol(data)
  mat.MI <- array(NA,dim=c(n.snp,n.snp,n.snp,n.snp))
  nPairs <- n.snp*(n.snp-1)/2
  compt <- 1
  for (s1 in 1:(n.snp-1)){
    for (s2 in (s1+1):n.snp){
      cat("Pairs", compt, "/",nPairs,"running")
      for (t.s1 in 1:(n.snp-1)){
        for (t.s2 in (t.s1+1):n.snp){
        	if (allelic){
        		mat.MI[s1,s2,t.s1,t.s2] <- mat.MI[t.s1,t.s2,s1,s2] <- getMI.allelic(data[,s1],data[,s2],data[,t.s1],data[,t.s2])	
        	}
        	else{
	          mat.MI[s1,s2,t.s1,t.s2] <- mat.MI[t.s1,t.s2,s1,s2] <- getMI.geno(data[,s1],data[,s2],data[,t.s1],data[,t.s2])
	          }
        }
      }
      cat("....Done \n")
      compt <- compt+1
    }
  }
  dimnames(mat.MI) <- list(names(data),names(data),names(data),names(data))
  mat.MI <- list(mat.MI=mat.MI,nb.locus=1)
  class(mat.MI) <- "MIMatrix"
  return(mat.MI)
}

getMatMI.2Regions <- function(data1,data2,allelic=FALSE){
  if (nrow(data2) != nrow(data1)){
  	print("The two datasets does not have the same number of individuals")
  	return(NA)
  }
  n.snp1 <- ncol(data1)
  n.snp2 <- ncol(data2)
  mat.MI <- array(NA,dim=c(n.snp1,n.snp2,n.snp1,n.snp2))
  nPairs <- n.snp1*n.snp2
  compt <- 1
  for (s1 in 1:n.snp1){
    for (s2 in 1:n.snp2){
      cat("Pairs", compt, "/",nPairs,"running")
      for (t.s1 in 1:n.snp1){
        for (t.s2 in 1:n.snp2){
        	if (allelic){
        		mat.MI[s1,s2,t.s1,t.s2] <- getMI.allelic(SNP1=data1[,s1],SNP2=data2[,s2],tSNP1=data1[,t.s1],tSNP2=data2[,t.s2])
        	}
        	else{
	          mat.MI[s1,s2,t.s1,t.s2] <- getMI.geno(data1[,s1],data2[,s2],data1[,t.s1],data2[,t.s2])
	          }
        }
      }
      cat("....Done \n")
      compt <- compt+1
    }
  }
  dimnames(mat.MI) <- list(names(data1),names(data2),names(data1),names(data2))
  mat.MI <- list(mat.MI=mat.MI,nb.locus=2)
  class(mat.MI) <- "MIMatrix"
  return(mat.MI)
}


getMatMI <- function(Region1,Region2=NULL,allelic=FALSE){
	if (is.null(Region2)){
		return(getMatMI.1Region(data=Region1,allelic=allelic))
	} else {
		return(getMatMI.2Regions(data1=Region1,data2=Region2,allelic=allelic))	
	}
}