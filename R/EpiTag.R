EpiTag.1Region <- function(mat.MI,threshold=0.8){
  n.snp=dim(mat.MI)[1]
  mat.MI.tmp <- mat.MI
  tag.list <- list()
  bin.list <- list()
  compt <- 1
  while(sum(unlist(lapply(bin.list,nrow))) < (n.snp*(n.snp-1)/2)){
    cat("n.Pair.Tag:",compt,"- n.Pair.Tagged",sum(unlist(lapply(bin.list,nrow))),"\n")
    tag.snp1 <- NA
    tag.snp2 <- NA
    max <- 1
    mean.MI <- 0
    for (s1 in 1:(n.snp-1)){
      for (s2 in (s1+1):n.snp){
        tmp.count <- sum(mat.MI.tmp[s1,s2,,] >= threshold,na.rm=TRUE)
        #print(tmp.count)
        if (tmp.count > max){
          max <- tmp.count
          w <- which(mat.MI.tmp[s1,s2,,] >= threshold)
          mean.MI <- mean(as.matrix(mat.MI.tmp[s1,s2,,])[w])
          tag.snp1 <- s1
          tag.snp2 <- s2
          tagged <- lapply(1:n.snp,FUN=function(i){which(mat.MI.tmp[s1,s2,i,] > threshold)})
          tagged.pair <- NULL
          for (ts1 in 1:n.snp){
            tmp <- cbind(rep(ts1,length(tagged[[ts1]])),tagged[[ts1]])
            tagged.pair <- rbind(tagged.pair,tmp)
          }
        }
        if (tmp.count==max){
          w <- which(mat.MI.tmp[s1,s2,,] >= threshold)
          mean.MI.tmp <- mean(as.matrix(mat.MI.tmp[s1,s2,,])[w])
          if (mean.MI.tmp > mean.MI){
            mean.MI <- mean.MI.tmp
            tag.snp1 <- s1
            tag.snp2 <- s2
            tagged <- lapply(1:n.snp,FUN=function(i){which(mat.MI.tmp[s1,s2,i,] > threshold)})
            tagged.pair <- NULL
            for (ts1 in 1:n.snp){
              tmp <- cbind(rep(ts1,length(tagged[[ts1]])),tagged[[ts1]])
              tagged.pair <- rbind(tagged.pair,tmp)
            }
          }
        }
        #	cat("tag1",tag.snp1,"tag2",tag.snp2,"mean",mean.MI,"\n")
      }
    }
    for (ii in 1:nrow(tagged.pair)){
      mat.MI.tmp[tagged.pair[ii,1],tagged.pair[ii,2],,] <- 0
      mat.MI.tmp[,,tagged.pair[ii,1],tagged.pair[ii,2]] <- 0
    }
    tag.list[[compt]] <- c(tag.snp1,tag.snp2)
    bin.list[[compt]] <- tagged.pair
    compt <- compt+1
  }
  tags <- unique(unlist(tag.list))
  return(list(tagSNPs.R1=dimnames(mat.MI)[[1]][tags],pairs.tag=tag.list,pairs.bins=bin.list))
}


###
## EpiTag.2Regions: compute the tagSNPs and the bins from a MutualInformation matrix and a threshold
###
EpiTag.2Regions <- function(mat.MI,threshold=0.8){
  n.snp1=dim(mat.MI)[1]
  n.snp2=dim(mat.MI)[2]
  mat.MI.tmp <- mat.MI
  tag.list.R1 <- list()
  tag.list.R2 <- list()
  tag.list <- list()
  bin.list <- list()
  compt <- 1
  while(sum(unlist(lapply(bin.list,nrow))) < (n.snp1*n.snp2)){
    cat("n.Pair.Tag:",compt,"- n.Pair.Tagged",sum(unlist(lapply(bin.list,nrow))),"\n")
    tag.snp1 <- NA
    tag.snp2 <- NA
    max <- 1
    mean.MI <- 0
    for (s1 in 1:n.snp1){
      for (s2 in 1:n.snp2){
        tmp.count <- sum(mat.MI.tmp[s1,s2,,] >= threshold,na.rm=TRUE)
        #print(tmp.count)
        if (tmp.count > max){
          max <- tmp.count
          w <- which(mat.MI.tmp[s1,s2,,] >= threshold)
          mean.MI <- mean(as.matrix(mat.MI.tmp[s1,s2,,])[w])
          tag.snp1 <- s1
          tag.snp2 <- s2
          tagged <- lapply(1:n.snp1,FUN=function(i){which(mat.MI.tmp[s1,s2,i,] > threshold)})
          tagged.pair <- NULL
          for (ts1 in 1:n.snp1){
            tmp <- cbind(rep(ts1,length(tagged[[ts1]])),tagged[[ts1]])
            tagged.pair <- rbind(tagged.pair,tmp)
          }
        }
        if (tmp.count==max){
          w <- which(mat.MI.tmp[s1,s2,,] >= threshold)
          mean.MI.tmp <- mean(as.matrix(mat.MI.tmp[s1,s2,,])[w])
          if (mean.MI.tmp > mean.MI){
            mean.MI <- mean.MI.tmp
            tag.snp1 <- s1
            tag.snp2 <- s2
            tagged <- lapply(1:n.snp1,FUN=function(i){which(mat.MI.tmp[s1,s2,i,] > threshold)})
            tagged.pair <- NULL
            for (ts1 in 1:n.snp1){
              tmp <- cbind(rep(ts1,length(tagged[[ts1]])),tagged[[ts1]])
              tagged.pair <- rbind(tagged.pair,tmp)
            }
          }
        }
        #	cat("tag1",tag.snp1,"tag2",tag.snp2,"mean",mean.MI,"\n")
      }
    }
    for (ii in 1:nrow(tagged.pair)){
      mat.MI.tmp[tagged.pair[ii,1],tagged.pair[ii,2],,] <- 0
      mat.MI.tmp[,,tagged.pair[ii,1],tagged.pair[ii,2]] <- 0
    }
    tag.list.R1[[compt]] <- tag.snp1
    tag.list.R2[[compt]] <- tag.snp2
    tag.list[[compt]] <- c(tag.snp1,tag.snp2)
    bin.list[[compt]] <- tagged.pair
#    cat(tag.snp1,"-",tag.snp2,"\n")
 #   print(bin.list)
  #  print(mat.MI.tmp[2,2,,])
    compt <- compt+1
  }
  tags.R1 <- unique(unlist(tag.list.R1))
  tags.R2 <- unique(unlist(tag.list.R2))
  return(list(tagSNPs.R1=dimnames(mat.MI)[[1]][tags.R1],tagSNPs.R2=dimnames(tmpMatMI$mat.MI)[[1]][tags.R2],pairs.tag=tag.list,pairs.bins=bin.list))
}


EpiTag <- function(Region1=NULL,Region2=NULL,mat.MI=NULL,threshold=0.8,allelic=FALSE){
	if (!is.null(Region1) & !is.null(Region2) & !is.null(mat.MI)){
		stop("At least of one the three arguments Region1, Region2 and mat.MI must be supplied")
	}
	if (!is.null(mat.MI)){
		if (class(mat.MI)!="MIMatrix"){
			stop("mat.MI is not of class MIMatrix")
		} else{
			if (mat.MI$nb.locus==1){
				return(EpiTag.1Region(mat.MI$mat.MI,threshold= threshold))
			} else {
				return(EpiTag.2Regions(mat.MI$mat.MI,threshold= threshold))
			}
		}
	} else {
		if (is.null(Region2)){
			mat.MI <- getMatMI(Region1=Region1,allelic=allelic)
			return(EpiTag.1Region(mat.MI$mat.MI,threshold= threshold))
		} else {
			mat.MI <- getMatMI(Region1=Region1,Region2=Region2,allelic=allelic)
			return(EpiTag.2Regions(mat.MI$mat.MI,threshold= threshold))
		}
	}
	
}

