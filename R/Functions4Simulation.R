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


getI.allelic.4simu <- function(proba.jointe){
	
  # /n.obs et pas /n 
#	proba.jointe <- table(SNP1,SNP2,tSNP1,tSNP2)/n.obs
	proba.marg <- getProba.Paire(proba.jointe)
	proba.SNP <- proba.marg$p.AB
	proba.tSNP <- proba.marg$p.CD

	I <- 0
	for (s1 in 1:2){
		for (s2 in 1:2){
			for (ts1 in 1:2){
				for (ts2 in 1:2){
					if (proba.jointe[2*(ts1-1)+ts2,2*(s1-1)+s2] != 0){
						I <- I+proba.jointe[2*(ts1-1)+ts2,2*(s1-1)+s2]*log(proba.jointe[2*(ts1-1)+ts2,2*(s1-1)+s2]/(proba.SNP[s1,s2]*proba.tSNP[ts1,ts2]))
					}
				}
			}
		}
	}
	return(I)
}

getH.allelic.4simu <- function(proba.Paire){
	
	H <- 0
	for (s1 in 1:2){
		for (s2 in 1:2){
			if (proba.Paire[s1,s2] != 0){
				H <- H-proba.Paire[s1,s2]*log(proba.Paire[s1,s2])
			}
		}
	}
	return(H)
}


getMI.allelic.4simu <- function(proba.jointe,Hi=NULL,Hj=NULL){
	I <- getI.allelic.4simu(proba.jointe)
	pP <- getProba.Paire(proba.jointe)
	if (is.null(Hi)){
		Hi <- getH.allelic.4simu(pP$p.AB)
		}
	if (is.null(Hj)){
		Hj <- getH.allelic.4simu(pP$p.CD)
		}
	return(I/sqrt(Hi*Hj))
}



get.pj.MaxMin <- function(p.A=0.6,p.B=0.6,p.C=0.6,p.D=0.6,p.E=0.6,p.F=0.6,r2.AB=0,r2.CD=0,r2.EF=0,r2.AC=0.9,r2.BD=0.9,r2.AE=0.9,r2.BF=0.9,n.sim.pj=100){
	pj.max <- NULL
	pj.min <- NULL
	my.max <- 0
	my.min <- 1

	myMI <- rep(NA,times=n.sim.pj)
	pj.list <- list()

	for (i in 1:n.sim.pj){
		if (i%%50 ==0){
			cat("Sim",i,"sur",n.sim.pj,"\n")
		}

		#pj <- get.probajointe.2(p.A=p.A,p.B=p.B,p.C=p.C,p.D=p.D,r2AC=r2AC,r2BD=r2BD,r2AB=r2AB,r2CD=r2CD,r.names=c("cd","cD","Cd","CD"),c.names=c("ab","aB","Ab","AB"))
		# modif Chlo:
		pj <- get.probajointe.2(p.A=p.A,p.B=p.B,p.C=p.C,p.D=p.D,r2AC=r2.AC,r2BD=r2.BD,r2AB=r2.AB,r2CD=r2.CD,r.names=c("cd","cD","Cd","CD"),c.names=c("ab","aB","Ab","AB"))
		#check.constraints(pj)
		tmp <- getMI.allelic.4simu(pj)
		myMI[i] <- tmp
		pj.list[[i]] <- pj
		
		if (tmp > my.max){
			pj.max <- pj
			my.max <- tmp
		}
		if (tmp < my.min){
			pj.min <- pj
			my.min <- tmp
		}
	}	
	return(list(MI=myMI,pj.list=pj.list,pj.min=pj.min,pj.max=pj.max))
}


get.probajointe.2 <- function(p.A=0.6,p.B=0.6,p.C=0.6,p.D=0.6,r2AC=0.8,r2BD=0.8,r2AB=0.8,r2CD=0.8,first=0,c.names=c("ab","aB","Ab","AB"),r.names=c("cd","cD","Cd","CD")){
	res <- list(pj=NULL,control=FALSE)
	compt <- 0
	while(!res$control){
		res <- generate.probajointe.2(p.A=p.A,p.B=p.B,p.C=p.C,p.D=p.D,r2AC=r2AC,r2BD=r2BD,r2AB=r2AB,r2CD=r2CD,first=0,c.names=c.names,r.names=r.names)
			compt <- compt+1
			#cat("Essai",compt,"\n")
	}
	return(res$pj)
}

generate.probajointe.2 <- function(p.A=0.6,p.B=0.6,p.C=0.6,p.D=0.6,r2AC=0.8,r2BD=0.8,r2AB=0.8,r2CD=0.8,first=0,c.names=c("ab","aB","Ab","AB"),r.names=c("cd","cD","Cd","CD")){
	proba.jointe <- matrix(0,ncol=4,nrow=4)
	#print(r.names)
	proba.jointe <- data.frame(proba.jointe,row.names=r.names)
	names(proba.jointe) <- c.names

	p.a <- 1-p.A
	p.b <- 1-p.B
	p.c <- 1-p.C
	p.d <- 1-p.D

	p.AC <- p.A*p.C+sqrt(r2AC*p.A*p.C*p.a*p.c)
	p.Ac <- p.A*p.c-sqrt(r2AC*p.A*p.C*p.a*p.c)
	p.aC <- p.a*p.C-sqrt(r2AC*p.A*p.C*p.a*p.c)
	p.ac <- p.a*p.c+sqrt(r2AC*p.A*p.C*p.a*p.c)

	p.BD <- p.B*p.D+sqrt(r2BD*p.B*p.D*p.b*p.d)
	p.Bd <- p.B*p.d-sqrt(r2BD*p.B*p.D*p.b*p.d)
	p.bD <- p.b*p.D-sqrt(r2BD*p.B*p.D*p.b*p.d)
	p.bd <- p.b*p.d+sqrt(r2BD*p.B*p.D*p.b*p.d)

	p.AB <- p.A*p.B+sqrt(r2AB*p.A*p.B*p.a*p.b)
	p.Ab <- p.A*p.b-sqrt(r2AB*p.A*p.B*p.a*p.b)
	p.aB <- p.a*p.B-sqrt(r2AB*p.A*p.B*p.a*p.b)
	p.ab <- p.a*p.b+sqrt(r2AB*p.A*p.B*p.a*p.b)

	p.CD <- p.C*p.D+sqrt(r2CD*p.C*p.D*p.c*p.d)
	p.Cd <- p.C*p.d-sqrt(r2CD*p.C*p.D*p.c*p.d)
	p.cD <- p.c*p.D-sqrt(r2CD*p.C*p.D*p.c*p.d)
	p.cd <- p.c*p.d+sqrt(r2CD*p.C*p.D*p.c*p.d)

#	cat("0 - p.ac",p.ac,"-p.BD:",p.BD,"-p.aB:",p.aB,"-p.cD:",p.cD,"\n")

	# p.abcd - 1
#	tmp.min <- max(0,p.ac-p.aB-p.ac,p.ac-p.aB-p.bd,p.ac-p.aB-p.ab,p.ac-p.aB-p.cD)
	tmp.min <- first
	tmp.max <- min(c(p.ac,p.bd,p.ab,p.cd))
	tmp.p <- runif(1,tmp.min,tmp.max)
	proba.jointe[1,1] <- tmp.p
	p.ac <- p.ac-tmp.p
	p.bd <- p.bd-tmp.p
	p.ab <- p.ab-tmp.p
	p.cd <- p.cd-tmp.p

#	cat("1 (abcd -",tmp.p,tmp.min,tmp.max,")- p.ac",p.ac,"-p.BD:",p.BD,"-p.aB:",p.aB,"-p.cD:",p.cD,"\n")
	#cat("1 (abcd) - p.ac",p.ac,"- p.bd",p.bd,"- p.ab",p.ab,"- p.cd",p.cd,"\n")
	ppp <- c(p.ac,p.bd,p.ab,p.cd)
	if (any(ppp < -1e-10)){
		return(list(control=FALSE))
	}

	# p.abcD - 2
#	tmp.min <- max(c(p.ac-p.aB))
	tmp.min <- 0
	tmp.max <- min(c(p.ac,p.bD,p.ab,p.cD))
	tmp.p <- runif(1,tmp.min,tmp.max)
	proba.jointe[2,1] <- tmp.p
	p.ac <- p.ac-tmp.p
	p.bD <- p.bD-tmp.p
	p.ab <- p.ab-tmp.p
	p.cD <- p.cD-tmp.p

#	cat("2 (abcD -",tmp.p,tmp.min,tmp.max,")- p.ac:",p.ac,"-p.BD:",p.BD,"-p.aB:",p.aB,"-p.cD:",p.cD,"\n")
	#cat("2 (abcD) - p.ac",p.ac,"- p.bD",p.bD,"- p.ab",p.ab,"- p.cD",p.cD,"\n")
	ppp <- c(p.ac,p.bD,p.ab,p.cD)
	if (any(ppp < -1e-10)){
		return(list(control=FALSE))
	}

	# p.aBcd - 3
#	tmp.min <- max(0,c(p.ac-p.ac,p.ac-p.BD,p.ac-p.aB,p.ac-p.cD))
	tmp.min <- 0
	tmp.max <- min(c(p.ac,p.Bd,p.aB,p.cd))
	tmp.p <- runif(1,tmp.min,tmp.max)
	proba.jointe[1,2] <- tmp.p
	p.ac <- p.ac-tmp.p
	p.Bd <- p.Bd-tmp.p
	p.aB <- p.aB-tmp.p
	p.cd <- p.cd-tmp.p
#	cat("3 (aBcd -",tmp.p,tmp.min,tmp.max,")- p.ac:",p.ac,"-p.BD:",p.BD,"-p.aB:",p.aB,"-p.cD:",p.cD,"\n")
	#cat("3 (aBcd) - p.ac",p.ac,"- p.Bd",p.Bd,"- p.aB",p.aB,"- p.cd",p.cd,"\n")
	ppp <- c(p.ac,p.Bd,p.aB,p.cd)
	if (any(ppp < -1e-10)){
		return(list(control=FALSE))
	}

	# p.aBcD - 4 - Fixed
	tmp.p <- p.ac
	proba.jointe[2,2] <- tmp.p
	p.ac <- p.ac-tmp.p
	p.BD <- p.BD-tmp.p
	p.aB <- p.aB-tmp.p
	p.cD <- p.cD-tmp.p
#	cat("4 (aBcD -",tmp.p,")- p.ac:",p.ac,"-p.BD:",p.BD,"-p.aB:",p.aB,"-p.cD:",p.cD,"\n")
	#cat("4 (aBcD) - p.ac",p.ac,"- p.BD",p.BD,"- p.aB",p.aB,"- p.cD",p.cD,"\n")
	ppp <- c(p.ac,p.BD,p.aB,p.cD)
	if (any(ppp < -1e-10)){
		return(list(control=FALSE))
	}

	# p.abCd - 5
	tmp.max <- min(c(p.aC,p.bd,p.ab,p.Cd))
	tmp.min <- 0
	tmp.p <- runif(1,tmp.min,tmp.max)
	proba.jointe[3,1] <- tmp.p
	p.aC <- p.aC-tmp.p
	p.bd <- p.bd-tmp.p
	p.ab <- p.ab-tmp.p
	p.Cd <- p.Cd-tmp.p
	#cat("5 (abCd) - p.aC",p.aC,"- p.bd",p.bd,"- p.ab",p.aB,"- p.Cd",p.cD,"\n")
	ppp <- c(p.aC,p.bd,p.ab,p.Cd)
	if (any(ppp < -1e-10)){
		return(list(control=FALSE))
	}

	# p.abCD - 6 - Fixed

	tmp.p <- p.ab
	proba.jointe[4,1] <- tmp.p
	p.aC <- p.aC-tmp.p
	p.bD <- p.bD-tmp.p
	p.ab <- p.ab-tmp.p
	p.CD <- p.CD-tmp.p
	#cat("6 (abCD) - p.aC",p.aC,"- p.bD",p.bD,"- p.ab",p.ab,"- p.CD",p.CD,"\n")
	ppp <- c(p.aC,p.bD,p.ab,p.CD)
	if (any(ppp < -1e-10)){
		return(list(control=FALSE))
	}

	# p.aBCd - 7
	tmp.max <- min(c(p.aC,p.Bd,p.aB,p.Cd))
	tmp.min <- 0
	tmp.p <- runif(1,tmp.min,tmp.max)
	proba.jointe[3,2] <- tmp.p
	p.aC <- p.aC-tmp.p
	p.Bd <- p.Bd-tmp.p
	p.aB <- p.aB-tmp.p
	p.Cd <- p.Cd-tmp.p
	#cat("7 (aBCd) - p.aC",p.aC,"- p.Bd",p.Bd,"- p.aB",p.aB,"- p.Cd",p.Cd,"\n")
	ppp <- c(p.aC,p.Bd,p.aB,p.Cd)
	if (any(ppp < -1e-10)){
		return(list(control=FALSE))
	}

	# p.aBCD - 8 - Fixed
	tmp.p <- p.aC
	proba.jointe[4,2] <- tmp.p
	p.aC <- p.aC-tmp.p
	p.BD <- p.BD-tmp.p
	p.aB <- p.aB-tmp.p
	p.CD <- p.CD-tmp.p
	#cat("8 (aBCD) - p.aC",p.aC,"- p.BD",p.BD,"- p.aB",p.aB,"- p.CD",p.CD,"\n")
	ppp <- c(p.aC,p.BD,p.aB,p.CD)
	if (any(ppp < -1e-10)){
		return(list(control=FALSE))
	}

	# p.Abcd - 9
	tmp.max <- min(c(p.Ac,p.bd,p.Ab,p.cd))
	tmp.min <- 0
	tmp.p <- runif(1,tmp.min,tmp.max)
	proba.jointe[1,3] <- tmp.p
	p.Ac <- p.Ac-tmp.p
	p.bd <- p.bd-tmp.p
	p.Ab <- p.Ab-tmp.p
	p.cd <- p.cd-tmp.p
	#cat("9 (Abcd) - p.Ac",p.Ac,"- p.bd",p.bd,"- p.Ab",p.Ab,"- p.cd",p.cd,"\n")
	ppp <- c(p.Ac,p.bd,p.Ab,p.cd)
	if (any(ppp < -1e-10)){
		return(list(control=FALSE))
	}

	# p.AbCd - 10 - Fixed
	tmp.p <- p.bd
	proba.jointe[3,3] <- tmp.p
	p.AC <- p.AC-tmp.p
	p.bd <- p.bd-tmp.p
	p.Ab <- p.Ab-tmp.p
	p.Cd <- p.Cd-tmp.p
	#cat("10 (AbCD) - p.AC",p.AC,"- p.bd",p.bd,"- p.Ab",p.Ab,"- p.Cd",p.Cd,"\n")
	ppp <- c(p.AC,p.bd,p.Ab,p.Cd)
	if (any(ppp < -1e-10)){
		return(list(control=FALSE))
	}

	# p.ABcd - 11 - Fixed
	tmp.p <- p.cd
	proba.jointe[1,4] <- tmp.p
	p.Ac <- p.Ac-tmp.p
	p.Bd <- p.Bd-tmp.p
	p.AB <- p.AB-tmp.p
	p.cd <- p.cd-tmp.p
	#cat("11 (ABcd) - p.Ac",p.Ac,"- p.Bd",p.Bd,"- p.AB",p.AB,"- p.cd",p.cd,"\n")
	ppp <- c(p.Ac,p.Bd,p.AB,p.Cd)
	if (any(ppp < -1e-10)){
		return(list(control=FALSE))
	}

	# p.ABCd - 12 - Fixed
	tmp.p <- p.Bd
	proba.jointe[3,4] <- tmp.p
	p.AC <- p.AC-tmp.p
	p.Bd <- p.Bd-tmp.p
	p.AB <- p.AB-tmp.p
	p.Cd <- p.Cd-tmp.p
	#cat("12 (ABCd) - p.AC",p.AC,"- p.Bd",p.Bd,"- p.AB",p.AB,"- p.Cd",p.Cd,"\n")
	ppp <- c(p.AC,p.Bd,p.Ab,p.Cd)
	if (any(ppp < -1e-10)){
		return(list(control=FALSE))
	}

	# p.AbcD - 13
	tmp.max <- min(c(p.Ac,p.bD,p.Ab,p.cD))
	tmp.min <- 0
	tmp.p <- runif(1,tmp.min,tmp.max)
	proba.jointe[2,3] <- tmp.p
	p.Ac <- p.Ac-tmp.p
	p.bD <- p.bD-tmp.p
	p.Ab <- p.Ab-tmp.p
	p.cD <- p.cD-tmp.p
	#cat("13 (AbcD) - p.Ac",p.Ac,"- p.bD",p.bD,"- p.Ab",p.Ab,"- p.cD",p.cD,"\n")
	ppp <- c(p.Ac,p.bD,p.Ab,p.cD)
	if (any(ppp < -1e-10)){
		return(list(control=FALSE))
	}

	# p.ABcD - 14 - Fixed
	tmp.p <- p.Ac
	proba.jointe[2,4] <- tmp.p
	p.Ac <- p.Ac-tmp.p
	p.BD <- p.BD-tmp.p
	p.AB <- p.AB-tmp.p
	p.cD <- p.cD-tmp.p
	#cat("14 (ABcD) - p.Ac",p.Ac,"- p.BD",p.BD,"- p.AB",p.AB,"- p.cD",p.cD,"\n")
	ppp <- c(p.Ac,p.BD,p.AB,p.cD)
	if (any(ppp < -1e-10)){
		return(list(control=FALSE))
	}

	# p.AbCD - 15 - Fixed
	tmp.p <- p.bD
	proba.jointe[4,3] <- tmp.p
	p.AC <- p.AC-tmp.p
	p.bD <- p.bD-tmp.p
	p.Ab <- p.Ab-tmp.p
	p.CD <- p.CD-tmp.p
	#cat("15 (AbCD) - p.AC",p.AC,"- p.bD",p.bD,"- p.Ab",p.Ab,"- p.CD",p.CD,"\n")
	ppp <- c(p.AC,p.bD,p.Ab,p.CD)
	if (any(ppp < -1e-10)){
		return(list(control=FALSE))
	}

	# p.ABCD - 16 - Fixed
	tmp.p <- p.AC
	proba.jointe[4,4] <- tmp.p
	p.AC <- p.AC-tmp.p
	p.BD <- p.BD-tmp.p
	p.AB <- p.AB-tmp.p
	p.CD <- p.CD-tmp.p
	#cat("16 (ABCD) - p.AC",p.AC,"- p.BD",p.BD,"- p.AB",p.AB,"- p.CD",p.CD,"\n")
	ppp <- c(p.AC,p.BD,p.AB,p.CD)
	if (any(ppp < -1e-10)){
		return(list(control=FALSE))
	}

	return(list(pj=proba.jointe,control=TRUE))
}

generate.pj.triplet <- function(pj.ab.cd,pj.ab.ef){
	proba.jointe <- array(NA,c(4,4,4))
	all.proba <- list(pj.ab.cd=pj.ab.cd,pj.ab.ef=pj.ab.ef,proba.jointe=proba.jointe)
	while(sum(is.na(all.proba$proba.jointe))){
		my.min <- get.Min(all.proba)
		all.proba <- respect.constraint(all.proba=all.proba,id=my.min$id,row=my.min$row,col=my.min$col)
	}
	return(all.proba)
}

get.Min <- function(all.proba){
	all.proba.ori <- all.proba
	all.proba$pj.ab.cd[all.proba$pj.ab.cd<1e-15] <- 1.1
	all.proba$pj.ab.ef[all.proba$pj.ab.ef<1e-15] <- 1.1
	#print(all.proba.ori$pj.ab.cd)
	#print(all.proba$pj.ab.cd)
	#print(all.proba.ori$pj.ab.ef)
	#print(all.proba$pj.ab.ef)
	id <- which.min(c(min(all.proba$pj.ab.cd,na.rm=TRUE),min(all.proba$pj.ab.ef,na.rm=TRUE)))
	if (id==1){
		col1 <- which.min(apply(all.proba$pj.ab.cd,2,min))
		row1 <- which.min(all.proba$pj.ab.cd[,col1])
		res <- list(id="ab.cd",row=row1,col=col1)
	}
	if (id==2){
		col2 <- which.min(apply(all.proba$pj.ab.ef,2,min))
		row2 <- which.min(all.proba$pj.ab.ef[,col2])
		res <- list(id="ab.ef",row=row2,col=col2)
	}
	return(res)
}

respect.constraint <- function(all.proba,id,row,col){
	if (id=="ab.cd"){
		p <- all.proba$pj.ab.cd[row,col]
		id.b <- (col-1)%%2
		id.a <- as.numeric(col>2)
		id.c <- (row-1)%%2
		id.d <- as.numeric(row>2)
		id.ab <- col
		id.cd <- row
		#cat(id,"-",id.ab,"-",id.cd,"\n")
		id.e <- 0;id.f <- 0;id.ef <- 2*id.e+id.f+1
		w.ef <- all.proba$proba.jointe[id.ef,id.cd,id.ab]
		if (is.na(w.ef)){
			min.ef <- min(all.proba$pj.ab.cd[id.cd,id.ab],all.proba$pj.ab.ef[id.ef,id.ab])
		}
		else{min.ef <- 1.1}
		id.e <- 0;id.f <- 1;id.ef <- 2*id.e+id.f+1
		w.eF <- all.proba$proba.jointe[id.ef,id.cd,id.ab]
		if (is.na(w.eF)){
			min.eF <- min(all.proba$pj.ab.cd[id.cd,id.ab],all.proba$pj.ab.ef[id.ef,id.ab])
		}
		else{min.eF <- 1.1}
		id.e <- 1;id.f <- 0;id.ef <- 2*id.e+id.f+1
		w.Ef <- all.proba$proba.jointe[id.ef,id.cd,id.ab]
		if (is.na(w.Ef)){
			min.Ef <- min(all.proba$pj.ab.cd[id.cd,id.ab],all.proba$pj.ab.ef[id.ef,id.ab])
		}
		else{min.Ef <- 1.1}
		id.e <- 1;id.f <- 1;id.ef <- 2*id.e+id.f+1
		w.EF <- all.proba$proba.jointe[id.ef,id.cd,id.ab]
		if (is.na(w.EF)){
			min.EF <- min(all.proba$pj.ab.cd[id.cd,id.ab],all.proba$pj.ab.ef[id.ef,id.ab])
		}
		else{min.EF <- 1.1}
		mins <- c(min.ef,min.eF,min.Ef,min.EF)
		#print(mins)
		n.mins <- sum(mins <=1)
		if (n.mins > 1){
			for (i in 1:(n.mins-1)){
				ww <- order(mins)[i]
				id.ef <- ww
				vec <- c(all.proba$pj.ab.cd[id.cd,id.ab],all.proba$pj.ab.ef[id.ef,id.ab])
				tmp.max <- min(vec)
				tmp.p <- runif(1,0,tmp.max)
				#cat("[",0,"-",tmp.max,"] - ",tmp.p,"\n")
				all.proba$proba.jointe[id.ef,id.cd,id.ab] <- tmp.p
				all.proba$pj.ab.cd[id.cd,id.ab] <- all.proba$pj.ab.cd[id.cd,id.ab]-tmp.p
				all.proba$pj.ab.ef[id.ef,id.ab] <- all.proba$pj.ab.ef[id.ef,id.ab]-tmp.p
			}
		}
		ww <- order(mins)[n.mins]
		id.ef <- ww
		tmp.p <- all.proba$pj.ab.cd[id.cd,id.ab]
		#cat("[",0,"-",tmp.p,"] - ",tmp.p,"\n")
		all.proba$proba.jointe[id.ef,id.cd,id.ab] <- tmp.p
		all.proba$pj.ab.cd[id.cd,id.ab] <- all.proba$pj.ab.cd[id.cd,id.ab]-tmp.p
		all.proba$pj.ab.ef[id.ef,id.ab] <- all.proba$pj.ab.ef[id.ef,id.ab]-tmp.p
	}
	
	##
	if (id=="ab.ef"){
		p <- all.proba$pj.ab.ef[row,col]
		id.b <- (col-1)%%2
		id.a <- as.numeric(col>2)
		id.f <- (row-1)%%2
		id.e <- as.numeric(row>2)
		id.ab <- col
		id.ef <- row
		id.c <- 0;id.d <- 0;id.cd <- 2*id.c+id.d+1
		#cat(id,"-",id.ab,"-",id.cd,"\n")
		w.cd <- all.proba$proba.jointe[id.ef,id.cd,id.ab]
		if (is.na(w.cd)){
			min.cd <- min(all.proba$pj.ab.cd[id.cd,id.ab],all.proba$pj.ab.ef[id.ef,id.ab])
		}
		else{min.cd <- 1.1}
		id.c <- 0;id.d <- 1;id.cd <- 2*id.c+id.d+1
		w.cD <- all.proba$proba.jointe[id.ef,id.cd,id.ab]
		if (is.na(w.cD)){
			min.cD <- min(all.proba$pj.ab.cd[id.cd,id.ab],all.proba$pj.ab.ef[id.ef,id.ab])
		}
		else{min.cD <- 1.1}
		id.c <- 1;id.d <- 0;id.cd <- 2*id.c+id.d+1
		w.Cd <- all.proba$proba.jointe[id.ef,id.cd,id.ab]
		if (is.na(w.Cd)){
			min.Cd <- min(all.proba$pj.ab.cd[id.cd,id.ab],all.proba$pj.ab.ef[id.ef,id.ab])
		}
		else{min.Cd <- 1.1}
		id.c <- 1;id.d <- 1;id.cd <- 2*id.c+id.d+1
		w.CD <- all.proba$proba.jointe[id.ef,id.cd,id.ab]
		if (is.na(w.CD)){
			min.CD <- min(all.proba$pj.ab.cd[id.cd,id.ab],all.proba$pj.ab.ef[id.ef,id.ab])
		}
		else{min.CD <- 1.1}
		mins <- c(min.cd,min.cD,min.Cd,min.CD)
		#print(mins)
		n.mins <- sum(mins <=1)
		if (n.mins > 1){
			for (i in 1:(n.mins-1)){
				ww <- order(mins)[i]
				id.cd <- ww
				vec <- c(all.proba$pj.ab.cd[id.cd,id.ab],all.proba$pj.ab.ef[id.ef,id.ab])
				tmp.max <- min(vec)
				tmp.p <- runif(1,0,tmp.max)
				#cat("[",0,"-",tmp.max,"] - ",tmp.p,"\n")
				all.proba$proba.jointe[id.ef,id.cd,id.ab] <- tmp.p
				all.proba$pj.ab.cd[id.cd,id.ab] <- all.proba$pj.ab.cd[id.cd,id.ab]-tmp.p
				all.proba$pj.ab.ef[id.ef,id.ab] <- all.proba$pj.ab.ef[id.ef,id.ab]-tmp.p
			}
		}
		ww <- order(mins)[n.mins]
		id.cd <- ww
		tmp.p <- all.proba$pj.ab.ef[id.ef,id.ab]
#		cat("[",0,"-",tmp.p,"] - ",tmp.p,"\n")
		all.proba$proba.jointe[id.ef,id.cd,id.ab] <- tmp.p
		all.proba$pj.ab.cd[id.cd,id.ab] <- all.proba$pj.ab.cd[id.cd,id.ab]-tmp.p
		all.proba$pj.ab.ef[id.ef,id.ab] <- all.proba$pj.ab.ef[id.ef,id.ab]-tmp.p
	}
	return(all.proba)
}


make.data.triplet <- function(pj.ab.cd.ef,n.obs=1000){
	#mat <- array(rmultinom(1,n.obs,all.proba$proba.jointe),dim=c(4,4,4))
  mat <- array(rmultinom(1,n.obs,pj.ab.cd.ef),dim=c(4,4,4))
  
	geno.R1 <- matrix(NA,ncol=3,nrow=n.obs)
	geno.R2 <- matrix(NA,ncol=3,nrow=n.obs)
	compt <- 1
	#print(mat)
	for (a in 1:2){
		for (b in 1:2){
			for (c in 1:2){
					for (d in 1:2){
						for (e in 1:2){
							for (f in 1:2){
								counts <- mat[2*(e-1)+f,2*(c-1)+d,2*(a-1)+b]
								if (counts > 0){
									geno.R1[(compt:(compt+counts-1)),1] <- a-1
									geno.R2[(compt:(compt+counts-1)),1] <- b-1
									geno.R1[(compt:(compt+counts-1)),2] <- c-1
									geno.R2[(compt:(compt+counts-1)),2] <- d-1
									geno.R1[(compt:(compt+counts-1)),3] <- e-1
									geno.R2[(compt:(compt+counts-1)),3] <- f-1
									compt <- compt+counts
								}
							}
						}
					}	
				}	
			}
		}
	return(list(geno.R1=geno.R1,geno.R2=geno.R2))
}


add.pair <- function(data,pj.ab.cd.ef,id.c.R1=2,id.d.R2=2,id.e.R1=3,id.f.R2=3){
  n.obs <- nrow(data$geno.R1)
  geno.R1 <- data$geno.R1[,c(id.c.R1,id.e.R1)]
  geno.R2 <- data$geno.R2[,c(id.d.R2,id.f.R2)]
  snp.X <- rep(NA,times=n.obs)
  snp.Y <- rep(NA,times=n.obs)
  for (i in 1:n.obs){
    proba <- pj.ab.cd.ef[(2*geno.R1[i,2]+geno.R2[i,2]+1),(2*geno.R1[i,1]+geno.R2[i,1]+1),]
    new  <- which(rmultinom(1,1,proba)==1)
    id.a <- as.numeric(new >=3)
    id.b <- (new+1)%%2
    snp.X[i] <- id.a
    snp.Y[i] <- id.b
  }
  data$geno.R1 <- cbind(data$geno.R1,snp.X)
  data$geno.R2 <- cbind(data$geno.R2,snp.Y)
  return(data)
}


