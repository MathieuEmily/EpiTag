source("/Users/memily/Github/EpiTag/R/getMatMI.R")
source("/Users/memily/Github/EpiTag/R/EpiTag.R")
source("/Users/memily/Github/EpiTag/R/Functions4EpiTag.R")
source("/Users/memily/Github/EpiTag/R/Functions4Simulation.R")

##########
## SCRIPT
##########

n.obs <- 1000
n.sim.pj <- 200
pj.max <- NULL
pj.min <- NULL
my.max <- 0
my.min <- 1


resMaxMin.1 <- get.pj.MaxMin(p.A=0.6,p.B=0.6,p.C=0.6,p.D=0.6,p.E=0.6,p.F=0.6,r2.AB=0,r2.CD=0,r2.EF=0,r2.AC=0.9,r2.BD=0.9,r2.AE=0.9,r2.BF=0.9,n.sim.pj=n.sim.pj)

resMaxMin.2 <- get.pj.MaxMin(p.A=0.6,p.B=0.6,p.C=0.6,p.D=0.6,p.E=0.6,p.F=0.6,r2.AB=0,r2.CD=0,r2.EF=0,r2.AC=0.85,r2.BD=0.85,r2.AE=0.85,r2.BF=0.85,n.sim.pj=n.sim.pj)

### matrice de proba jointes des 3 premieres paires de SNPs
proba.jointe.3.1 <- generate.pj.triplet(pj.ab.cd=resMaxMin.1$pj.min,pj.ab.ef=resMaxMin.1$pj.max)$proba.jointe 

### matrice de proba jointes entre un "nouveau" SNP (le 4eme, le 5eme par exemple) 
# et le 2eme et le 3 eme SNPs
proba.jointe.3.2 <- generate.pj.triplet(pj.ab.cd=resMaxMin.2$pj.min,pj.ab.ef=resMaxMin.2$pj.list[[order(resMaxMin.2$MI)[100]]])$proba.jointe 


data <- make.data.triplet(pj.ab.cd.ef=proba.jointe.3.1,n.obs=n.obs) ## Simulation des 3 premieres paires SNPs
data <- add.pair(data=data,pj.ab.cd.ef=proba.jointe.3.2) ## Rajout d'une 4eme paire de SNP
data <- add.pair(data=data,pj.ab.cd.ef=proba.jointe.3.2) ## Rajout d'une 5eme paire de SNP


RegionA <- as.data.frame(data$geno.R1)
names(RegionA) <- paste("SNPA.",1:5,sep="")
RegionB <- as.data.frame(data$geno.R2)
names(RegionB) <- paste("SNPA.",1:5,sep="")
#package.skeleton("titi",list=c("RegionA","RegionB"),force=TRUE)

MatMI <- getMatMI(Region1=RegionA,Region2=RegionB,allelic=TRUE)
my.EpiTag <- EpiTag(Region1=RegionA,Region2=RegionA,threshold=0.81,allelic=TRUE)
my.EpiTag2 <- EpiTag(mat.MI=MatMI,threshold=0.81,allelic=TRUE)


tmpRA <- RegionA
tmpRB <- RegionB[,1:4]
tmpMatMI <- getMatMI(Region1=tmpRA,Region2=tmpRB,allelic=TRUE)
# MatMI <- getMatMI.2Regions(data1=data$geno.R1,data2=data$geno.R2,allelic=TRUE)
tmp.EpiTag <- EpiTag.2Regions(mat.MI=tmpMatMI$mat.MI,threshold=0.81)

# MatMI <- getMatMI.2Regions(data1=RegionA,data2=RegionB,allelic=TRUE)
# my.EpiTag <- EpiTag.2Regions(mat.MI=MatMI,threshold=0.90)


######
######

chr15 <- read.table("/Users/memily/RennesII/Professionel/Recherche/Seminaires-Conferences/2013/JDS2013/simulations/JDS13/chr15.dat")
dim(chr15)


MC <- cor(chr15[1:100,1:20])
cor.plot(MC)
test <- chr15[1:100,1:20]
dim(test)
Chr15.20S <- test

sort(sample(1:100,20))

names(Chr15.20S) <- paste("SNP",sort(sample(1:100,20)),sep="")
MatMI15 <- getMatMI(Region1=Chr15.20S)
res.EpiTag <- EpiTag(mat.MI=MatMI15,threshold=0.6,allelic=FALSE)
save(file="/Users/memily/Github/EpiTag/inst/extdata/MatMIChr15.Rdata",MatMI15)
save(file="/Users/memily/Github/EpiTag/inst/extdata/Chr15_20S.Rdata",Chr15.20S)
