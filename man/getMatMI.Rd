\name{getMatMI}
\alias{getMatMI}
\title{Computation of the mutual information between SNP pairs in genomic regions.}
\usage{
getMatMI(Region1,Region2=NULL,allelic=FALSE)
}
\arguments{
\item{Region1}{A nxp1 \code{\link{matrix}} or \code{\link{data.frame}} containing allele or genotype of the p1 SNPs for the n individuals.}
\item{Region2}{A nxp2 \code{\link{matrix}} or \code{\link{data.frame}} containing allele or genotype of the p2 SNPs for the n individuals. Region 1 and Region 2 must have the same number of rows.}
\item{allelic}{A flag indicating whether data are allelic \code{allelic=TRUE} or genotypic \code{allelic=FALSE}.}
}

\value{
	An object of class \code{MIMatrix} which is a list containing the following components:
	\item{mat.MI}{A a p1xp2xp1xp2 \code{\link{array}} that stores the mutual information between all pairs of SNPs}
	\item{nb.locus}{An \code{\link{integer}} with two possible value: \code{nb.locus=1} corresponding to the computation of the mutual information within 1 genomic region (i.e. when \code{Region2=NULL}) and \code{nb.locus=2} for the computation of the mutual information between 2 genomic regions (i.e. when \code{Region2!=NULL}).}
}
\description{
	\code{getMatMI} allows the computation of the mutual information between SNP pairs within 1 genomic regions or between 2 genomic regions.
	}
\details{
	\code{getMatMI} is an implementation of the statistical procedure proposed in Friguet and  Emily (2013). When considering only 1 genomic region (i.e. when \code{Region2=NULL}), \code{getMatMI} computes the pairwise mutual information between all pairs of SNPs : MI((i,j),(k,l)) where MI is the Mutual Information, (i,j) two SNPs within the genomic region and (k,l) two SNPs within the same genomic regions.  When considering only 2 genomic regions (i.e. when \code{Region2!=NULL}), \code{getMatMI} computes the pairwise mutual information between pairs of SNPs across the two regions: MI((i,j),(k,l)) where MI is the Mutual Information, (i,j) two SNPs within the first genomic region and (k,l) two SNPs within the second genomic regions.
}
\examples{
############
## Example with simulated allelic data
############
data(Region1)
data(Region2)
MatMI <- getMatMI(Region1=RegionA,Region2=RegionB,allelic=TRUE)

############
## Example with a real subset of genotypic data
############
## Loading of the genotypic data
load(system.file("extdata/Chr15_20S.Rdata", package="EpiTag"))
## Computation of the mutual information
MatMI15 <- getMatMI(Region1=Chr15.20S)

}
\references{
Friguet, C. and Emily, M. Sélection de marqueurs pour la détection d'interactions de gènes. Proceedings of 45èmes Journées de Statistique, May 2013, Toulouse, France. pp.178, 2013.\cr
}

\seealso{
\code{\link{EpiTag}}
}
