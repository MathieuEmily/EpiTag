\name{EpiTag}
\alias{EpiTag}
\title{Selection of TagSNPs }
\usage{
EpiTag(Region1=NULL,Region2=NULL,mat.MI=NULL,threshold=0.8,allelic=FALSE)
}
\arguments{
\item{Region1}{A nxp1 \code{\link{matrix}} or \code{\link{data.frame}} containing allele or genotype of the p1 SNPs for the n individuals.}
\item{Region2}{A nxp2 \code{\link{matrix}} or \code{\link{data.frame}} containing allele or genotype of the p2 SNPs for the n individuals. Region 1 and Region 2 must have the same number of rows.}
\item{mat.MI}{An object of the class \code{MIMatrix} where the \code{mat.MI} element is that contains a p1xp2xp1xp2 \code{\link{array}}. For example, an output of the \code{\link{getMatMI}} function.}
\item{threshold}{Threshold applied to the Mutual Information to create the successive bins of SNP pairs.}
\item{allelic}{A flag indicating whether data are allelic \code{allelic=TRUE} or genotypic \code{allelic=FALSE}.}
}

\value{
A list containing the following components:
	\item{tagSNPs.R1}{A \code{\link{vector}} that stores tagSNPs selected in \code{Region1}.}
	\item{tagSNPs.R2}{A \code{\link{vector}} that stores tagSNPs selected in \code{Region2}. Returned only if \code{Region2} is provided or if \code{mat.MI} contains the Mutual Information between two regions.}
	\item{pairs.tag}{A list of pairs of SNPs selected as a tag pair.}
	\item{pairs.bins}{A list of bins of SNPs pairs. Each bin is ix2 matrix where each row is a SNP pair tagged by the corresponding tag pair in the object \code{pairs.tag} and i is the number of such SNP pairs.}
}
\description{
	\code{EpiTag} allows the selection of tag SNPs in 1 or 2 genomic regions. The selection aims at optimizing the power to detect epistasis between SNPs.
	}
\details{
	\code{SMILE} is an implementation of the statistical procedure proposed in Friguet and  Emily (2013). The EpiTag method is based on the computation of the "matrix" of mutual information between all pairs of SNP pairs within a single genomic region and between 2 genomic regions. Mutual information can be computed before the selection, by using the argument \code{mat.MI} or during the selection process if \code{mat.MI=NULL}. Based on the mutual information, EpiTag performs a greedy algorithm to iteratively select the subset of most informative SNPs to optimize power to detect epistasis.
}
\examples{
############
## Example with simulated allelic data
############
data(Region1)
data(Region2)

## Computation of the mutual information and of the tag SNP selection simultaneously
res.EpiTag1 <- EpiTag(Region1=RegionA,Region2=RegionA,threshold=0.81,allelic=TRUE)

## Same results with a two steps computation of the mutual information and of the tag SNP selection
MatMI <- getMatMI(Region1=RegionA,Region2=RegionB,allelic=TRUE)
res.EpiTag2 <- EpiTag(mat.MI=MatMI,threshold=0.81,allelic=TRUE)


############
## Example with a real subset of genotypic data
############
## Loading of the genotypic data
load(system.file("extdata/Chr15_20S.Rdata", package="EpiTag"))
## Loading of the precomputed array of Mutual Information
load(system.file("extdata/MatMIChr15.Rdata", package="EpiTag"))

## Computation of the mutual information and of the tag SNP selection simultaneously
MatMI15 <- getMatMI(Region1=Chr15.20S)
res.EpiTag.real1 <- EpiTag(mat.MI=MatMI15,threshold=0.6,allelic=FALSE)

## Same results with the  precomputation of the mutual information
res.EpiTag.real2 <- EpiTag(mat.MI=MatMI15,threshold=0.6,allelic=FALSE)

}
\references{
Friguet, C. and Emily, M. Sélection de marqueurs pour la détection d'interactions de gènes. Proceedings of 45èmes Journées de Statistique, May 2013, Toulouse, France. pp.178, 2013.\cr
}

\seealso{
\code{\link{getMatMI}}
}
