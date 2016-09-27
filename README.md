# EpiTag
The aim of this package is to propose a method for selecting TagSNPs to optimize power for detecting epistasis.

## Installation
To install and load the package in R

```ruby
library(devtools)
install_github("MathieuEmily/EpiTag")
library(EpiTag)
```

Function EpiTag can take as input two data.frame Region1 and Region2, where Region1 has n rows (individuals) and p1 columns (number of SNPs) and Region2 has n rows (individuals) and p2 columns (number of SNPs). Cell (i,j) of the input matrices contains the allele or the genotype of the individual i for the SNP j.


## Example:

```ruby
data(RegionA)
data(RegionB)
```

By typing 
```ruby
res.EpiTag <- EpiTag(Region1=RegionA,Region2=RegionB,threshold=0.81,allelic=TRUE)
```
you can get the set of tag SNPs selected in each genomic region:
```ruby
> res.EpiTag$tagSNPs.R1
[1] "SNPA.4"
> res.EpiTag$tagSNPs.R2
[1] "SNPB.3"
```