# EpiTag

To install and load the package in R

```ruby
library(devtools)
install_github("MathieuEmily/EpiTag")
library(EpiTag)
```

Function EpiTag can take as input two data.frame Region1 and Region2, where Region1 has n rows (individuals) and p1 columns (number of SNPs) and Region2 has n rows (individuals) and p2 columns (number of SNPs). Cell (i,j) of the input matrices contains the allele or the genotype of the individual i for the SNP j.

data(Region1)
data(Region2)

## Computation of the mutual information and of the tag SNP selection simultaneously
res.EpiTag1 <- EpiTag(Region1=RegionA,Region2=RegionA,threshold=0.81,allelic=TRUE)



# Example:

```ruby
data(Region1)
data(Region2)
```

By typing 
```ruby
res.EpiTag1 <- EpiTag(Region1=RegionA,Region2=RegionA,threshold=0.81,allelic=TRUE)
```
you can get the graphical display of the clustering and the following responses:
```ruby
$selected
[1] "ShP"

$hc

Call:
hclust(d = as.dist(MatDis), method = "single")

Cluster method   : single 
Number of objects: 30 


$d2s
           BMD        BeT        Bgl        BoC        BoT        BrS        CoS        Dac
BeT 0.58281942                                                                             
Bgl 0.52995976 0.63365377                                                                  
BoC 0.33822047 0.25937500 0.33874194                                                       
BoT 0.53028294 0.38070451
```