
<link href="http://joey711.github.com/phyloseq/markdown.css" rel="stylesheet"></link>

# Simulations Design, Tiny Demo Code

The following two examples are the exact code used to create the values used in the simplified overview of the two simulation types in the article. As you can see, even these simplified examples are based on real data, and are fully reproducible.

---

## Sample Clustering Simulations Example

### Define templates from real microbiomes


```r
set.seed(20140206)
library("phyloseq")
data("GlobalPatterns")
OceanReal = subset_samples(GlobalPatterns, SampleType == "Ocean")
OceanTop = names(sort(taxa_sums(OceanReal), decreasing = TRUE)[1:3])
FecesReal = subset_samples(GlobalPatterns, SampleType == "Feces")
FecesTop = names(sort(taxa_sums(FecesReal), decreasing = TRUE)[1:3])
# Get matrix with both
BothOF = subset_samples(GlobalPatterns, SampleType %in% c("Ocean", "Feces"))
BothOF = prune_taxa(c(FecesTop, OceanTop), BothOF)
BothOF = as(otu_table(BothOF), "matrix")[, c(sample_names(OceanReal), sample_names(FecesReal))]
BothOF
```

```
##          NP2    NP3    NP5 M31Fcsw M11Fcsw  TS28   TS29
## 12812  15527  15682 161687       7      13    28     14
## 557211 87667   4013  72549      14      23    30     29
## 534609 10712 148400  15722       9      10    17     16
## 158660    13     39     28   82686  244551  7978  24572
## 331820    24    101    105  354695  452219 92915   1180
## 189047     4     33     29   14353    9419 33789 251215
```

```r
# Divide by 1000, for example simplicity
BothOF = floor(BothOF/1000)
BothOF
```

```
##        NP2 NP3 NP5 M31Fcsw M11Fcsw TS28 TS29
## 12812   15  15 161       0       0    0    0
## 557211  87   4  72       0       0    0    0
## 534609  10 148  15       0       0    0    0
## 158660   0   0   0      82     244    7   24
## 331820   0   0   0     354     452   92    1
## 189047   0   0   0      14       9   33  251
```

```r
OceanVec = rowSums(BothOF[, sample_names(OceanReal)])
OceanVec = matrix(OceanVec, nrow = length(OceanVec), dimnames = list(names(OceanVec)))
OceanVec
```

```
##        [,1]
## 12812   191
## 557211  163
## 534609  173
## 158660    0
## 331820    0
## 189047    0
```

```r
FecesVec = rowSums(BothOF[, sample_names(FecesReal)])
FecesVec = matrix(FecesVec, nrow = length(FecesVec), dimnames = list(names(FecesVec)))
FecesVec
```

```
##        [,1]
## 12812     0
## 557211    0
## 534609    0
## 158660  357
## 331820  899
## 189047  307
```

```r
# Write the example table and the example multinomial
write.csv(BothOF, file = "sim_cluster_ocean_feces_otu_matrix.csv", col.names = FALSE, 
    row.names = FALSE)
write.csv(OceanVec, file = "ocean-multinomial.csv", col.names = FALSE, row.names = FALSE)
write.csv(FecesVec, file = "feces-multinomial.csv", col.names = FALSE, row.names = FALSE)
```


### Deterministic mixing

Mix each multinomial by adding total/EffectSize counts from the other. Create mixed multinomial by adding counts from the other in precise proportion, a total of `Library Size / Effect Size`.


```r
EffectSize = 10
addToOcean = round(sum(OceanVec) * FecesVec/(sum(FecesVec) * EffectSize), 0)
addToFeces = round(sum(FecesVec) * OceanVec/(sum(OceanVec) * EffectSize), 0)
# Add them together to create 'dirty' multinomial
dirtyOcean = addToOcean + OceanVec
dirtyOcean
```

```
##        [,1]
## 12812   191
## 557211  163
## 534609  173
## 158660   12
## 331820   30
## 189047   10
```

```r
dirtyFeces = addToFeces + FecesVec
dirtyFeces
```

```
##        [,1]
## 12812    57
## 557211   48
## 534609   51
## 158660  357
## 331820  899
## 189047  307
```

```r
# Write the dirty multinomials
write.csv(dirtyOcean, file = "dirty-ocean-multinomial.csv", col.names = FALSE, 
    row.names = FALSE)
write.csv(dirtyFeces, file = "dirty-feces-multinomial.csv", col.names = FALSE, 
    row.names = FALSE)
# Example of 'simulated' count matrix with 5 columns/samples each class
J = 5
NLOcean = sample(colSums(BothOF), size = J, replace = TRUE)
NLFeces = sample(colSums(BothOF), size = J, replace = TRUE)
# Simulate Ocean
OceanSim = sapply(NLOcean, function(NL, dirtyOcean) {
    table(sample(rownames(dirtyOcean), NL, replace = TRUE, prob = dirtyOcean))[rownames(dirtyOcean)]
}, dirtyOcean)
# Simulate Feces
FecesSim = sapply(NLFeces, function(NL, dirtyFeces) {
    table(sample(rownames(dirtyFeces), NL, replace = TRUE, prob = dirtyFeces))[rownames(dirtyFeces)]
}, dirtyFeces)
# Convert NA to zero
OceanSim[is.na(OceanSim)] <- 0L
FecesSim[is.na(FecesSim)] <- 0L
# Write simulated table
write.csv(cbind(OceanSim, FecesSim), file = "cluster_ocean_feces-ex-sim.csv", 
    col.names = FALSE, row.names = FALSE)
```



---

## Differential Abundance Simulations Example


```r
# Reset package and data call, for modularity and protection
library("phyloseq")
data("GlobalPatterns")
SkinReal = subset_samples(GlobalPatterns, SampleType == "Skin")
SkinTop = names(sort(taxa_sums(SkinReal), decreasing = TRUE)[1:6])
SkinMat = round(as(otu_table(prune_taxa(SkinTop, SkinReal)), "matrix")/1000, 
    0)
SkinMat
```

```
##        M31Plmr M11Plmr F21Plmr
## 64396       34       1      15
## 589787       4      20       4
## 94166       29       1       6
## 484436       1      85       3
## 98605      161       6      13
## 332405      42       2       3
```

```r
# Write to csv
write.csv(SkinMat, file = "sim_diff_abund_matrix_ex_real_matrix.csv", col.names = FALSE, 
    row.names = FALSE)
# Create the rowsum version (multinomial)
SkinMatRS = matrix(rowSums(SkinMat), nrow = nrow(SkinMat), dimnames = list(rownames(SkinMat)))
SkinMatRS
```

```
##        [,1]
## 64396    50
## 589787   28
## 94166    36
## 484436   89
## 98605   180
## 332405   47
```

```r
write.csv(SkinMatRS, file = "sim_diff_abund_matrix_ex_real_matrix_RS.csv", col.names = FALSE, 
    row.names = FALSE)
# Simulate by sampling from multinomial Example of 'simulated' count matrix
# with 4 columns/samples each class
J = 4
NL = sample(round(sample_sums(GlobalPatterns)/10000, 0), size = J * 2, replace = TRUE)
NL
```

```
##     NP3 M31Tong M11Fcsw M11Plmr   Even2  AQC4cm     SV1 M31Plmr 
##     148     200     208      43      97     236      70      72
```

```r
# Simulate Test and NULL
nullmat = sapply(NL, function(NL, SkinMatRS) {
    table(sample(rownames(SkinMatRS), NL, replace = TRUE, prob = SkinMatRS))[rownames(SkinMatRS)]
}, SkinMatRS)
nullmat[is.na(nullmat)] <- 0L
nullmat
```

```
##        NP3 M31Tong M11Fcsw M11Plmr Even2 AQC4cm SV1 M31Plmr
## 64396   14      21      29       7     8     17   9       6
## 589787  14      11       9       3     5     14   2       4
## 94166    8      19      17       5    10     19   7       7
## 484436  33      42      37      10    18     51  19      15
## 98605   60      86      94      14    48    100  26      32
## 332405  19      21      22       4     8     35   7       8
```

```r
write.csv(nullmat, file = "sim_diff_abund_matrix_before_effect.csv", col.names = FALSE, 
    row.names = FALSE)
# Apply 'effect' to random (arbitrary in this case) rows
testmat = nullmat
EffectSize = 10
effectrows = c(1, 4)
effectcols = 1:J
testmat[effectrows, effectcols] <- EffectSize * testmat[effectrows, effectcols]
testmat
```

```
##        NP3 M31Tong M11Fcsw M11Plmr Even2 AQC4cm SV1 M31Plmr
## 64396  140     210     290      70     8     17   9       6
## 589787  14      11       9       3     5     14   2       4
## 94166    8      19      17       5    10     19   7       7
## 484436 330     420     370     100    18     51  19      15
## 98605   60      86      94      14    48    100  26      32
## 332405  19      21      22       4     8     35   7       8
```

```r
# write to csv
write.csv(testmat, file = "sim_diff_abund_matrix_after_effect.csv", col.names = FALSE, 
    row.names = FALSE)
```

