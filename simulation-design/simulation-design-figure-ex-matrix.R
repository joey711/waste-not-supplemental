################################################################################
#
# Sample Clustering, Mixing
#
################################################################################
set.seed(20140206)
library("phyloseq")
data("GlobalPatterns")
OceanReal = subset_samples(GlobalPatterns, SampleType=="Ocean")
OceanTop = names(sort(taxa_sums(OceanReal), decreasing=TRUE)[1:3])
FecesReal = subset_samples(GlobalPatterns, SampleType=="Feces")
FecesTop = names(sort(taxa_sums(FecesReal), decreasing=TRUE)[1:3])
# Get matrix with both
BothOF = subset_samples(GlobalPatterns, SampleType %in% c("Ocean", "Feces"))
BothOF = prune_taxa(c(FecesTop, OceanTop), BothOF)
BothOF = as(otu_table(BothOF), "matrix")[, c(sample_names(OceanReal), sample_names(FecesReal))]
# Divide by 1000, for example simplicity
BothOF = floor(BothOF/1000)
OceanVec = rowSums(BothOF[, sample_names(OceanReal)])
OceanVec = matrix(OceanVec, nrow=length(OceanVec), dimnames=list(names(OceanVec)))
FecesVec = rowSums(BothOF[, sample_names(FecesReal)])
FecesVec = matrix(FecesVec, nrow=length(FecesVec), dimnames=list(names(FecesVec)))
# Write the example table and the example multinomial
write.csv(BothOF,   file="sim_cluster_ocean_feces_otu_matrix.csv", col.names=FALSE, row.names=FALSE)
write.csv(OceanVec, file="ocean-multinomial.csv", col.names=FALSE, row.names=FALSE)
write.csv(FecesVec, file="feces-multinomial.csv", col.names=FALSE, row.names=FALSE)
# Defined artificial mixing.
# Mix each multinomial by adding total/EffectSize counts from the other.
# Create mixed multinomial by adding counts from the other in precise proportion,
# a total of Library Size / Effect Size
EffectSize = 10
addToOcean = round(sum(OceanVec) * FecesVec / (sum(FecesVec) * EffectSize), 0)
addToFeces = round(sum(FecesVec) * OceanVec / (sum(OceanVec) * EffectSize), 0)
# Add them together to create "dirty" multinomial
dirtyOcean = addToOcean + OceanVec
dirtyFeces = addToFeces + FecesVec
# Write the dirty multinomials
write.csv(dirtyOcean, file="dirty-ocean-multinomial.csv", col.names=FALSE, row.names=FALSE)
write.csv(dirtyFeces, file="dirty-feces-multinomial.csv", col.names=FALSE, row.names=FALSE)
# Example of "simulated" count matrix with 5 columns/samples each class
J = 5
NLOcean = sample(colSums(BothOF), size=J, replace=TRUE)
NLFeces = sample(colSums(BothOF), size=J, replace=TRUE)
# Simulate Ocean
OceanSim = sapply(NLOcean, function(NL, dirtyOcean){
  table(sample(rownames(dirtyOcean), NL, replace=TRUE, prob=dirtyOcean))[rownames(dirtyOcean)]
}, dirtyOcean)
# Simulate Feces
FecesSim = sapply(NLFeces, function(NL, dirtyFeces){
  table(sample(rownames(dirtyFeces), NL, replace=TRUE, prob=dirtyFeces))[rownames(dirtyFeces)]
}, dirtyFeces)
# Convert NA to zero
OceanSim[is.na(OceanSim)] <- 0L
FecesSim[is.na(FecesSim)] <- 0L
# Write simulated table
write.csv(cbind(OceanSim, FecesSim), file="cluster_ocean_feces-ex-sim.csv", col.names=FALSE, row.names=FALSE)
################################################################################
#
# Differential Abundance
#
################################################################################
library("phyloseq")
data("GlobalPatterns")
SkinReal = subset_samples(GlobalPatterns, SampleType=="Skin")
SkinTop = names(sort(taxa_sums(SkinReal), decreasing=TRUE)[1:6])
SkinMat = round(as(otu_table(prune_taxa(SkinTop, SkinReal)), "matrix")/1000, 0)
# Write to csv
write.csv(SkinMat, file="sim_diff_abund_matrix_ex_real_matrix.csv", col.names=FALSE, row.names=FALSE)
# Create the rowsum version (multinomial)
SkinMatRS = matrix(rowSums(SkinMat), nrow=nrow(SkinMat), dimnames=list(rownames(SkinMat)))
write.csv(SkinMatRS, file="sim_diff_abund_matrix_ex_real_matrix_RS.csv", col.names=FALSE, row.names=FALSE)
# Simulate by sampling from multinomial
# Example of "simulated" count matrix with 4 columns/samples each class
J = 4
NL = sample(round(sample_sums(GlobalPatterns)/10000, 0), size=J*2, replace=TRUE)
# Simulate Test and NULL
nullmat = sapply(NL, function(NL, SkinMatRS){
  table(sample(rownames(SkinMatRS), NL, replace=TRUE, prob=SkinMatRS))[rownames(SkinMatRS)]
}, SkinMatRS)
nullmat[is.na(nullmat)] <- 0L
write.csv(nullmat, file="sim_diff_abund_matrix_before_effect.csv", col.names=FALSE, row.names=FALSE)
# Apply "effect" to random (arbitrary in this case) rows
testmat = nullmat
EffectSize = 10
effectrows = c(1, 4)
effectcols = 1:J
testmat[effectrows, effectcols] <- EffectSize * testmat[effectrows, effectcols]
# write to csv
write.csv(testmat,  file="sim_diff_abund_matrix_after_effect.csv", col.names=FALSE, row.names=FALSE)