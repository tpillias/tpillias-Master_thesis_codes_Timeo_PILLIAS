##### QTL mapping ####
# using multifam map and analysing family 106 only
rm(list = ls())
library(qtl)
library(dplyr)
library(LinkageMapView)

#### IMPORT DATA ####

# 1. make sure that everything runs properly. we'll test the pipeline on one chromosome only
# f.ex we'll use chromosome 2.
# we prepare a folder with:
# the genotype file for chromosome 2 in family 106 (map2genotypes.awk and manual tweaking)
# the phenotype file for family 106
# make sure that the phenotype and the genotype file have the same individuals with matching IDs
# the mapfile with marker position information
# make sure .csv files are COMMA (,) separated and not SEMICOLON (;)
# we use the qtl function read.cross
# with cross.type = "" corresponding the 106 cross type
#?read.cross

setwd("~/Desktop/Alexis_qtl/qtl_111")
setwd("C:/Users/Admin/Desktop/qtl_111")
f111 <- read.cross(format = "csvs", dir = ".", 
                   genfile = "genofile_denovo_preparation.csv", 
                   mapfile = "female_map_111.csv", 
                   phefile = "phenofile_111.csv",
                   genotypes = c("A","H","B"),
                   sep = ";", dec = ",")

# jittermap due to multiple markers mapping to the same position
f111 <- jittermap(f111)

# check
summary(f111)
#plotMap(f101)

#outfile = file.path("LG_viz_01.pdf")
# dupnbr to manage different markers in the same position
#lmv.linkage.plot(f101,outfile,dupnbr = TRUE, ruler = TRUE, denmap=TRUE)


#### DATA QUALITY CHECK ####

# check segregation distortion
# to see whether genotypes appear in the expected proportions
gt_111 <- geno.table(f111)
seg.dist_111<-gt_111[ gt_111$P.value < 1e-7, ]
write.csv(seg.dist_111, "seg.dist.female.csv")
# The last column is a p-value for a χ2 test of Mendelian proportions (1:2:1 in an intercross)

# Compare individual's genotypes
# compare the genotype data for each pair of individuals from a cross, to identify pairs that have unusually similar genotypes. 
# These may indicate sample mix-ups of some kind.
# We can first look at proportion of markers with identical genotypes for each pair of individuals within this family.
cg_111 <- comparegeno(f111)
hist(cg_111, breaks=200, xlab="Proportion of identical genotypes")
which(cg_111 > 0.9, arr.ind=TRUE)

# Pairwise recombination fractions
# check H19_103_rqtl, makes sense when there are more chromosomes
f111 <- est.rf(f111)
checkAlleles(f111)                   

# assessment of marker placement.
# plot.rf is used to plot the pairwise recombination fractions and LOD scores
plotRF(f111, alternate.chrid=TRUE)
# Purple corresponds to a large LOD score or a small recombination fraction, while yellow is the reverse. 
# Missing values would appear in light gray.

# check the order of markers on a chromosome using the *ripple function*. 
# check Rmd - takes a very long time to run even if it only counts number of crossovers (i.e. no maximum likelihood method)

# Identifying genotyping errors
f111 <- calc.errorlod(f111)
top <- top.errorlod(f111, cutoff=5)
top

# Counting number of crossovers(XOs)
nxo <- countXO(f111)
which(countXO (f111) == 0)
plot(nxo [1:82], ylab="No. crossovers")
mean(nxo[1:82])

# these are all inds for which we have phenotype and no genotype - all good

# Missing genotype information
plotInfo(f111, col="dodgerblue")

out.hk_111 <- scanone(f111, method="hk", model="binary")
plot(out.hk_111, ylab="LOD score", col = "dodgerblue")

out.all_hk_111<- scanone(f111, pheno.col=2:7, verbose = F, method = "hk")
plot(out.all_hk_111, ylab="LOD score", main = "method = HK")
plot(out.all_hk_111[, c("chr", "pos", "stripes_relaxed")],
     ylab = "LOD score",
     main = "QTL mapping for selected traits",
     col= "orange")

lod_sub <- out.all_hk_111 %>%
  filter(chr %in% c("18", "22"))

library(dplyr)
chr_to_plot <- "22"

lod_chr <- out.all_hk_111 %>% 
  as.data.frame() %>%
  filter(chr == chr_to_plot)
marker_pos <- lod_chr$pos

ggplot(lod_chr) +
  geom_line(aes(x = pos, y = stripes_strict), color = "darkorange", size = 0.6) +
  geom_line(aes(x = pos, y = stripes_relaxed), color = "steelblue", size = 0.6)+
  geom_rug(sides = "b", color = "black", alpha = 0.7, length = unit(0.02, "npc")) +
  theme_minimal(base_size = 14) +
  labs(title = "",
       x = "Position (cM)",
       y = "LOD score") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(), 
    axis.ticks = element_line(color = "black")
  )



summary(out.all_hk_111, threshold=3, format="allpheno")
summary(out.all_hk_111, threshold=3, format="allpeaks")

operm_all_hk_111 <- scanone(f111, pheno.col=2:7, n.perm = 1000, verbose = F, method = "hk")
summary(operm_all_hk_111, alpha=c(0.20, 0.1, 0.05), pvalues=TRUE)

hk_alpha020_111<-summary(out.all_hk_111, format="allpeaks", perms=operm_all_hk_111, alpha=0.20, pvalues=TRUE)
hk_alpha020_pheno111<-summary(out.all_hk_111, format="allpheno", perms=operm_all_hk_111, alpha=0.20, pvalues=TRUE)
hk_alpha010_111<-summary(out.all_hk_111, format="allpeaks", perms=operm_all_hk_111, alpha=0.10, pvalues=TRUE)
hk_alpha005_111<-summary(out.all_hk_111, format="allpeaks", perms=operm_all_hk_111, alpha=0.05, pvalues=TRUE)

write.csv(hk_alpha020_111, "./hk_alpha020_111_denovo.tsv") #PC8 (chr36:19.7), V6 (chr3:141.3, chr30:55.17), V1 (chr9:124.92) and PC4(chr37:9.28) 
write.csv(hk_alpha010_111, "./hk_alpha010_111_denovo.tsv") #PC8 (chr36:19.7), V6 (chr3:141.3) and V1 (chr9:124.92)
write.csv(hk_alpha005_111, "./hk_alpha005_111_denovo.tsv") #PC8 (chr36:19.7) and V6 (chr3:141.3)

plot(out.all_hk_111, lodcolumn=1:5, col=rainbow(5))

##### QTL LOCATION and QTL EFFECT #####
# using results from HK, in order of larger to lower significance
# sign codes (***) = 0.05; (**) = 0.1; (*) = 0.2


## BAYES INTERVAL AND BOOTSTRAPPING
bayesint(out.all_hk_111, c(18), 0.95, expandtomarkers=TRUE) #method 1
out.boot <- scanoneboot(f111, chr=18, n.boot=1000) # method 2
summary(out.boot)



set.seed(123)
out.boot <- scanoneboot(f111, chr = 22, n.boot = 1000)
boot_pos <- summary(out.boot)$pos
quantile(boot_pos, c(0.025, 0.975))


##### QTL EFFECT #####

# for stripes strict
find.marker(f111, 18, 25.3) #args: family file, chromosome and position, you'll get the name of a marker, e.g. "3307".
effectplot(f111, mname1="433", pheno.col = 4) # args: family file, marker name obtained in previous function, pheno.col corresponding to stripes strict.
eff <- effectplot(f111, mname1="433", draw=FALSE, pheno.col = 4) # same as above
eff
plotPXG(f111, marker="433", pheno.col = 4)
