setwd("~/Desktop")

library(readxl)
library(ggplot2)
library(dplyr) 
library(micropan)
library(ggdendro)
library(phyloseq)
library(ranacapa)
library(igraph)
library(knitr)
library(vegan)

data <- t(read.table("MYCOPLASMA-GCs-occurrence-frequency.txt"))
h.est <- heaps(data, n.perm = 1000)
# If alpha < 1 it indicates an open pan-genome, according to https://www.sciencedirect.com/science/article/pii/S1369527408001239?via%3Dihub#fig1.
summary(h.est)

print(h.est)
print(chao(data))
# 35954 total number of gene clusters for Mycoplasma. When based on present genomes, this is a rather uncertain estimate!
fitted <- binomixEstimate(data, K.range = 2:40)
print(fitted$BIC.tbl) # Lowest Bayesian information criterion (BIC) estimates best fit to model. 
fitted$BIC.tbl$K.range[fitted$BIC.tbl$BIC==min(fitted$BIC.tbl$BIC)]
# K 5 have the lowest BIC, indicating best fit.
# We also see that in this row the estimate of pan-genome size is 38635, and the size of the core-genome is 10.
plot(fitted$BIC.tbl$BIC~fitted$BIC.tbl$K.range,
     col="dark blue", lwd = 2, 
     xlab="K value", ylab="Fitted model BIC", 
     main="BICs across range of K for fitted model")
axis(side=1,at=seq(0,49,5))

sp <- specaccum(data, method = "random", permutations=1000)
summary(sp)
#pdf("GC_Accumulation.pdf", height = 6, width = 8)
plot(sp, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col="grey80", xlab="Genomes", ylab="Gene Clusters", main="Gene Cluster accumulation plot")
boxplot(sp, col="red", add=TRUE)
## Include information of heaps' law alpha and gamma to the plot, based on h.est.
text(32,2000,expression(paste(gamma, " = 1- ",alpha," = 0.719, where ",alpha," = 0.281", sep ="\n")), cex = 0.7)
#dev.off()
save.image(file = "Heaps_Summary_Mycoplasma_PAN.RData")
