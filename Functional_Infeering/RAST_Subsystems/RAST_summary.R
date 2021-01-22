library(ggplot2)
library(plotly)
library(readxl)
library(corrplot)
library(reshape2)
library(plotly)
library(dplyr)
library(RColorBrewer)
library(psych)
library(ggdendro)

### Setwd
data <- read_excel("Rast_Cat_Comparison.xlsx")
Subcategory <- read_excel("Rast_Subcategory_Comparison.xlsx")
ANI <- read.csv("sourmash_kmer_51_mash_similarity.txt",sep = "\t", row.names = 1)
ANI <- as.matrix(ANI)

gut_keep <- c("Mycoplasma_pneumoniae","Mycoplasma_pirum",
              "Mycoplasma_mykiss","Mycoplasma_salar",
              "Mycoplasma_lavaretus","Mycoplasma_iowae",
              "Mycoplasma_fermentans", "Mycoplasma_arthritidis_158L3_1",
              "Mycoplasma_arginini","Mycoplasma_cloacale")

ANI <- ANI[gut_keep,gut_keep]
# Plot ANI to get order
order.plot <- corrplot(ANI,title = "Correlation Plot", 
                  method = "color", outline = T, 
                  addgrid.col = "darkgray", 
                  order="hclust",
                  cl.pos = "b", tl.col = "indianred4", 
                  tl.cex = 1.5, cl.cex = 1.5)

order <- row.names(order.plot)
order.1 <- order

## Original order with all genomes included
#order[1] <- "Mycoplasma_wenyonii"
#order[9] <- "Mycoplasma_iowae_695"
#order[11] <- "Mycoplasma_haemofelis_Langford1"
#order[14] <- "Mycoplasma_mobile_163K"
#order[15] <- "Mycoplasma_genitalium_G37"
#order[23] <- "Mycoplasma_capricolum"
#order[27] <- "Mycoplasma_cynos"
#order[34] <- "Mycoplasma_conjunctivae"
#order[36] <- "Mycoplasma_bovis_PG45"
#order[38] <- "Mycoplasma_arthritidis"
#order <- order[-35]
#order <- order[-10]
#order <- order[-3]

order[5] <- "Mycoplasma_iowae_695"
order[1] <- "Mycoplasma_arthritidis"
ANI_dist <- ANI
#ANI_dist <- ANI_dist[-c(35,21,1),-c(35,21,1)]

dist <- dist(ANI_dist, method = "maximum")^2
hc <- hclust(dist, method =  "complete")
dendro <- ggdendrogram(hc, 
                       rotate = FALSE, 
                       size = 10,
                       theme_dendro = TRUE)

pdf("ANI_maximum_dendrogram_small.pdf", width = 10, height = 5)
dendro  + scale_y_reverse()
dev.off()

## create color palette
colors <- c(brewer.pal(n = 8, name = 'Dark2'),
            brewer.pal(n = 8, name = 'Set2'),
            brewer.pal(n = 8, name = 'Paired'),
            c("white"))


colnames(data) <- c("No.", "Category", "Genome", "Subcategory", "Subsystem", "Role","Features","Total","Prob_Subs")

## Sort genomes for only gut related genomes
plot <- ggplot(data, aes(fill=Category, y=Prob_Subs, x=Genome)) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic() +scale_fill_manual(values=colors) +
  theme(axis.text.x = element_text(angle = 90))


summary <- 
  data %>%
  group_by(Category,Genome) %>%
  summarise(
    per=paste0(round(100*sum(Prob_Subs),2),'%'),
    Percent=round(100*sum(Prob_Subs),2))


test <- subset(Subcategory, New == "Mycoplasma_arthritidis" | New == "Mycoplasma_arginini" | New == "Mycoplasma_cloacale" | New == "Mycoplasma_fermentans" | New == "Mycoplasma_iowae_695" | New == "Mycoplasma_lavaretus" | New ==  "Mycoplasma_mykiss" | New == "Mycoplasma_salar" | New == "Mycoplasma_pneumoniae" | New == "Mycoplasma_pirum")
test_2 <- test[order(test$Category),]
y_order <- c(unique(test_2$Subsystem))

pdf("Subcategory_wo_legend.pdf", width = 10, height = 20)
ggplot(data = test_2, aes(x = New, y =factor(Subsystem, levels=c(y_order)))) +
  stat_bin_2d(aes(fill = Category), 
              binwidth = c(5,15),
              colour = 'white',
              size = 1.05)  +
  theme_classic()  +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "RAST Overview of Mycoplasma genomes") +
  scale_x_discrete(limits=order) + theme(legend.position = "none")
dev.off()
