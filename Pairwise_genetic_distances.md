# Explore pairwise genetic distances between individuals

#### Download R packages
```
library(ade4)
library(adegenet)
library(pegas)
library(ggplot2)
library(poppr)
library(reshape)
library(cowplot)
```
### 1. Load a VCF file and convert to GENIND object
```
setwd("file path")

vcf <- read.vcf("NAME.vcf", to = 130572)
X<- loci2genind(vcf)
```

### 2. Compute the distance matrix reflecting the percentage of allelic differences between two individuals
```
dist_matrice <- diss.dist(X, percent = TRUE, mat = TRUE)
pairwise_distance_ld [dist_matrice == 0] <- NA
pairwise_distance_ld <- melt(dist_matrice)
write.table(pairwise_distance,file="./distmatrice_pairs.txt", sep = "\t", quote=FALSE, row.names = TRUE, col.names = TRUE)
```
### 3. Visualize the frequency of pairwise genetic distances from distances matrices of both species

* Distance matrix 'distmatrice_pairs_lr.txt' and 'distmatrice_pairs_ld.txt' are available in the repository.
* Clonal assessment: we chose the maximum value of technical replicates (d = 0.022) to separate disctinct Multi-locus lineages (MLLs) 
* Technical replicates were excluded from the following analyses but pairwise distances of replicates are reported below

##### *Laminaria rodriguezii (n=46)*

Replicate 1 | Replicate 2 | distance (d)
------------|------------|------------
LRBM2P27_ir | LRBM2P27_r | 0.022
LRBO2_ir | LRBO2_r | 0.014

```
setwd("file path")
distmatrice_pairs_lr <- read.table("distmatrice_pairs_lr.txt", header = TRUE)

distance_lr <- ggplot(distmatrice_pairs_lr, aes(x=value)) +  
  geom_histogram(color="black", fill="firebrick1", alpha =0.6, bins  = 100 , position = "dodge",aes(x=value,y=(..count..)/length(distmatrice_pairs_lr$value))) +
  scale_x_continuous(breaks = seq(0,0.3,0.01), name = '% of difference',limits = c(0,0.30))+
  scale_y_continuous(name = 'Frequency of pairwise comparaison',labels = scales::number_format(accuracy = 0.01))+
  geom_vline(xintercept=0.022, linetype="dashed", color = "black",size=1)+
  geom_text(aes(x=0.035, label="d = 0.022", y=0.05), colour="black", angle=0, vjust = 1.2, text=element_text(size=11))+
  theme_bw() +
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        text=element_text(family="Tahoma"),
        axis.title = element_text(size = 11),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10)) 
```
##### *Laminaria digitata (n=116)*

Replicate 1 | Replicate 2 | distance (d)
------------|------------|------------
LD_ALL_2018_01_ir | LD_ALL_2018_01_r | 0.020
LD_BODO_2018_8_ir | LD_BODO_2018_8_r | 0.0009
LD_CLA_2018_01_ir | LD_CLA_2018_01_r | 0.001
LD_PER_2018_01_ir | LD_PER_2018_01_r | 0.0004
LD_SBOB_2018_01_ir | LD_SBOB_2018_01_r | 0.0009

```
setwd("file path")
distmatrice_pairs_ld <- read.table("distmatrice_pairs_ld.txt", header = TRUE)

distance_ld <- ggplot(distmatrice_pairs_ld, aes(x=value)) +  
  geom_histogram(color="black", fill="aquamarine3", alpha =0.6, bins  = 100 , position = "dodge",aes(x=value,y=(..count..)/length(distmatrice_pairs_ld$value))) +
  scale_x_continuous(breaks = seq(0,0.3,0.01), name = '% of difference',limits = c(0,0.30))+
   scale_y_continuous(name = 'Frequency of pairwise comparaison',labels = scales::number_format(accuracy = 0.01),limits= c(0,0.05))+
  geom_vline(xintercept=0.022, linetype="dashed", color = "black",size=1)+
  geom_text(aes(x=0.035, label="d = 0.022", y=0.05), colour="black", angle=0, vjust = 1.2, text=element_text(size=11))+
  theme_bw() +
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        text=element_text(family="Tahoma"),
        axis.title = element_text(size = 11),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10)) 
```
### 4. Produce the figures as they appear in the manuscript and save in image format (300 dpi)
```
plot_grid(distance_lr,distance_ld,
          labels =c("A","B"),label_fontface = "bold", ncol=1, nrow = 2,
          label_size = 17,
          label_x = 0.12, label_y = 0.9,
          hjust = 1, vjust = 1,
          label_fontfamily = "Tahoma")
          
ggsave("Figure1.jpg",width=40 ,height=15,dpi=300,units="cm")
```
          



