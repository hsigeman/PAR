# 5. R code for plotting and statistics

```bash
#########################################################
# PAR GENE ANALYSES
#
# Analyses of ancestral PAR genes across Sylvioidea and
# outgroup species, including exon completeness, Z–W
# divergence, singleton variation and coverage analyses.
#########################################################

library(doBy)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggthemes)
library(RColorBrewer)
library(ggExtra)
library(ggpubr)
library(phytools) 
library(ggtree)
library(ape)
library(tidyr)
library(viridis)
library(extrafont)
library(fontcm)
font_import()
font_install('fontcm')
library(RColorBrewer)
library(gridExtra)
library(cowplot)
setwd("~/work/PAR/results/")

options(scipen=999)

#########################################################
# EXON COMPLETENESS ANALYSIS
#
# Calculate the proportion of each PAR exon recovered in
# species-specific genome assemblies relative to the
# zebra finch reference annotation.
#########################################################

PAR.genes.bed=read.table("datasets/allSp.PAR.genes.bed",header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
PAR.genes.bed <- plyr::rename(PAR.genes.bed, c("V1"="contig", "V2"="contig_start","V3"="contig_end", "V4"="Gene", "V5"="Trans", "V6"="exon_nr", "V7"="species"))
PAR.genes.bed$species <- sub("\\_1EV02922", "", PAR.genes.bed$species)

ZF.PAR=read.table("datasets/ZF.PAR.genes.bed",header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
ZF.PAR <- plyr::rename(ZF.PAR, c("V1"="contig", "V2"="contig_start","V3"="contig_end", "V4"="Gene", "V5"="Trans", "V6"="exon_nr", "V7"="species"))
PAR.genes.bed <- rbind(PAR.genes.bed, ZF.PAR)

ZF.PAR.genes.bed=read.table("datasets/ZF.PAR.genes.bed",header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
ZF.PAR.genes.bed <- plyr::rename(ZF.PAR.genes.bed, c("V1"="chr", "V2"="chr_start","V3"="chr_end", "V4"="Gene", "V5"="Trans", "V6"="exon_nr", "V7"="ZF"))

PAR.genes.bed <- merge(PAR.genes.bed, ZF.PAR.genes.bed, by=c("Trans", "Gene", "exon_nr"))

PAR.genes.bed$perc_found <- (PAR.genes.bed$contig_end-PAR.genes.bed$contig_start)/(PAR.genes.bed$chr_end-PAR.genes.bed$chr_start)

geneName=read.table("datasets/PAR_genes_geneID_transID_geneName.list",header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
geneName <- plyr::rename(geneName, c("V1"="name", "V2"="Gene","V3"="Trans"))
PAR.genes.bed <- merge(PAR.genes.bed, geneName, by=c("Trans", "Gene"))

pdf("~/work/PAR/results/plots/allSp_exon_completeness.pdf", height = 10, width=20)
ggplot(subset(PAR.genes.bed), aes(y=species, x=exon_nr, fill=as.factor(round(perc_found, digits = 1)))) + 
  geom_point(pch =22, size = 3) +facet_wrap(~name) + theme_tufte() + scale_fill_startrek()
dev.off()

#########################################################
# Z CHROMOSOME DATASETS
#
# Load Z-linked gene annotations, sample metadata and
# species-specific gene coordinates used in downstream
# analyses.
#########################################################

## Zebra finch gene position
ZF.Z.genes.bed=read.table("datasets/ZF.Z.genes.bed",header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
ZF.Z.genes.bed <- plyr::rename(ZF.Z.genes.bed, c("V1"="chr", "V2"="chr_start","V3"="chr_end", "V4"="Gene", "V5"="Trans", "V6"="exon_nr", "V7"="ZF"))
ZF.Z.genes.startPos <- aggregate(chr_start ~ chr + Gene + Trans, data = ZF.Z.genes.bed, max)

ZF.PAR=read.table("datasets/ZF.PAR.genes.bed",header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
ZF.PAR <- plyr::rename(ZF.PAR, c("V1"="contig", "V2"="contig_start","V3"="contig_end", "V4"="Gene", "V5"="Trans", "V6"="exon_nr", "V7"="species"))
PAR.genes <- unique(ZF.PAR$Gene)

## Sylvioidea bed file for all Z chromosome genes
allSp.Z.genes.bed=read.table("datasets/allSp.Z.genes.bed",header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
allSp.Z.genes.bed <- plyr::rename(allSp.Z.genes.bed, c("V1"="contig", "V2"="contig_start","V3"="contig_end", "V4"="Gene", "V5"="Trans", "V6"="exon_nr", "V7"="species"))

## Sample information for all Sylvioidea species
samples.sex=read.table("datasets/samples_sex.tsv",header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
samples.sex <- plyr::rename(samples.sex, c("V1"="indv", "V2"="species","V3"="sex"))

## List of Z genes with species name - used to fill in missing values for the singleton data
fillInMissing=read.table("datasets/Z.genes.allSp.list",header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
fillInMissing <- plyr::rename(fillInMissing, c("V1"="chr", "V2"="Gene","V3"="Trans", "V4"="species"))
fillInMissing <- merge(fillInMissing, ZF.Z.genes.startPos, by=c("chr", "Gene", "Trans"))

#########################################################
# PAR GENE PHYLOGENIES
#
# Extract within-species Z–W branch distances from PAR
# gene trees and summarize divergence patterns across
# genes and species.
#########################################################

setwd("~/work/PAR/mafft_ZF/")
#uncharacterized1  <- read.tree("uncharacterized1.realn.fasta.auto.gt.0.8.trim.treefile")
LMAN1  <- read.tree("LMAN1.realn.fasta.auto.gt.0.8.trim.treefile")
uncharacterized2  <- read.tree("uncharacterized2.realn.fasta.auto.gt.0.8.trim.treefile")
RAX  <- read.tree("RAX.realn.fasta.auto.gt.0.8.trim.treefile")
GRP  <- read.tree("GRP.realn.fasta.auto.gt.0.8.trim.treefile")
SEC11C  <- read.tree("SEC11C.realn.fasta.auto.gt.0.8.trim.treefile")
ZNF532 <- read.tree("ZNF532.realn.fasta.auto.gt.0.8.trim.treefile")
MALT1  <- read.tree("MALT1.realn.fasta.auto.gt.0.8.trim.treefile")
ALPK2  <- read.tree("ALPK2.realn.fasta.auto.gt.0.8.trim.treefile")
uncharacterized3 <- read.tree("uncharacterized3.realn.fasta.auto.gt.0.8.selectseqs.trim.treefile")
#MIR122 <- read.tree("MIR122.realn.fasta.auto.gt.0.8.trim.treefile")
NEDD4L <- read.tree("NEDD4L.realn.fasta.auto.gt.0.8.trim.treefile")
uncharacterized4 <- read.tree("uncharacterized4.realn.fasta.auto.gt.0.8.selectseqs.trim.treefile")
#uncharacterized5 <- read.tree("uncharacterized5.realn.fasta.auto.gt.0.8.trim.treefile")
ATP8B1  <- read.tree("ATP8B1.realn.fasta.auto.gt.0.8.trim.treefile")
NARS1 <- read.tree("NARS1.realn.fasta.auto.gt.0.8.trim.treefile")
FECH  <- read.tree("FECH.realn.fasta.auto.gt.0.8.trim.treefile")
#ONECUT2 <- read.tree("ONECUT2.realn.fasta.auto.gt.0.8.trim.treefile")
ST8SIA3  <- read.tree("ST8SIA3.realn.fasta.auto.gt.0.8.trim.treefile")
uncharacterized6 <- read.tree("uncharacterized6.realn.fasta.auto.gt.0.8.trim.treefile")
WDR7  <- read.tree("WDR7.realn.fasta.auto.gt.0.8.trim.treefile")
#TXNL1 <- read.tree("TXNL1.realn.fasta.auto.gt.0.8.trim.treefile")
#uncharacterized7 <- read.tree("uncharacterized7.realn.fasta.auto.gt.1.trim.treefile")

plot(WDR7)

dist.RAX <- reshape2::melt(cophenetic.phylo(x = RAX))
dist.RAX$Gene <- "RAX"
dist.ATP8B1 <-  reshape2::melt(cophenetic.phylo(x = ATP8B1))
dist.ATP8B1$Gene <- "ATP8B1"
dist.LMAN1 <-  reshape2::melt(cophenetic.phylo(x = LMAN1))
dist.LMAN1$Gene <- "LMAN1"
dist.uncharacterized2 <-  reshape2::melt(cophenetic.phylo(x = uncharacterized2))
dist.uncharacterized2$Gene <- "uncharacterized2"
dist.GRP <-  reshape2::melt(cophenetic.phylo(x = GRP))
dist.GRP$Gene <- "GRP"
dist.SEC11C <-  reshape2::melt(cophenetic.phylo(x = SEC11C))
dist.SEC11C$Gene <- "SEC11C"
dist.ZNF532<-  reshape2::melt(cophenetic.phylo(x = ZNF532))
dist.ZNF532$Gene <- "ZNF532"
dist.MALT1 <-  reshape2::melt(cophenetic.phylo(x = MALT1))
dist.MALT1$Gene <- "MALT1"
dist.ALPK2 <-  reshape2::melt(cophenetic.phylo(x = ALPK2))
dist.ALPK2$Gene <- "ALPK2"
dist.uncharacterized3 <-  reshape2::melt(cophenetic.phylo(x = uncharacterized3))
dist.uncharacterized3$Gene <- "uncharacterized3"
dist.NEDD4L <-  reshape2::melt(cophenetic.phylo(x = NEDD4L))
dist.NEDD4L$Gene <- "NEDD4L"
dist.uncharacterized4 <-  reshape2::melt(cophenetic.phylo(x = uncharacterized4))
dist.uncharacterized4$Gene <- "uncharacterized4"
dist.NARS1 <-  reshape2::melt(cophenetic.phylo(x = NARS1))
dist.NARS1$Gene <- "NARS1"
dist.ST8SIA3 <-  reshape2::melt(cophenetic.phylo(x = ST8SIA3))
dist.ST8SIA3$Gene <- "ST8SIA3"
dist.uncharacterized6 <-  reshape2::melt(cophenetic.phylo(x = uncharacterized6))
dist.uncharacterized6$Gene <- "uncharacterized6"
dist.WDR7 <-  reshape2::melt(cophenetic.phylo(x = WDR7))
dist.WDR7$Gene <- "WDR7"
dist.FECH <-  reshape2::melt(cophenetic.phylo(x = FECH))
dist.FECH$Gene <- "FECH"

dist.data <- rbind(dist.RAX, dist.ATP8B1, dist.LMAN1, dist.uncharacterized2, dist.GRP, dist.SEC11C, dist.ZNF532, dist.MALT1, 
                   dist.ALPK2, dist.uncharacterized3, dist.NEDD4L, dist.uncharacterized4,
                   dist.NARS1, dist.ST8SIA3, dist.uncharacterized6, dist.WDR7, dist.FECH)

# Remove self-comparisons and retain only pairwise distances.
dist.data <- subset(dist.data, dist.data$Var1!=dist.data$Var2)
dist.data <- dist.data %>%
  separate(Var1, c("Species1", "Species1_sex"), "_")
dist.data <- dist.data %>%
  separate(Var2, c("Species2", "Species2_sex"), "_")

# Keep within-species Z–W comparisons.
sp.comp.dist.data <- subset(dist.data, dist.data$Species1==dist.data$Species2)
sp.comp.dist.data <- sp.comp.dist.data %>% select(Species1, Gene, value)
head(sp.comp.dist.data)
sp.comp.dist.data.wide <- reshape(sp.comp.dist.data, idvar = "Species1", timevar = "Gene", direction = "wide")

sp.comp.dist.data$Species1 <- sub("HirDau", "CecDau", sp.comp.dist.data$Species1)
sp.comp.dist.data$Species1 <- ordered(sp.comp.dist.data$Species1,
                                      levels = c("TaeGut", "FicAlb", "CetCet", "AegCau", "PhyCol", "PycBar", "SylAtr", "TurAlt", "CecDau", "AcrSch", "LocLus", "CisJun", 
                                                 "SylBra", "PanBia", "EreAlp", "AlaArv"))

# Heatmap of Z–W branch distances across PAR genes.
pdf("~/work/PAR/results/PAR_heatmap.pdf")
ggplot(sp.comp.dist.data, aes(Gene, Species1, fill= value)) + 
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  geom_tile( colour = "black")  +  #scale_fill_gradient2(limits=c(0, 0.1),high = "firebrick4", low = "dodgerblue4") +
  # scale_fill_gradient(low="dodgerblue2", high="firebrick3", na.value = "grey") +
  theme_tufte()+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(axis.text.x=element_text(size=14, angle=90), axis.text.y=element_text(size=14)) + labs(title="", x="", y="", fill="") 
dev.off()

sp.comp.dist.data$type <- "Z + 4A"
sp.comp.dist.data$type[which(sp.comp.dist.data$Species1 =="AlaArv")] <- "Z + 4A + 3 + 5"
sp.comp.dist.data$type[which(sp.comp.dist.data$Species1 =="EreAlp")] <- "Z + 4A + 3"
sp.comp.dist.data$type[which(sp.comp.dist.data$Species1 =="PanBia")] <- "Z + 4A + 3"
sp.comp.dist.data$type[which(sp.comp.dist.data$Species1 =="CisJun")] <- "Z + 4A + 4"
sp.comp.dist.data$type[which(sp.comp.dist.data$Species1 =="SylBra")] <- "Z + 4A + 8"
sp.comp.dist.data$type[which(sp.comp.dist.data$Species1 =="TaeGut")] <- "Z"
sp.comp.dist.data$type[which(sp.comp.dist.data$Species1 =="FicAlb")] <- "Z"

myColors <- ifelse(levels(sp.comp.dist.data$type)=="Z" , rgb(0.1,0.1,0.7,0.5) , 
                   ifelse(levels(sp.comp.dist.data$type)=="Z + 4A", rgb(0.8,0.1,0.3,0.6),
                          "grey90" ) )

sp.comp.dist.data <- unique(sp.comp.dist.data)
sp.comp.dist.data
summaryBy(value ~ type + Species1, data=sp.comp.dist.data, keep.names=TRUE, FUN=c(median, sd, IQR))

# Split data into the reference species and other species
reference_species <- "TaeGut"
ref_data <- sp.comp.dist.data %>%
  filter(Species1 == reference_species)

# Compare species-specific divergence values against the
# zebra finch reference using Wilcoxon tests.

# Perform a two-sample t-test for each species
results <- sp.comp.dist.data %>%
  filter(Species1 != reference_species) %>%
  group_by(Species1) %>%
  summarise(
    median_value = median(value, na.rm = TRUE),  # Median across all genes
    sd_value = sd(value, na.rm = TRUE),                   # Standard deviation
    iqr_value = IQR(value, na.rm = TRUE), 
    p_value = if (n() > 1) {
      wilcox.test(value, ref_data$value)$p.value  # Wilcoxon
    } else {
      NA  # Not enough data for a test
    },
    higher_than_ref = median(value, na.rm = TRUE) > median(ref_data$value, na.rm = TRUE)  # Comparison of means
  ) %>%
  mutate(
    bonferroni_p = p.adjust(p_value, method = "bonferroni"),
    bonferroni_p_formatted = format(bonferroni_p, digits = 4, scientific = FALSE)
  )

# View the results
print(results)

###### 

long.pairwise_stats <- melt(pairwise_stats[[3]])
long.pairwise_stats <- subset(long.pairwise_stats, long.pairwise_stats$value!="NA")
outname <- sprintf("~/work/results/pairwise_distances_aug2020.tsv")
write.table(sp.comp.dist.data, file = outname, sep = "\t", quote = FALSE, row.names = F)

outname <- sprintf("~/work/PAR/results/pairwise_distances_wilcoxon_pairwise_aug2020.tsv")
write.table(long.pairwise_stats, file = outname, sep = "\t", quote = FALSE, row.names = F)

PAR.genes.bed
PAR.genes.bed.small <- unique(PAR.genes.bed %>% select(name, chr_start))
PAR.genes.bed.small <- aggregate(chr_start ~ name, PAR.genes.bed.small, function(x) min(x))
PAR.genes.bed.small <- PAR.genes.bed.small[order(PAR.genes.bed.small$chr_start),] 
PAR.genes.bed.small$order <- seq.int(nrow(PAR.genes.bed.small))

# Associate divergence estimates with PAR gene positions
# on the zebra finch Z chromosome.

sp.comp.dist.data.pos <- merge(sp.comp.dist.data, PAR.genes.bed.small, by.x="Gene", by.y = "name")
sp.comp.dist.data.pos <- sp.comp.dist.data.pos[order(sp.comp.dist.data.pos$order),] 
sp.comp.dist.data.pos$rownr <- seq.int(nrow(sp.comp.dist.data.pos))

sc.genes <- c("LMAN1","uncharacterized2", "NEDD4L", "ATP8B1", "NARS1", "FECH", "ST8SIA3", "uncharacterized6", "WDR7") 

sp.comp.dist.data.pos$symbol <- " "
sp.comp.dist.data.pos.sc <- sp.comp.dist.data.pos[ sp.comp.dist.data.pos$Gene %in% sc.genes, ]
sp.comp.dist.data.pos.sc.larks <- subset(sp.comp.dist.data.pos.sc, c(sp.comp.dist.data.pos.sc$Species1=="AlaArv"))

sp.comp.dist.data.pos.sc.larks$symbol <- "X"
sp.comp.dist.data.pos.sc.notlarks <- subset(sp.comp.dist.data.pos.sc, c(sp.comp.dist.data.pos.sc$Species1!="AlaArv")

sp.comp.dist.data.pos.nosc <- sp.comp.dist.data.pos[ ! sp.comp.dist.data.pos$Gene %in% sc.genes, ]
sp.comp.dist.data.pos.symbol <- rbind(sp.comp.dist.data.pos.sc.larks, sp.comp.dist.data.pos.sc.notlarks, sp.comp.dist.data.pos.nosc)
outname <- sprintf("~/work/PAR/results/pairwise_distances_aug2020.tsv")
write.table(sp.comp.dist.data.pos.symbol, file = outname, sep = "\t", quote = FALSE, row.names = F)

setwd("~/work/PAR/results")
pdf("PAR_heatmap_flip.pdf", height=10, width=18)
ggplot(sp.comp.dist.data.pos.symbol, aes(Species1, reorder(Gene, order), fill= value)) + 
  #  scale_fill_distiller(palette = "YlOrRd", direction = 1, limits=c(0, 0.03), na.value = "purple") +
  geom_tile( colour = "black", size = 0.2)  +  scale_fill_viridis(direction = -1, option = "viridis") + #geom_text(aes(label=symbol)) +
  #scale_fill_gradient2(limits=c(0, 0.1),high = "firebrick4", low = "dodgerblue4") +
  # scale_fill_gradient(low="dodgerblue2", high="firebrick3", na.value = "grey") +
  theme_tufte(base_size = 16)+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(axis.text.x=element_text(size=24, angle=90), axis.text.y=element_text(size=24)) + labs(title=" ", x="", y="", fill="") 
dev.off()
embed_fonts("PAR_heatmap_flip.pdf", outfile="PAR_heatmap_flip_font.pdf")

myColors <- ifelse(levels(sp.comp.dist.data$Species1)=="AlaArv" , rgb(163/255,213/255,91/255) , 
                          ifelse(levels(sp.comp.dist.data$Species1)=="CisJun", rgb(163/255,213/255,91/255),
                                 ifelse(levels(sp.comp.dist.data$Species1)=="SylBra", rgb(163/255,213/255,91/255),
                                        "gold" ) ))

# Build the plot
boxplot(data$value ~ data$names , 
        col=myColors , 
        ylab="disease" , xlab="- variety -")

setwd("~/work/PAR/results")

pdf("pairwise_distance_fusion_types3_new.pdf", width= 10, height = 5)
boxplot(value~Species1, data=sp.comp.dist.data, notch=FALSE, 
        col=myColors ,
        #col=(c("gold")),
        main=" ", xlab=" ")
dev.off()
embed_fonts("pairwise_distance_fusion_types3_new.pdf", outfile="pairwise_distance_fusion_types3_new_font.pdf")

setwd("~/work/PAR/results/")

extra.sc.sp.comp.dist.data.pos <- subset(sp.comp.dist.data.pos, c(sp.comp.dist.data.pos$Species1=="AlaArv" | sp.comp.dist.data.pos$Species1=="PanBia"
                                                                  | sp.comp.dist.data.pos$Species1=="EreAlp" | sp.comp.dist.data.pos$Species1=="SylBra" | sp.comp.dist.data.pos$Species1=="CisJun"))
extra.sc.sp.comp.dist.data.pos$type <- "Extra fusion"

no.extra.sc.sp.comp.dist.data.pos <- subset(sp.comp.dist.data.pos, c(sp.comp.dist.data.pos$Species1!="AlaArv" & sp.comp.dist.data.pos$Species1!="PanBia"
                                                                     & sp.comp.dist.data.pos$Species1!="EreAlp" &  sp.comp.dist.data.pos$Species1!="SylBra" & sp.comp.dist.data.pos$Species1!="CisJun"))
no.extra.sc.sp.comp.dist.data.pos$type <- "No extra fusion"

sex_linked_PAR.sp.comp.dist.data.pos <- subset(sp.comp.dist.data.pos, c(sp.comp.dist.data.pos$Species1=="AlaArv" |
                                                                          sp.comp.dist.data.pos$Species1=="SylBra" | sp.comp.dist.data.pos$Species1=="CisJun"))
sex_linked_PAR.sp.comp.dist.data.pos$PAR_type <- "Extra fusion"

autosomal_PAR.sp.comp.dist.data.pos <- subset(sp.comp.dist.data.pos, c(sp.comp.dist.data.pos$Species1!="AlaArv" & 
                                                                         sp.comp.dist.data.pos$Species1!="SylBra" & sp.comp.dist.data.pos$Species1!="CisJun"))
autosomal_PAR.sp.comp.dist.data.pos$PAR_type <- "No extra fusion"

pearson.a <- subset(sp.comp.dist.data.pos, c(sp.comp.dist.data.pos$Species1=="TaeGut" | sp.comp.dist.data.pos$Species1=="PhyCol" | sp.comp.dist.data.pos$Species1=="AegCau"
                                             | sp.comp.dist.data.pos$Species1=="AcrSch" | sp.comp.dist.data.pos$Species1=="PanBia" |  sp.comp.dist.data.pos$Species1=="CecDau" |
                                               sp.comp.dist.data.pos$Species1=="LocLus"))

pearson.a <-unique(pearson.a[ , 1:5 ] )

pearson.b <- subset(sp.comp.dist.data.pos, c(sp.comp.dist.data.pos$Species1=="PycBar" | sp.comp.dist.data.pos$Species1=="CetCet" | sp.comp.dist.data.pos$Species1=="SylAtr" 
                                             | sp.comp.dist.data.pos$Species1=="TurAlt" | sp.comp.dist.data.pos$Species1=="EreAlp" | sp.comp.dist.data.pos$Species1=="FicAlb"))
pearson.b <-unique(pearson.b[ , 1:5 ] )

pearson.c <- subset(sp.comp.dist.data.pos, c(sp.comp.dist.data.pos$Species1=="AlaArv" | sp.comp.dist.data.pos$Species1=="CisJun" | sp.comp.dist.data.pos$Species1=="SylBra"))
pearson.c <-unique(pearson.c[ , 1:5 ] )

sp.comp.dist.data.pos$species_name <- "null"
sp.comp.dist.data.pos$species_name[sp.comp.dist.data.pos$Species1=="TaeGut"]="Taeniopygia guttata (Group I: Z)"
sp.comp.dist.data.pos$species_name[sp.comp.dist.data.pos$Species1=="AlaArv"]="Alauda arvensis (Group VII: Z;4A,3,5)*"
sp.comp.dist.data.pos$species_name[sp.comp.dist.data.pos$Species1=="EreAlp"]="Eremophila alpestris (Group VI: Z;4A,3)"
sp.comp.dist.data.pos$species_name[sp.comp.dist.data.pos$Species1=="PanBia"]="Panurus biarmicus  (Group V: Z;4A;3;5)"
sp.comp.dist.data.pos$species_name[sp.comp.dist.data.pos$Species1=="SylBra"]="Sylvietta brachyura (Group IV: Z;4A;8)*"
sp.comp.dist.data.pos$species_name[sp.comp.dist.data.pos$Species1=="AegCau"]="Aegithalos caudatus (Group II: Z;4A)"
sp.comp.dist.data.pos$species_name[sp.comp.dist.data.pos$Species1=="AcrSch"]="Acrocephalus schoenobaenus (Group II: Z;4A)"
sp.comp.dist.data.pos$species_name[sp.comp.dist.data.pos$Species1=="CetCet"]="Cettia cetti (Group II: Z;4A)"
sp.comp.dist.data.pos$species_name[sp.comp.dist.data.pos$Species1=="CisJun"]="Cisticola juncidis (Group III: Z;4A;4)*"
sp.comp.dist.data.pos$species_name[sp.comp.dist.data.pos$Species1=="PhyCol"]="Phylloscopus collybita (Group II: Z;4A)"
sp.comp.dist.data.pos$species_name[sp.comp.dist.data.pos$Species1=="CecDau"]="Cecropis daurica (Group II: Z;4A)"
sp.comp.dist.data.pos$species_name[sp.comp.dist.data.pos$Species1=="FicAlb"]="Ficedula albicollis (Group I: Z)"
sp.comp.dist.data.pos$species_name[sp.comp.dist.data.pos$Species1=="PycBar"]="Pycnonotus barbatus (Group II: Z;4A)"
sp.comp.dist.data.pos$species_name[sp.comp.dist.data.pos$Species1=="SylAtr"]="Sylvia atricapilla (Group II: Z;4A)"
sp.comp.dist.data.pos$species_name[sp.comp.dist.data.pos$Species1=="LocLus"]="Locustella luscinioides (Group II: Z;4A)"
sp.comp.dist.data.pos$species_name[sp.comp.dist.data.pos$Species1=="TurAlt"]="Argya altirostris (Group II: Z;4A)"

# Test for correlations between genomic position and
# Z–W divergence within species.

setwd("~/work/PAR/results")
pdf("pairwise_distance_spearman.allSp.pdf", width= 14, height = 12)
ggscatter(sp.comp.dist.data.pos, x = "chr_start", y = "value", fill = "Species1", 
          size = 3, shape = 21, facet.by = "species_name", combine = T,xlab = "Z chromosome position", 
          add = "reg.line", conf.int = F,ylab = "pairwise distance", repel = T,
          add.params = list(color = "Species1", fill = "grey", size = 0.75),
          title = NULL, show.legend.text = FALSE, 
          cor.coef = TRUE,
          scales = "free_y",
          cor.coeff.args = list(method = "spearman", cor.coef.name = "rho"))
dev.off()
embed_fonts("pairwise_distance_spearman.allSp.pdf", outfile="pairwise_distance_spearman.allSp_font.pdf")

setwd("~/work/PAR/results/")

sp_name <- unique(sp.comp.dist.data.pos$species_name)
sp.comp.dist.data.pos$chr_start_mb <- sp.comp.dist.data.pos$chr_start/1000000
cols <- viridis(17, alpha = 1, begin = 0, end = 1, option = "D")

scaleFUN <- function(x) sprintf("%.3f", x)
plot_list = list()
for (i in 1:16) {
  data_i <- subset(sp.comp.dist.data.pos, sp.comp.dist.data.pos$species_name==sp_name[i])
  y_limit <- max(data_i$value)+0.001
  p = ggscatter(data_i, x = "chr_start_mb", y = "value", 
                size = 9, shape = 21,xlab = " ", fill = "#A3D55B", 
                add = "reg.line", conf.int = F,ylab = " ",
                title = sp_name[i], show.legend.text = FALSE, show.legend = FALSE, 
                cor.coef = TRUE,
                cor.coeff.args = list(method = "spearman", cor.coef.name = "rho", size=6, family="Helvetica"), 
                cor.coef.coord = c(0.05, y_limit),
                repel = T) + scale_y_continuous(labels=scaleFUN)  
  p = p +  theme(plot.title=element_text(size=20,face="bold", hjust = 0.5),
                 axis.text=element_text(size=20),
                 axis.title=element_text(size=8), text = element_text(size=22)) + border() 
  plot_list[[i]] = p
}

# Save plots to tiff. Makes a separate file for each plot.
for (i in 1:16) {
  file_name = paste("spearman_position_", sp_name[i], ".pdf", sep="")
  file_name_font = paste("font_iris_plot_", i, ".pdf", sep="")
  pdf(file_name, width= 7, height = 5, family="Times")
  print(plot_list[[i]])
  dev.off()
  embed_fonts(file_name, outfile=file_name_font)
}

# Another option: create pdf where each page is a separate plot.
pdf("plots.pdf")
for (i in 1:3) {
  print(plot_list[[i]])
}
dev.off()

plot_list[1] # Acr
plot_list[2] # Syl
plot_list[3] # AlaArv
plot_list[4] # Ere
plot_list[5] # Pan
plot_list[6] # ZF
plot_list[7] # Fic
plot_list[8] # Aeg
plot_list[9] #Cet
plot_list[10] #Cis
plot_list[11] #Phy
plot_list[12] # Loc
plot_list[13] #Cec
plot_list[14] # SylAtr
plot_list[15] # Argya
plot_list[16] # Pyc

#ZF, FC + 1
pdf("~/work/PAR/results/plots/Fig5_row1.pdf", height=4, width = 28)
grid.arrange(grobs=c(plot_list[6], plot_list[7], plot_list[1]), 
             nrow = 1, 
             labels = c("A","B", "C"))
dev.off()

# CetCet, Aeg, Phy
pdf("~/work/PAR/results/plots/Fig5_row2.pdf", height=4, width = 28)
grid.arrange(grobs=c(plot_list[9], plot_list[8], plot_list[11]), 
             nrow = 1, 
             labels = c("A","B", "C"))
dev.off()

#Pyc SylAtr, Argya
pdf("~/work/PAR/results/plots/Fig5_row3.pdf", height=4, width = 28)
grid.arrange(grobs=c(plot_list[16], plot_list[14], plot_list[15]), 
             nrow = 1, 
             labels = c("A","B", "C"))
dev.off()

#Pyc SylAtr, Argya
pdf("~/work/PAR/results/plots/Fig5_row4.pdf", height=4, width = 28)
grid.arrange(grobs=c(plot_list[13], plot_list[1], plot_list[12]), 
             nrow = 1, 
             labels = c("A","B", "C"))
dev.off()

# Cis SylBra Pan
pdf("~/work/PAR/results/plots/Fig5_row5.pdf", height=4, width = 28)
grid.arrange(grobs=c(plot_list[10], plot_list[2], plot_list[5]), 
             nrow = 1, 
             labels = c("A","B", "C"))
dev.off()

# Ere Alauda
pdf("~/work/PAR/results/plots/Fig5_row6.pdf", height=4, width = 28)
grid.arrange(grobs=c(plot_list[4], plot_list[3], plot_list[4]), 
             nrow = 1, 
             labels = c("A","B", "C"))
dev.off()

#########################################################
# SINGLETON ANALYSES
#
# Compare female and male singleton counts among PAR
# genes within species.
#########################################################

Z.singletons=read.table("~/work/PAR/results/PAR_FicAlb/allSp.Z.genes.singletons.S.sum.krakenSNPs.out",header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
Z.singletons <- plyr::rename(Z.singletons, c("V1"="nr.singletons", "V2"="Gene","V3"="indv", "V4"="species", "V5"="sex"))
specieslist <- unique(Z.singletons$species)
samples.sex.15 <- samples.sex[ samples.sex$species %in% specieslist, ]
samples.Z.gene.combination <- merge(samples.sex.15, ZF.Z.genes.startPos)

Z.singletons <- merge(Z.singletons, samples.Z.gene.combination, all=TRUE, by=c("Gene", "sex", "indv", "species"))
Z.singletons$species <- sub("\\_1EV02922", "", Z.singletons$species)
Z.singletons$species <- sub("HirDau", "CecDau", Z.singletons$species)
Z.singletons[is.na(Z.singletons)] <- 0
PAR.genes.bed.small2 <- unique(PAR.genes.bed %>% select(Gene, name, chr_start))
PAR.genes.bed.small2 <- aggregate(chr_start ~ name + Gene, PAR.genes.bed.small2, function(x) min(x))
PAR.genes.bed.small2 <- PAR.genes.bed.small2[order(PAR.genes.bed.small2$chr_start),] 
PAR.genes.bed.small2$order <- seq.int(nrow(PAR.genes.bed.small2))
Z.singletons.PAR <- merge(Z.singletons, PAR.genes.bed.small2, by=c("Gene"))
Z.singletons.PAR$species <- ordered(Z.singletons.PAR$species,
                                    levels = c("TaeGut", "FicAlb", "CetCet", "AegCau", "PhyCol", "PycBar", "SylAtr", "TurAlt", "CecDau", "LocLus", "AcrSch", "CisJun", 
                                               "SylBra", "PanBia", "EreAlp", "AlaArv"))

Z.singletons.PAR$species_name <- "null"
Z.singletons.PAR$species_name[Z.singletons.PAR$species=="TaeGut"]="Taeniopygia guttata"
Z.singletons.PAR$species_name[Z.singletons.PAR$species=="AlaArv"]="Alauda arvensis"
Z.singletons.PAR$species_name[Z.singletons.PAR$species=="EreAlp"]="Eremophila alpestris"
Z.singletons.PAR$species_name[Z.singletons.PAR$species=="PanBia"]="Panurus biarmicus"
Z.singletons.PAR$species_name[Z.singletons.PAR$species=="SylBra"]="Sylvietta brachyura"
Z.singletons.PAR$species_name[Z.singletons.PAR$species=="AegCau"]="Aegithalos caudatus"
Z.singletons.PAR$species_name[Z.singletons.PAR$species=="AcrSch"]="Acrocephalus schoenobaenus"
Z.singletons.PAR$species_name[Z.singletons.PAR$species=="CetCet"]="Cettia cetti"
Z.singletons.PAR$species_name[Z.singletons.PAR$species=="CisJun"]="Cisticola juncidis"
Z.singletons.PAR$species_name[Z.singletons.PAR$species=="PhyCol"]="Phylloscopus collybita"
Z.singletons.PAR$species_name[Z.singletons.PAR$species=="CecDau"]="Cecropis daurica"
Z.singletons.PAR$species_name[Z.singletons.PAR$species=="FicAlb"]="Ficedula albicollis"
Z.singletons.PAR$species_name[Z.singletons.PAR$species=="PycBar"]="Pycnonotus barbatus"
Z.singletons.PAR$species_name[Z.singletons.PAR$species=="SylAtr"]="Sylvia atricapilla"
Z.singletons.PAR$species_name[Z.singletons.PAR$species=="LocLus"]="Locustella luscinioides"
Z.singletons.PAR$species_name[Z.singletons.PAR$species=="TurAlt"]="Argya altirostris"

sc.genes2 <- unique(sp.comp.dist.data.pos.symbol$Gene)

Z.singletons.PAR <- Z.singletons.PAR[ Z.singletons.PAR$name %in% sc.genes2, ]

pdf("~/work/PAR/results/PAR_heatmap_singletons_per_sex.pdf", height = 10, width =14)
ggplot(Z.singletons.PAR, aes(reorder(name, order),species_name, fill= nr.singletons)) + 
  geom_tile( colour = "black", size = 0.5)  +  scale_fill_viridis(direction = -1, option = "viridis") +
  theme_tufte(base_family="Helvetica", base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(axis.text.x=element_text(size=18, angle=90), axis.text.y=element_text(size=18)) + labs(title="Per sex number of singletons", x="", y="", fill="") + facet_wrap(~sex, ncol = 1)
dev.off()

# Convert singleton counts to paired female/male format.
Z.singletons.PAR_wide <- reshape2::dcast(Z.singletons.PAR, species + Gene + Trans + order + name ~ sex, value.var="nr.singletons")

pdf("~/work/PAR/results/PAR_heatmap_singletons_sex_diff.pdf", height = 10, width =10)
ggplot(Z.singletons.PAR_wid, aes(reorder(name, order),species, fill= female-male)) + 
  geom_tile( colour = "black", size = 0.5)  +  scale_fill_viridis(direction = -1, option = "viridis") +
  theme_tufte(base_size = 18, base_family = "serif")+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(axis.text.x=element_text(size=18, angle=90), axis.text.y=element_text(size=18)) + labs(title="Female-male difference in singletons", x="", y="", fill="") 
dev.off()

# Paired Wilcoxon tests comparing female and male
# singleton counts within species.

results <- Z.singletons.PAR_wide %>%
  group_by(species) %>%
  summarise(
    p_value = if (n() > 1) {
      wilcox.test(female, male, paired = TRUE)$p.value
    } else {
      NA  # Not enough data for the test
    },
    female_higher = mean(female > male)  # Proportion of cases where female > male
  )

results <- Z.singletons.PAR_wide %>%
  group_by(species) %>%
  summarise(
    median_female = median(female, na.rm = TRUE),  # Mean value for females
    median_male = median(male, na.rm = TRUE),      # Mean value for males
    p_value = if (n() > 1) {
      wilcox.test(female, male, paired = TRUE)$p.value
    } else {
      NA  # Not enough data for the test
    },
    female_higher = mean(female > male, na.rm = TRUE)  # Proportion of cases where female > male
  ) %>%
  mutate(
    bonferroni_p = p.adjust(p_value, method = "bonferroni")
  )
summaryBy(value ~ type + Species1, data=sp.comp.dist.data, keep.names=TRUE, FUN=c(median, sd, IQR))

new_table <- summaryBy(nr.singletons ~ sex + species_name, data=Z.singletons.PAR, keep.names=TRUE, FUN=c(median, sd, IQR, min, max))
outname <- sprintf("~/work/PAR/results/singleton_stats.tsv")
write.table(new_table, file = outname, sep = "\t", quote = FALSE, row.names = F)

Z.singletons.PAR_wide
kruskal.test(value~Species1, data = sp.comp.dist.data) # where y1 is numeric and A is a factor
pairwise.wilcox.test(log2(sp.comp.dist.data$value), sp.comp.dist.data$Species1,
                     p.adjust.method = "BH")

pairwise_stats <- pairwise.wilcox.test(sp.comp.dist.data$value, sp.comp.dist.data$Species1,
                                       p.adjust.method = "bonf")

long.pairwise_stats <- melt(pairwise_stats[[3]])
long.pairwise_stats <- subset(long.pairwise_stats, long.pairwise_stats$value!="NA")
outname <- sprintf("~/work/results/pairwise_distances_aug2020.tsv")
write.table(sp.comp.dist.data, file = outname, sep = "\t", quote = FALSE, row.names = F)

#########################################################
# COVERAGE ANALYSES
#
# Analyse female-to-male coverage ratios for Z-linked and
# PAR genes across species.
#########################################################

# Depth stats
Z.depth=read.table("~/work/PAR/results/datasets/allSp.Z.genes.bamstat04.geneInfo.out",header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
Z.depth <- plyr::rename(Z.depth, c("V1"="indv", "V2"="min.depth","V3"="max.depth","V4"="avg.depth", "V5"="median.depth","V6"="species","V7"="Gene","V8"="Trans"))
Z.depth.mean <- aggregate(avg.depth ~ species + Gene + indv + Trans + indv, Z.depth, mean)
Z.depth.mean.combination <- merge(Z.depth.mean, samples.Z.gene.combination)
Z.depth.mean.combination_wide <- reshape2::dcast(Z.depth.mean.combination, species + Gene + Trans + chr_start ~ sex, value.var="avg.depth")

Z.depth.mean.combination_wide <- subset(Z.depth.mean.combination_wide, c(Z.depth.mean.combination_wide$female>5 & Z.depth.mean.combination_wide$female<80))
Z.depth.mean.combination_wide <- subset(Z.depth.mean.combination_wide, c(Z.depth.mean.combination_wide$male>5 & Z.depth.mean.combination_wide$male<80))
aggregate(male ~ species, data = Z.depth.mean.combination_wide, mean)

# Female-to-male coverage ratio per gene.

### Calculate coverage ratio and remove outliers
Z.depth.mean.combination_wide$cov_ratio <- Z.depth.mean.combination_wide$female/Z.depth.mean.combination_wide$male
Z.depth.mean.combination_wide <- subset(Z.depth.mean.combination_wide, c(Z.depth.mean.combination_wide$cov_ratio < 2.5))
Z.depth.mean.combination_wide.median <- aggregate(cov_ratio ~ species, data = Z.depth.mean.combination_wide, median)
Z.depth.mean.combination_wide.median <- plyr::rename(Z.depth.mean.combination_wide.median, c("cov_ratio"="species_median"))
Z.depth.mean.combination_wide <- merge(Z.depth.mean.combination_wide, Z.depth.mean.combination_wide.median, by=c("species"))

# Normalize coverage ratios within species and centre
# values around zero.
Z.depth.mean.combination_wide$cov_ratio_scaled <- Z.depth.mean.combination_wide$cov_ratio / Z.depth.mean.combination_wide$species_median
Z.depth.mean.combination_wide <- na.omit(Z.depth.mean.combination_wide) 
Z.depth.mean.combination_wide$cov_ratio_scaled <- scale(Z.depth.mean.combination_wide$cov_ratio_scaled, center = TRUE, scale = FALSE)

Z.depth.mean.combination_wide.PAR <- merge(Z.depth.mean.combination_wide, PAR.genes.bed.small2, by=c("Gene"))

Z.depth.mean.combination_wide.PAR <- Z.depth.mean.combination_wide.PAR[ Z.depth.mean.combination_wide.PAR$name %in% sc.genes2, ]
Z.depth.mean.combination_wide.PAR$species <- sub("HirDau", "CecDau", Z.depth.mean.combination_wide.PAR$species)

Z.depth.mean.combination_wide$species_name <- "null"
Z.depth.mean.combination_wide$species_name[sp.comp.dist.data.pos$species=="TaeGut"]="Taeniopygia guttata (Group I: Z)"
Z.depth.mean.combination_wide$species_name[Z.depth.mean.combination_wide$species=="AlaArv"]="Alauda arvensis (Group VII: Z;4A,3,5)"
Z.depth.mean.combination_wide$species_name[Z.depth.mean.combination_wide$species=="EreAlp"]="Eremophila alpestris (Group VI: Z;4A,3)"
Z.depth.mean.combination_wide$species_name[Z.depth.mean.combination_wide$species=="PanBia"]="Panurus biarmicus  (Group V: Z;4A;3;5)"
Z.depth.mean.combination_wide$species_name[Z.depth.mean.combination_wide$species=="SylBra"]="Sylvietta brachyura (Group IV: Z;4A;8)"
Z.depth.mean.combination_wide$species_name[Z.depth.mean.combination_wide$species=="AegCau"]="Aegithalos caudatus (Group II: Z;4A)"
Z.depth.mean.combination_wide$species_name[Z.depth.mean.combination_wide$species=="AcrSch"]="Acrocephalus schoenobaenus (Group II: Z;4A)"
Z.depth.mean.combination_wide$species_name[Z.depth.mean.combination_wide$species=="CetCet"]="Cettia cetti (Group II: Z;4A)"
Z.depth.mean.combination_wide$species_name[Z.depth.mean.combination_wide$species=="CisJun"]="Cisticola juncidis (Group III: Z;4A;4)"
Z.depth.mean.combination_wide$species_name[Z.depth.mean.combination_wide$species=="PhyCol"]="Phylloscopus collybita (Group II: Z;4A)"
Z.depth.mean.combination_wide$species_name[Z.depth.mean.combination_wide$species=="HirDau"]="Cecropis daurica (Group II: Z;4A)"
Z.depth.mean.combination_wide$species_name[Z.depth.mean.combination_wide$species=="FicAlb"]="Ficedula albicollis (Group I: Z)"
Z.depth.mean.combination_wide$species_name[Z.depth.mean.combination_wide$species=="PycBar"]="Pycnonotus barbatus (Group II: Z;4A)"
Z.depth.mean.combination_wide$species_name[Z.depth.mean.combination_wide$species=="SylAtr_1EV02922"]="Sylvia atricapilla (Group II: Z;4A)"
Z.depth.mean.combination_wide$species_name[Z.depth.mean.combination_wide$species=="LocLus"]="Locustella luscinioides (Group II: Z;4A)"
Z.depth.mean.combination_wide$species_name[Z.depth.mean.combination_wide$species=="TurAlt"]="Argya altirostris (Group II: Z;4A)"

# Coverage ratios along the Z chromosome.
# Grey shading indicates the ancestral PAR.

options(scipen=999)
pdf("~/work/PAR/results/plots/FigureS1_PAR_genes_depth.pdf", height = 10, width =15)
ggplot(subset(Z.depth.mean.combination_wide, chr_start<10000000), aes(x = chr_start/1000000, y = cov_ratio_scaled, fill = species_name)) +
  geom_rect(mapping=aes(xmin=0, xmax=0.450000, ymin=-Inf, ymax=Inf), fill="grey" , alpha=0.5) +
  geom_point(pch = 21) + facet_wrap(~species_name) + theme_minimal() +geom_hline(yintercept = 0.0, linetype = "dashed") + 
  labs(x = "Z chromosome position (Mb)", y = "Female-to-male genome coverage ratio") + 
  theme_bw(base_family="Helvetica", base_size = 12) + theme(legend.position="none") 
dev.off()

pdf("~/work/PAR/results/plots/PAR_genes_depth.female.pdf", height = 10, width =15)
ggplot(subset(Z.depth.mean.combination_wide, chr_start<10000000), aes(x = chr_start, y = female, fill = species)) +
  geom_rect(mapping=aes(xmin=0, xmax=450000, ymin=-Inf, ymax=Inf), fill="grey" , alpha=0.5) +
  geom_point(pch = 21) + facet_wrap(~species) + theme_minimal() +geom_hline(yintercept = 0.5)
dev.off()

pdf("~/work/PAR/results/plots/PAR_genes_depth.male.pdf", height = 10, width =15)
ggplot(subset(Z.depth.mean.combination_wide, chr_start<10000000), aes(x = chr_start, y = female, fill = species)) +
  geom_rect(mapping=aes(xmin=0, xmax=450000, ymin=-Inf, ymax=Inf), fill="grey" , alpha=0.5) +
  geom_point(pch = 21) + facet_wrap(~species) + theme_minimal() +geom_hline(yintercept = 0.5)
dev.off()

Z.depth.mean.combination_wide.PAR.small <- Z.depth.mean.combination_wide.PAR %>% select(species, name, female, male)
Z.depth.mean.combination_wide.PAR.small_wide <- Z.depth.mean.combination_wide.PAR.small %>%
  pivot_wider(
    id_cols = name,
    names_from = species,
    values_from = c(female, male),
    names_glue = "{species}_{.value}"
  )

gene_counts <- Z.depth.mean.combination_wide %>%
  group_by(species) %>%
  summarise(
    total_genes = n(),
    genes_first_10Mb = sum(chr_start <= 10000000, na.rm = TRUE)
  )

gene_counts_PAR <- Z.depth.mean.combination_wide.PAR %>%
  group_by(species) %>%
  summarise(
    PAR_genes = n()
  )

merge(gene_counts, gene_counts_PAR, by = "species")

outname <- sprintf("~/work/PAR/results/cov_values_PAR_genes.tsv")
write.table(Z.depth.mean.combination_wide.PAR.small_wide, file = outname, sep = "\t", quote = FALSE, row.names = F)

```