
############## Bile acid analysis ##############
library(ggplot2)
library(ggpubr)
library(tidyr)

bile <- read.table("bile_acid_profiles.txt")

#subset samples to include only those with bile acid profiles
ps2.bile.ibd <- subset_samples(ps2.ibd, rownames(sample_data(ps2.ibd)) %in% bile$Sample)
bile <- bile[bile$Sample %in% rownames(sample_data(ps2.bile.ibd)),]

all(rownames(sample_data(ps2.bile.ibd)) %in% bile$Sample)
all(bile$Sample %in% rownames(sample_data(ps2.bile.ibd)))

#match samples and ensure samples matched
bile <- bile[match(rownames(sample_data(ps2.bile.ibd)), bile$Sample),]
all(bile$Sample == rownames(sample_data(ps2.bile.ibd)))

#rename rows and remove categorical column
rownames(bile) <- bile$Sample
bile <- bile[,-c(1:2)]
#remove any columns with NA
bile.ibd <- bile.ibd[,!(is.na(colSums(bile.ibd)))]

#create data frame
bile.ibd.df <- cbind(sample_data(ps2.bile.ibd), bile.ibd)



##### 16S PCoA with bile acid gradients #####


if(!(is.null(rownames(bile.ibd))) & !(is.null(rownames(sample_data(ps2.bile.ibd))))){
  if(all(rownames(bile.ibd) == rownames(sample_data(ps2.bile.ibd)))){
    ps.transf <- ps2.bile.ibd
    sample_data(ps.transf) <- cbind(sample_data(ps.transf), bile.ibd)
    ps.transf <- transform_sample_counts(ps.transf, function(x) {x/sum(x)})
  }
}



####### Figure 5A #######
#Clusters PCoA

pal <- c("red3", "grey25", "chartreuse2")

out.donor.bray <- ordinate(ps.transf, method = "MDS", distance = "bray")

bray.df <- cbind(out.donor.bray$vectors[,1:2], sample_data(ps.transf))
bray.df$Axis.1 <- bray.df$Axis.1*-1 #to keep orientated with prior PCoA for IBD
eval.1 <- (out.donor.bray$values$Eigenvalues[1]/sum(out.donor.bray$values$Eigenvalues))*100
eval.2 <- (out.donor.bray$values$Eigenvalues[2]/sum(out.donor.bray$values$Eigenvalues))*100
p2 <- ggplot(bray.df, aes(x=Axis.1, y=Axis.2, colour=Cluster, fill=Cluster)) + geom_point(size=2) + scale_fill_manual("Cluster", values = pal) + 
  scale_colour_manual("Cluster", values = pal) + stat_ellipse(level = 0.8) + theme_classic() + 
  xlab(paste("PCoA 1 [", round(eval.1, 2), "%]", sep = "")) + ylab(paste("PCoA 2 [", round(eval.2, 2), "%]", sep = "")) + 
  theme(legend.position = "right", text = element_text(size=16), axis.ticks.length = unit(0.2, "cm"), 
        axis.text = element_text(size=14), plot.title = element_text(hjust = 0.5), 
        axis.line = element_line(colour = 'black', size = 1.05), axis.ticks = element_line(colour = 'black', size = 1.05)) +
  ggtitle("Cluster")

#plot
p2


####### Figure 5B #######
#Primary:Secondary BAs PCoA

p.bray <- ggplot(bray.df, aes(x=Axis.1, y=Axis.2, colour=log10(Primary_BA.Secondary_BA), group=Cluster)) + 
  geom_point(size=2, alpha=1) + stat_ellipse(level = 0.8, color="grey25", lty=2) + theme_classic() + 
  xlab(paste("PCoA 1 [", round(eval.1, 2), "%]", sep = "")) + 
  ylab(paste("PCoA 2 [", round(eval.2, 2), "%]", sep = "")) + 
  theme(legend.position = "right", text = element_text(size=16), axis.ticks.length = unit(0.2, "cm"), 
        axis.text = element_text(size=14), plot.title = element_text(hjust = 0.5), 
        axis.line = element_line(colour = 'black', size = 1.05), 
        axis.ticks = element_line(colour = 'black', size = 1.05)) +
  ggtitle("Primary:Secondary Bile Acid Ratio (log)")
p3 <- p.bray + scale_fill_gradient2(parse(text=paste(paste("1^o", "2^o~", sep = ":"), "BAs", sep = "")),
                                    low = "blue",
                                    mid = "white",
                                    high = "red",
                                    midpoint = 0,
                                    space = "Lab",
                                    na.value = "grey50",
                                    guide = "colourbar",
                                    aesthetics = "colour"
)

#plot
p3


####### Figure 5C #######
#Unconjugated BAs PCoA

p.bray <- ggplot(bray.df, aes(x=Axis.1, y=Axis.2, colour=BA_nonconj.AB_Tot, group=Cluster)) + 
  geom_point(size=2, alpha=1) + stat_ellipse(level = 0.8, color="grey25", lty=2) + theme_classic() + 
  xlab(paste("PCoA 1 [", round(eval.1, 2), "%]", sep = "")) + 
  ylab(paste("PCoA 2 [", round(eval.2, 2), "%]", sep = "")) + 
  theme(legend.position = "right", text = element_text(size=16), axis.ticks.length = unit(0.2, "cm"), 
        axis.text = element_text(size=14), plot.title = element_text(hjust = 0.5), 
        axis.line = element_line(colour = 'black', size = 1.05), 
        axis.ticks = element_line(colour = 'black', size = 1.05)) +
  ggtitle("Proportion of Unconjugated Bile Acids")
p4 <- p.bray + scale_fill_gradient("Proportion\nunconjugated",
                                   low = "darkorange",
                                   
                                   high = "grey89",
                                   
                                   space = "Lab",
                                   na.value = "grey50",
                                   guide = "colourbar",
                                   aesthetics = "colour"
)
#plot
p4

####### Figure 5D #######
#Sulphated BAs PCoA

p.bray <- ggplot(bray.df, aes(x=Axis.1, y=Axis.2, colour=Sulfo.BA_Total, group=Cluster)) + 
  geom_point(size=2, alpha=1) + stat_ellipse(level = 0.8, color="grey25", lty=2) + theme_classic() + 
  xlab(paste("PCoA 1 [", round(eval.1, 2), "%]", sep = "")) + 
  ylab(paste("PCoA 2 [", round(eval.2, 2), "%]", sep = "")) + 
  theme(legend.position = "right", text = element_text(size=16), axis.ticks.length = unit(0.2, "cm"), 
        axis.text = element_text(size=14), plot.title = element_text(hjust = 0.5), 
        axis.line = element_line(colour = 'black', size = 1.05), 
        axis.ticks = element_line(colour = 'black', size = 1.05)) +
  ggtitle("Proportion of Sulphated Bile Acids")
p5 <- p.bray + scale_fill_gradient("Proportion\nsulphated",
                                   low = "lightblue",
                                   
                                   high = "maroon",
                                   
                                   space = "Lab",
                                   na.value = "grey50",
                                   guide = "colourbar",
                                   aesthetics = "colour"
)

#plot
p5







####### Figure 5E #######

# Bile acid PCA #


#replace min
addMin <- function(df){
  cm <- apply(df, 2, function (x) min(x[x>0])/5)
  df2 <- df
  for(i in 1:length(cm)){
    df2[,i][df2[,i] == 0] <- cm[i]
  }
  return(df2)
}
bile.ibd.comp <- addMin(bile.ibd)

#log
bile.ibd.comp <- log2(bile.ibd.comp)

#center
bile.ibd.comp <- sweep(bile.ibd.comp,2,colMeans(bile.ibd.comp),"-")



library(ade4)
library(factoextra)
#all
sampleTable <- sample_data(ps2.bile.ibd)

gene.pac <- dudi.pca(bile.ibd.comp, nf=3, scale = FALSE, center = FALSE, scannf = FALSE)
pca.scrs <- gene.pac$li

if(!(is.null(rownames(bile.ibd.df))) & !(is.null(rownames(pca.scrs)))){
  if(all(rownames(bile.ibd.df) == rownames(pca.scrs))){
    scrs <- cbind(pca.scrs, bile.ibd.df)
  }
}

#Clusters PCA

pal <- c("red3", "grey25", "chartreuse2")



eval.1 <- round(gene.pac$eig[1]/sum(gene.pac$eig),3)*100
eval.2 <- round(gene.pac$eig[2]/sum(gene.pac$eig),3)*100

p2 <- ggplot(pca.scrs, aes(x=Axis1, y=Axis2, colour=Cluster, fill=Cluster)) + geom_point(size=2) + scale_fill_manual("Cluster", values = pal) + 
  scale_colour_manual("Cluster", values = pal) + stat_ellipse(level = 0.8) + theme_classic() + 
  xlab(paste("PCA 1 [", round(eval.1, 2), "%]", sep = "")) + ylab(paste("PCA 2 [", round(eval.2, 2), "%]", sep = "")) + 
  theme(legend.position = "right", text = element_text(size=16), axis.ticks.length = unit(0.2, "cm"), 
        axis.text = element_text(size=14), plot.title = element_text(hjust = 0.5), 
        axis.line = element_line(colour = 'black', size = 1.05), axis.ticks = element_line(colour = 'black', size = 1.05)) +
  ggtitle("Cluster")

p2


####### Figure 5F #######
#Primary:Secondary BAs PCA

p.bray <- ggplot(pca.scrs, aes(x=Axis1, y=Axis2, colour=log10(Primary_BA.Secondary_BA), group=Cluster)) + 
  geom_point(size=2, alpha=1) + stat_ellipse(level = 0.8, color="grey25", lty=2) + theme_classic() + 
  xlab(paste("PCA 1 [", round(eval.1, 2), "%]", sep = "")) + 
  ylab(paste("PCA 2 [", round(eval.2, 2), "%]", sep = "")) + 
  theme(legend.position = "right", text = element_text(size=16), axis.ticks.length = unit(0.2, "cm"), 
        axis.text = element_text(size=14), plot.title = element_text(hjust = 0.5), 
        axis.line = element_line(colour = 'black', size = 1.05), 
        axis.ticks = element_line(colour = 'black', size = 1.05)) +
  ggtitle("Primary:Secondary Bile Acid Ratio (log)")
p3 <- p.bray + scale_fill_gradient2(parse(text=paste(paste("1^o", "2^o~", sep = ":"), "BAs", sep = "")),
                                    low = "blue",
                                    mid = "white",
                                    high = "red",
                                    midpoint = 0,
                                    space = "Lab",
                                    na.value = "grey50",
                                    guide = "colourbar",
                                    aesthetics = "colour"
)

p3




####### Figure 5G #######
#Barplots

#Bile acid composition
bile.ibd.prop <- bile.ibd[,1:27]
bile.ibd.prop <- sweep(bile.ibd.prop, 1, rowSums(bile.ibd.prop), "/")

if(!(is.null(rownames(bile.ibd))) & !(is.null(rownames(sample_data(ps2.bile.ibd))))){
  if(all(rownames(bile.ibd) == rownames(sample_data(ps2.bile.ibd)))){
    bile.ibd.df <- cbind(SampleID = rownames(sample_data(ps2.bile.ibd)),sample_data(ps2.bile.ibd)[,c("Neoplasia", "Cluster")], bile.ibd.prop)
  }
}

bile.ibd.df.long <-  gather(bile.ibd.df, key = "Bile_acid", value = "Proportion", TUDCA:HCA)
bile.ibd.df.long$Bile_acid <- factor(bile.ibd.df.long$Bile_acid, 
                                     levels = c("CA", "CA.3S", "GCA", "TCA",
                                                "CDCA", "CDCA.3S", "GCDCA", "TCDCA",
                                                "HCA",
                                                "DCA", "DCA.3S", "GDCA", "TDCA",
                                                "LCA", "LCA.3S", "GLCA", "GLCA.3S", "TLCA", "TLCA.3S",
                                                "UDCA", "UDCA.3S", "GUDCA", "GUDCA.3S", "TUDCA", "TUDCA.3S",
                                                "HDCA", "THDCA"))
bile.ibd.df.long$SampleID <- factor(bile.ibd.df.long$SampleID, 
                                    levels = rev(rownames(bile.ibd.df)[order(bile.ibd.df$CA)]))
pal <- c("#FF3300", "#FF9966", "#CC3300", "#CC6666", "#990033", "#CC6600", "#FFCC33", "#FFFF66", "#CC9900",
         "#0000FF", "#99CCFF", "#003366", "#66FFFF", "#0099FF", "#33FF66", "#3399FF", "#006633", "#009966", 
         "#CCFF33", "#003300", "#9933FF", "#6600CC", "#CC33FF", "#9900CC", "#FF00FF", "#996699", "#9999CC")
p.bar <- ggplot(bile.ibd.df.long, aes(x=SampleID, y=Proportion, fill=Bile_acid)) + geom_bar(stat = "identity") + 
  facet_grid(~Cluster, scales = "free") + scale_fill_manual("Bile acids",values = pal) + xlab("") + 
  theme_classic() + 
  theme(text = element_text(size=20), axis.text.x=element_blank(), axis.text = element_text(size=16))

#plot
p.bar






#Concentration boxplots

#Cluster
#primary

bile.ibd.mod <- bile.ibd
bile.ibd.mod[,1:27] <- sweep(bile.ibd.mod[,1:27],1,bile.ibd.mod$Total_BA,"/")
#Normalise the grouped bile acids to total
bile.ibd.mod[,c(31:37,39:43,51:52)] <- sweep(bile.ibd.mod[,c(31:37,39:43,51:52)],1,bile.ibd.mod$Total_BA,"/")

bile.ibd.raw <- bile.ibd.mod[,c(32,33)]
colnames(bile.ibd.raw) <- gsub(".1", "", colnames(bile.ibd.raw))
if(!(is.null(rownames(bile.ibd))) & !(is.null(rownames(sample_data(ps2.bile.ibd))))){
  if(all(rownames(bile.ibd) == rownames(sample_data(ps2.bile.ibd)))){
    bile.ibd.df <- cbind(SampleID = rownames(sample_data(ps2.bile.ibd)),sample_data(ps2.bile.ibd)[,c("Neoplasia", "Cluster")], bile.ibd.raw)
  }
}

####### Figure 5H #######
bile.ibd.df.long <-  gather(bile.ibd.df, key = "Bile_acid", value = "Concentration", CA:CDCA)
bile.ibd.df.long$Bile_acid <- factor(bile.ibd.df.long$Bile_acid, 
                                     levels = c("CA", 
                                                "CDCA"))

pal <- c("red3", "grey25", "chartreuse2")
my_comparisons <- list(c("C.1", "C.3"), c("C.1", "C.2"), c("C.2", "C.3"))
p.clus <- ggplot(bile.ibd.df.long, aes(x=Cluster, y=Concentration, colour=Cluster)) + geom_boxplot(outlier.alpha = 0) +
  geom_point() +
  facet_grid(~Bile_acid, scales = "free") + scale_colour_manual("Cluster",values = pal) + xlab("") + 
  theme_classic() + 
  theme(text = element_text(size=16), axis.text.x=element_blank(), axis.text = element_text(size=12)) +
  stat_compare_means(comparisons = rev(my_comparisons), label = "p.signif") + ylab("Proportion of total")

#plot
p.clus

####### Figure 5I #######
#secondary
bile.ibd.raw <- bile.ibd.mod[,c(35,36,37)]
colnames(bile.ibd.raw) <- gsub(".1", "", colnames(bile.ibd.raw))
if(!(is.null(rownames(bile.ibd))) & !(is.null(rownames(sample_data(ps2.bile.ibd))))){
  if(all(rownames(bile.ibd) == rownames(sample_data(ps2.bile.ibd)))){
    bile.ibd.df <- cbind(SampleID = rownames(sample_data(ps2.bile.ibd)),sample_data(ps2.bile.ibd)[,c("Neoplasia", "Cluster")], bile.ibd.raw)
  }
}

bile.ibd.df.long <-  gather(bile.ibd.df, key = "Bile_acid", value = "Concentration", DCA:UDCA)
bile.ibd.df.long$Bile_acid <- factor(bile.ibd.df.long$Bile_acid, 
                                     levels = c("DCA",
                                                "LCA",
                                                "UDCA"))
pal <- c("red3", "grey25", "chartreuse2")
my_comparisons <- list(c("C.1", "C.3"), c("C.1", "C.2"), c("C.2", "C.3"))
p.clus <- ggplot(bile.ibd.df.long, aes(x=Cluster, y=Concentration, colour=Cluster)) + geom_boxplot(outlier.alpha = 0) +
  geom_point() +
  facet_grid(~Bile_acid, scales = "free") + scale_colour_manual("Cluster",values = pal) + xlab("") + 
  theme_classic() + 
  theme(text = element_text(size=16), axis.text.x=element_blank(), axis.text = element_text(size=12)) +
  stat_compare_means(comparisons = rev(my_comparisons), label = "p.signif") + ylab("Proportion of total")

pdf("Bile_acids/top_10_cluster_secondary.pdf", height = 3.2, width = 4.2)
p.clus
dev.off()


####### Figure supplementary #######
#conjugated


bile.ibd.raw <- bile.ibd.mod[,c(44:47, 53)]
colnames(bile.ibd.raw) <- c("UDCA", "Glycine-conjugated", "Taurine-conjugated", "Sulfo-conjugated", "Unconjugated")
if(!(is.null(rownames(bile.ibd))) & !(is.null(rownames(sample_data(ps2.bile.ibd))))){
  if(all(rownames(bile.ibd) == rownames(sample_data(ps2.bile.ibd)))){
    bile.ibd.df <- cbind(SampleID = rownames(sample_data(ps2.bile.ibd)),sample_data(ps2.bile.ibd)[,c("Neoplasia", "Cluster")], bile.ibd.raw)
  }
}

bile.ibd.df.long <-  gather(bile.ibd.df, key = "Bile_acid", value = "Concentration", 'UDCA':'Unconjugated')


pal <- c("red3", "grey25", "chartreuse2")
my_comparisons <- list(c("C.1", "C.3"), c("C.1", "C.2"), c("C.2", "C.3"))
p.clus <- ggplot(bile.ibd.df.long, aes(x=Cluster, y=Concentration, colour=Cluster)) + geom_boxplot(outlier.alpha = 0) +
  geom_point() +
  facet_grid(~Bile_acid, scales = "free") + scale_colour_manual("Cluster",values = pal) + xlab("") + 
  theme_classic() + 
  theme(text = element_text(size=16), axis.text.x=element_blank(), axis.text = element_text(size=12)) +
  stat_compare_means(comparisons = rev(my_comparisons), label = "p.signif") + ylab("Proportion of total")


p.clus



#Neoplasia

####### Figure 5J #######
#primary

bile.ibd.raw <- bile.ibd.mod[,c(32,33)]
colnames(bile.ibd.raw) <- gsub(".1", "", colnames(bile.ibd.raw))
if(!(is.null(rownames(bile.ibd))) & !(is.null(rownames(sample_data(ps2.bile.ibd))))){
  if(all(rownames(bile.ibd) == rownames(sample_data(ps2.bile.ibd)))){
    bile.ibd.df <- cbind(SampleID = rownames(sample_data(ps2.bile.ibd)),sample_data(ps2.bile.ibd)[,c("Neoplasia", "Cluster")], bile.ibd.raw)
  }
}

bile.ibd.df.long <-  gather(bile.ibd.df, key = "Bile_acid", value = "Concentration", CA:CDCA)
bile.ibd.df.long$Bile_acid <- factor(bile.ibd.df.long$Bile_acid, 
                                     levels = c("CA", 
                                                "CDCA"))

pal <- c("blue", "purple", "pink2", "salmon", "red")
my_comparisons <- list(c("N0", "N3"), c("N0", "N2"), c("N0", "N1"), c("N0", "A1"))
p.neo <- ggplot(bile.ibd.df.long, aes(x=Neoplasia, y=Concentration, colour=Neoplasia)) + geom_boxplot(outlier.alpha = 0) +
  geom_point() +
  facet_grid(~Bile_acid, scales = "free") + scale_colour_manual("Neoplasia",values = pal) + xlab("") + 
  theme_classic() + 
  theme(text = element_text(size=16), axis.text.x=element_blank(), axis.text = element_text(size=12)) +
  stat_compare_means(comparisons = rev(my_comparisons), label = "p.signif") + ylab("Proportion of total")

p.neo

####### Figure 5K #######
#secondary

bile.ibd.raw <- bile.ibd.mod[,c(35,36,37)]
colnames(bile.ibd.raw) <- gsub(".1", "", colnames(bile.ibd.raw))
if(!(is.null(rownames(bile.ibd))) & !(is.null(rownames(sample_data(ps2.bile.ibd))))){
  if(all(rownames(bile.ibd) == rownames(sample_data(ps2.bile.ibd)))){
    bile.ibd.df <- cbind(SampleID = rownames(sample_data(ps2.bile.ibd)),sample_data(ps2.bile.ibd)[,c("Neoplasia", "Cluster")], bile.ibd.raw)
  }
}

bile.ibd.df.long <-  gather(bile.ibd.df, key = "Bile_acid", value = "Concentration", DCA:UDCA)

pal <- c("blue", "purple", "pink2", "salmon", "red")
my_comparisons <- list(c("N0", "N3"), c("N0", "N2"), c("N0", "N1"), c("N0", "A1"))
p.neo <- ggplot(bile.ibd.df.long, aes(x=Neoplasia, y=Concentration, colour=Neoplasia)) + geom_boxplot(outlier.alpha = 0) +
  geom_point() +
  facet_grid(~Bile_acid, scales = "free") + scale_colour_manual("Neoplasia",values = pal) + xlab("") + 
  theme_classic() + 
  theme(text = element_text(size=16), axis.text.x=element_blank(), axis.text = element_text(size=12)) +
  stat_compare_means(comparisons = rev(my_comparisons), label = "p.signif") + ylab("Proportion of total")

p.neo


####### Supplementary figure #######
#conjugated

bile.ibd.raw <- bile.ibd.mod[,c(44:47, 53)] 
colnames(bile.ibd.raw) <- c("UDCA", "Glycine-conjugated", "Taurine-conjugated", "Sulfo-conjugated", "Unconjugated")
if(!(is.null(rownames(bile.ibd))) & !(is.null(rownames(sample_data(ps2.bile.ibd))))){
  if(all(rownames(bile.ibd) == rownames(sample_data(ps2.bile.ibd)))){
    bile.ibd.df <- cbind(SampleID = rownames(sample_data(ps2.bile.ibd)),sample_data(ps2.bile.ibd)[,c("Neoplasia", "Cluster")], bile.ibd.raw)
  }
}

bile.ibd.df.long <-  gather(bile.ibd.df, key = "Bile_acid", value = "Concentration", 'UDCA':'Unconjugated')


pal <- c("blue", "purple", "pink2", "salmon", "red")
my_comparisons <- list(c("N0", "N3"), c("N0", "N2"), c("N0", "N1"), c("N0", "A1"))
p.clus <- ggplot(bile.ibd.df.long, aes(x=Neoplasia, y=Concentration, colour=Neoplasia)) + geom_boxplot(outlier.alpha = 0) +
  geom_point() +
  facet_grid(~Bile_acid, scales = "free") + scale_colour_manual("Neoplasia",values = pal) + xlab("") + 
  theme_classic() + 
  theme(text = element_text(size=16), axis.text.x=element_blank(), axis.text = element_text(size=12)) +
  stat_compare_means(comparisons = rev(my_comparisons), label = "p.signif") + ylab("Proportion of total")

p.clus






