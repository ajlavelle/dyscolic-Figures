### Import data from Qiime2 pipeline ###


#Step 1 - load phyloseq and set working directory
library(phyloseq)
library(ggplot2)
library(ggpubr)

#Step 2 - make sure file in working directory using setwd(~/...) and give them names in R. The 'exported' folder should be in this directory
list.files() #lists the files in the current working directory

biomfilename  = "exported/otu_table_json.biom" #otu_table_json.biom should be the name of the file stored in the working directory
treefilename = "exported/tree.nwk";

#Step 3 - use import_biome function
ps <- import_biom(biomfilename, treefilename)
#change orientation of otu table
otu_table(ps) <- t(otu_table(ps))

#Step 4 - parse input
taxtable <- read.table("exported/taxonomy.tsv", header=T, sep = "\t") #Qiime2
if(dim(taxtable)[1] != dim(otu_table(ps))[2]){
  print("Error in importing! - Check for  ' in taxa names")
} else {
  taxtable$Taxon <- as.character(taxtable$Taxon)
  tax.mat <- matrix(ncol = 7, nrow = nrow(taxtable))
  for(i in 1:nrow(tax.mat)){
    ss.x <- strsplit(taxtable$Taxon[i], ";")
    for(j in 1:ncol(tax.mat)){
      if(is.na(ss.x[[1]][j]) == FALSE){
        tax.mat[i,j] <- ss.x[[1]][j]
      } else {
        tax.mat[i,j] <- NA
      }
    }
  }
  
  colnames(tax.mat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  rownames(tax.mat) <- taxtable$Feature.ID
  tax.mat[,1] <- gsub("D_0__", "", tax.mat[,1])
  tax.mat[,2] <- gsub("D_1__", "", tax.mat[,2])
  tax.mat[,3] <- gsub("D_2__", "", tax.mat[,3])
  tax.mat[,4] <- gsub("D_3__", "", tax.mat[,4])
  tax.mat[,5] <- gsub("D_4__", "", tax.mat[,5])
  tax.mat[,6] <- gsub("D_5__", "", tax.mat[,6])
  tax.mat[,7] <- gsub("D_6__", "", tax.mat[,7])
  
  tax.mat[,1] <- gsub(" ", "_", tax.mat[,1])
  tax.mat[,2] <- gsub(" ", "_", tax.mat[,2])
  tax.mat[,3] <- gsub(" ", "_", tax.mat[,3])
  tax.mat[,4] <- gsub(" ", "_", tax.mat[,4])
  tax.mat[,5] <- gsub(" ", "_", tax.mat[,5])
  tax.mat[,6] <- gsub(" ", "_", tax.mat[,6])
  tax.mat[,7] <- gsub(" ", "_", tax.mat[,7])
  
  tax.mat[,1] <- gsub("-", "_", tax.mat[,1])
  tax.mat[,2] <- gsub("-", "_", tax.mat[,2])
  tax.mat[,3] <- gsub("-", "_", tax.mat[,3])
  tax.mat[,4] <- gsub("-", "_", tax.mat[,4])
  tax.mat[,5] <- gsub("-", "_", tax.mat[,5])
  tax.mat[,6] <- gsub("-", "_", tax.mat[,6])
  tax.mat[,7] <- gsub("-", "_", tax.mat[,7])
  
  tax.mat[,1] <- gsub("\\[", "", tax.mat[,1])
  tax.mat[,2] <- gsub("\\[", "", tax.mat[,2])
  tax.mat[,3] <- gsub("\\[", "", tax.mat[,3])
  tax.mat[,4] <- gsub("\\[", "", tax.mat[,4])
  tax.mat[,5] <- gsub("\\[", "", tax.mat[,5])
  tax.mat[,6] <- gsub("\\[", "", tax.mat[,6])
  tax.mat[,7] <- gsub("\\[", "", tax.mat[,7])
  
  tax.mat[,1] <- gsub("\\]", "", tax.mat[,1])
  tax.mat[,2] <- gsub("\\]", "", tax.mat[,2])
  tax.mat[,3] <- gsub("\\]", "", tax.mat[,3])
  tax.mat[,4] <- gsub("\\]", "", tax.mat[,4])
  tax.mat[,5] <- gsub("\\]", "", tax.mat[,5])
  tax.mat[,6] <- gsub("\\]", "", tax.mat[,6])
  tax.mat[,7] <- gsub("\\]", "", tax.mat[,7])
  
  tax.mat[,7][grepl("uncultured", tax.mat[,7])] <- NA
  tax.mat[,7][grepl("Ambiguous", tax.mat[,7])] <- NA
  tax.mat[,7][grepl("unident", tax.mat[,7])] <- NA
  
}

#match tax table and otu table and ensure match
tax.mat <- tax.mat[match(colnames(otu_table(ps)), rownames(tax.mat)),]
all(colnames(otu_table(ps)) == rownames(tax.mat))

#import mapping
mapping_file <- read.table("exported/mapping_r.txt", header=T)

#match mapping file and otu table and ensure amtch
mapping_file <- mapping_file[match(rownames(otu_table(ps)), rownames(mapping_file)),]
all(rownames(mapping_file) == rownames(otu_table(ps)))

#add to phyloseq object
sample_data(ps) <- sample_data(mapping_file)
tax_table(ps) <- tax_table(tax.mat)


########################################################################

                      ###### Filtering #######

########################################################################

#keep unfiltered as ps
ps_new <- ps


# Create table, number of features for each phyla and remove NA phyla
table(tax_table(ps_new)[, "Phylum"], exclude = NULL)
ps_new <- subset_taxa(ps_new, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))


# Define prevalence of each taxa
# (in how many samples did each taxa appear at least once)
prev0 = apply(X = otu_table(ps_new),
              MARGIN = ifelse(taxa_are_rows(ps_new), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prev0,
                    TotalAbundance = taxa_sums(ps_new),
                    tax_table(ps_new))
keepPhyla = table(prevdf$Phylum)[(table(prevdf$Phylum) > 0)]
prevdf1 = subset(prevdf, Phylum %in% names(keepPhyla))
# Define prevalence threshold as 1 
prevalenceThreshold = 1
# Execute prevalence filter, using `prune_taxa()` function
ps1 = prune_taxa((prev0 > prevalenceThreshold), ps_new)
ps1
#ps2
ps2 <- ps1


# check retained seqs post filtering
sum(colSums(otu_table(ps2)))/sum(colSums(otu_table(ps)))

# check range of sequences per sample
rowSums(otu_table(ps2))[order(rowSums(otu_table(ps2)))]



#################################################################

                  ### Setting factors ###

#################################################################

#set factors
sample_data(ps2)$Subtype <- factor(sample_data(ps2)$Subtype, levels = c("Control", "CD", "UC"))
#set factors
sample_data(ps2)$Neoplasia <- factor(sample_data(ps2)$Neoplasia, levels = c("N0", "A1", "N1", "N2", "N3"))

#Combine neoplasia
sample_data(ps2)$Neoplasia2 <- as.character(sample_data(ps2)$Neoplasia)
sample_data(ps2)$Neoplasia2[sample_data(ps2)$Neoplasia == "N1" | 
                              sample_data(ps2)$Neoplasia == "N2" | 
                              sample_data(ps2)$Neoplasia == "N3"] <- "Nx"
sample_data(ps2)$Neoplasia2 <- factor(sample_data(ps2)$Neoplasia2, levels = c("N0", "Nx", "A1"))

#for differential abundance testing on each IBD subtype
sample_data(ps2)$Neoplasia.st <- paste(as.character(sample_data(ps2)$Neoplasia2),as.character(sample_data(ps2)$Subtype), sep = ".")
sample_data(ps2)$Neoplasia.st <- factor(sample_data(ps2)$Neoplasia.st, 
                                        levels = c("N0.Control", "A1.Control", "Nx.Control",
                                                   "N0.UC", "Nx.UC", "A1.UC", "N0.CD", "Nx.CD", "A1.CD"))
#for differential abundance testing on IBD
sample_data(ps2)$Neoplasia.st.ibd <- paste(as.character(sample_data(ps2)$Neoplasia2),as.character(sample_data(ps2)$Subtype), sep = ".")
sample_data(ps2)$Neoplasia.st.ibd[sample_data(ps2)$Neoplasia.st.ibd == "N0.UC"] <- "N0.IBD"
sample_data(ps2)$Neoplasia.st.ibd[sample_data(ps2)$Neoplasia.st.ibd == "N0.CD"] <- "N0.IBD"
sample_data(ps2)$Neoplasia.st.ibd[sample_data(ps2)$Neoplasia.st.ibd == "Nx.UC"] <- "Nx.IBD"
sample_data(ps2)$Neoplasia.st.ibd[sample_data(ps2)$Neoplasia.st.ibd == "Nx.CD"] <- "Nx.IBD"
sample_data(ps2)$Neoplasia.st.ibd[sample_data(ps2)$Neoplasia.st.ibd == "A1.UC"] <- "A1.IBD"
sample_data(ps2)$Neoplasia.st.ibd[sample_data(ps2)$Neoplasia.st.ibd == "A1.CD"] <- "A1.IBD"

sample_data(ps2)$Neoplasia.st.ibd <- factor(sample_data(ps2)$Neoplasia.st.ibd, 
                                            levels = c("N0.Control", "A1.Control", "Nx.Control",
                                                       "N0.IBD", "Nx.IBD", "A1.IBD"))

#for plotting diversity
sample_data(ps2)$Neoplasia3cat <- as.character(sample_data(ps2)$Neoplasia)
sample_data(ps2)$Neoplasia3cat[sample_data(ps2)$Neoplasia3cat == "N1"] <- "N1-N3"
sample_data(ps2)$Neoplasia3cat[sample_data(ps2)$Neoplasia3cat == "N2"] <- "N1-N3"
sample_data(ps2)$Neoplasia3cat[sample_data(ps2)$Neoplasia3cat == "N3"] <- "N1-N3"
sample_data(ps2)$Neoplasia3cat <- factor(sample_data(ps2)$Neoplasia3cat, levels = c("N0", "A1", "N1-N3"))

######### subset

ps2.ibd <- subset_samples(ps2, sample_data(ps2)$Subtype != "Control")
clusters <- read.table("Bile_acids/clusters.txt", header = T)

if(all(rownames(clusters) == rownames(sample_data(ps2.ibd)))){
  sample_data(ps2.ibd)$Cluster <- clusters$DMM_run_1
} else {
  print("ERROR! Rownames do not match!")
}


#################################################################

                  ### Alpha diversity ###

#################################################################



######alpha diversity#####
library(RColorBrewer)
library(ggpubr)
dir.create(paste(getwd(), "alpha_diversity", sep = "/"))




####### Figure 2A #######

# All combined
#set variable, palate and comparisons
var <- "Subtype"; pal <- c("dodgerblue", "orange", "coral3"); my_comparisons <- list(c("Control", "UC"), c("Control", "CD"), c("CD", "UC"))

set.seed(101)
ps_rare <- ps2
p <- plot_richness(ps_rare, x=var, color = var, measures = c("Shannon")) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_color_discrete(var) + xlab("") + ylab("Shannon Diversity") + geom_point()
p2 <- p + geom_point(size=2) + 
  theme(strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme_classic() + 
  theme(text = element_text(size=22), axis.text = element_text(size=16) ,legend.position = "right") + 
  scale_color_manual(var, values = pal)

#plot
p2 + theme(axis.text.x=element_blank()) + stat_compare_means(comparisons = rev(my_comparisons),label = "p.signif")



#stratified by neoplasia
#set variable, palate and comparisons
var <- "Subtype"; var1 <- "Subtype"; pal <- c("dodgerblue", "orange", "coral3"); my_comparisons <- list(c("Control", "UC"), c("Control", "CD"), c("CD", "UC"))

set.seed(101)
ps_rare <- ps2
p <- plot_richness(ps_rare, x=var, color = var, measures = c("Shannon")) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_color_discrete(var) + xlab("") + ylab("") + geom_point() + facet_grid(~Neoplasia3cat, scales = "free_y")
p2 <- p + geom_point(size=2) + 
  theme(strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme_classic() + 
  theme(text = element_text(size=22), axis.text = element_text(size=16) ,legend.position = "right") + 
  scale_color_manual(var1, values = pal)

#plot
p2 + theme(axis.text.x=element_blank()) + stat_compare_means(comparisons = rev(my_comparisons),label = "p.signif")






####### Figure 2B #######

# All combined
#set variable, palate and comparisons
var <- "Neoplasia3cat"; var1 <- "Neoplasia"; pal <- c("blue", "purple", "red3"); my_comparisons <- list(c("N0", "N1-N3"), c("N0", "A1"), c("A1", "N1-N3"))

set.seed(101)
ps_rare <- ps2
p <- plot_richness(ps_rare, x=var, color = var, measures = c("Shannon")) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_color_discrete("Neoplasia",var) + xlab("") + ylab("Shannon Diversity") + geom_point()
p2 <- p + geom_point(size=2) + 
  theme(strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme_classic() + 
  theme(text = element_text(size=22), axis.text = element_text(size=16) ,legend.position = "right") + 
  scale_color_manual(var1, values = pal)

#plot
p2 + theme(axis.text.x=element_blank()) + stat_compare_means(comparisons = rev(my_comparisons),label = "p.signif")

#3 cat-stratified by subtype
var <- "Neoplasia3cat"; var1 <- "Neoplasia"; pal <- c("blue", "purple", "red3"); my_comparisons <- list(c("N0", "N1-N3"), c("N0", "A1"), c("A1", "N1-N3"))

set.seed(101)
ps_rare <- ps2
p <- plot_richness(ps_rare, x=var, color = var, measures = c("Shannon")) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_color_discrete("Neoplasia",var) + xlab("") + ylab("") + geom_point() + facet_grid(~Subtype, scales = "free_y")
p2 <- p + geom_point(size=2) + 
  theme(strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme_classic() + 
  theme(text = element_text(size=22), axis.text = element_text(size=16) ,legend.position = "right") + 
  scale_color_manual(var1, values = pal)

#plot
p2 + theme(axis.text.x=element_blank()) + stat_compare_means(comparisons = rev(my_comparisons),label = "p.signif")






####### Figure 3D #######

#set variable, palate and comparisons
var <- "Cluster"; pal <- c("red3", "grey25", "chartreuse2"); my_comparisons <- list(c("C.1", "C.2"), c("C.2", "C.3"), c("C.1", "C.3"))

set.seed(101)
ps_rare <- ps2.ibd
p <- plot_richness(ps_rare, x=var, color = var, measures = c("Observed", "Shannon")) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_color_discrete(var) + xlab("") + ylab("") + geom_point()
p2 <- p + geom_point(size=2) + 
  theme(strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme_classic() + 
  theme(text = element_text(size=22), axis.text = element_text(size=16) ,legend.position = "right") + 
  scale_color_manual(var, values = pal)

#plot
p2 + theme(axis.text.x=element_blank()) + stat_compare_means(comparisons = my_comparisons,label = "p.signif")








#################################################################

                  ### Beta diversity ###

#################################################################


#TSS normalise data
ps.transf <- transform_sample_counts(ps2, function(x) {x/sum(x)})

#ordination Bray-Curtis/MDS
out.donor.bray <- ordinate(ps.transf, method = "MDS", distance = "bray")

####### Figure 2C #######

#set variable, legend name and palate
var <-"Subtype"; var1 <- "Subtype"; pal <- c("dodgerblue", "orange", "coral3"); 

bray.df <- cbind(out.donor.bray$vectors[,1:2], sample_data(ps.transf))
eval.1 <- (out.donor.bray$values$Eigenvalues[1]/sum(out.donor.bray$values$Eigenvalues))*100
eval.2 <- (out.donor.bray$values$Eigenvalues[2]/sum(out.donor.bray$values$Eigenvalues))*100
p.bray <- ggplot(bray.df, aes(x=Axis.1, y=Axis.2, colour=Subtype, fill=Subtype)) + 
  geom_point(size=2, alpha=1) + scale_fill_manual(values = pal) + 
  scale_colour_manual(values = pal) + stat_ellipse(level = 0.8) + theme_classic() + 
  xlab(paste("PCoA 1 [", round(eval.1, 2), "%]", sep = "")) + 
  ylab(paste("PCoA 2 [", round(eval.2, 2), "%]", sep = "")) + 
  theme(legend.position = "right", text = element_text(size=16), axis.ticks.length = unit(0.2, "cm"), 
        axis.text = element_text(size=14), plot.title = element_text(hjust = 0.5), 
        axis.line = element_line(colour = 'black', size = 1.05), 
        axis.ticks = element_line(colour = 'black', size = 1.05)) +
  ggtitle("Bray-Curtis divergence")



#now check PERMANOVA and add resulting values into plot
set.seed(101)
metadata <- as(sample_data(ps.transf), "data.frame")
x <- adonis(otu_table(ps.transf) ~ unlist(sample_data(ps.transf)[,var]), permutations = 999, method = "bray", data=metadata)

#calculate annotation co-ordinates
x.min <- ggplot_build(p.bray)$layout$panel_scales_x[[1]]$range$range[1]
x.max <- ggplot_build(p.bray)$layout$panel_scales_x[[1]]$range$range[2]
x.range <- x.max -x.min
y.min <- ggplot_build(p.bray)$layout$panel_scales_y[[1]]$range$range[1]
y.max <- ggplot_build(p.bray)$layout$panel_scales_y[[1]]$range$range[2]
y.range <- y.max -y.min

#add results to plot
p.bray2 <- p.bray + annotate("label", label = sprintf("italic(R)^2 == %s", round(x$aov.tab$R2[1],3)), x = x.min + x.range/5, y = y.max*1.3, parse = TRUE, label.size = NA, size = 4) + 
  annotate("label", label = sprintf("italic(P) == %s", round(x$aov.tab$`Pr(>F)`[1], 5)), x = x.min + x.range/5, y = y.max*1.1, parse = TRUE, label.size = NA, size = 4)

#subtype plot
p.bray2


####### Figure 2D #######

#set variable, legend name and palate
var <- "Neoplasia"; pal <- c("blue", "purple", "pink2", "salmon", "red");


p.bray <- ggplot(bray.df, aes(x=Axis.1, y=Axis.2, colour=Neoplasia, fill=Neoplasia)) + 
  geom_point(size=2, alpha=1) + scale_fill_manual(values = pal) + 
  scale_colour_manual(values = pal) + stat_ellipse(level = 0.8) + theme_classic() + 
  xlab(paste("PCoA 1 [", round(eval.1, 2), "%]", sep = "")) + 
  ylab(paste("PCoA 2 [", round(eval.2, 2), "%]", sep = "")) + 
  theme(legend.position = "right", text = element_text(size=16), axis.ticks.length = unit(0.2, "cm"), 
        axis.text = element_text(size=14), plot.title = element_text(hjust = 0.5), 
        axis.line = element_line(colour = 'black', size = 1.05), 
        axis.ticks = element_line(colour = 'black', size = 1.05)) +
  ggtitle("Bray-Curtis divergence")

#now check PERMANOVA and add resulting values into plot
set.seed(101)
metadata <- as(sample_data(ps.transf), "data.frame")
x <- adonis(otu_table(ps.transf) ~ unlist(sample_data(ps.transf)[,var]), permutations = 999, method = "bray", data=metadata)

#calculate annotation co-ordinates
x.min <- ggplot_build(p.bray)$layout$panel_scales_x[[1]]$range$range[1]
x.max <- ggplot_build(p.bray)$layout$panel_scales_x[[1]]$range$range[2]
x.range <- x.max -x.min
y.min <- ggplot_build(p.bray)$layout$panel_scales_y[[1]]$range$range[1]
y.max <- ggplot_build(p.bray)$layout$panel_scales_y[[1]]$range$range[2]
y.range <- y.max -y.min

#add results to plot
p.bray2 <- p.bray + annotate("label", label = sprintf("italic(R)^2 == %s", round(x$aov.tab$R2[1],3)), x = x.min + x.range/5, y = y.max*1.3, parse = TRUE, label.size = NA, size = 4) + 
  annotate("label", label = sprintf("italic(P) == %s", round(x$aov.tab$`Pr(>F)`[1], 5)), x = x.min + x.range/5, y = y.max*1.1, parse = TRUE, label.size = NA, size = 4)

#plot
p.bray2



####### Figure 3A #######

#Use IBD only

#TSS normalise data
ps.transf <- transform_sample_counts(ps2.ibd, function(x) {x/sum(x)})

#ordination Bray-Curtis/MDS
out.donor.bray <- ordinate(ps.transf, method = "MDS", distance = "bray")
bray.df <- cbind(out.donor.bray$vectors[,1:2], sample_data(ps.transf))
eval.1 <- (out.donor.bray$values$Eigenvalues[1]/sum(out.donor.bray$values$Eigenvalues))*100
eval.2 <- (out.donor.bray$values$Eigenvalues[2]/sum(out.donor.bray$values$Eigenvalues))*100

#set variable, legend name and palate
var <- "Neoplasia"; pal <- c("blue", "purple", "pink2", "salmon", "red");


p.bray <- ggplot(bray.df, aes(x=Axis.1, y=Axis.2, colour=Neoplasia, fill=Neoplasia)) + 
  geom_point(size=2, alpha=1) + scale_fill_manual(values = pal) + 
  scale_colour_manual(values = pal) + stat_ellipse(level = 0.8) + theme_classic() + 
  xlab(paste("PCoA 1 [", round(eval.1, 2), "%]", sep = "")) + 
  ylab(paste("PCoA 2 [", round(eval.2, 2), "%]", sep = "")) + 
  theme(legend.position = "right", text = element_text(size=16), axis.ticks.length = unit(0.2, "cm"), 
        axis.text = element_text(size=14), plot.title = element_text(hjust = 0.5), 
        axis.line = element_line(colour = 'black', size = 1.05), 
        axis.ticks = element_line(colour = 'black', size = 1.05)) +
  ggtitle("Bray-Curtis divergence")

#now check PERMANOVA and add resulting values into plot
set.seed(101)
metadata <- as(sample_data(ps.transf), "data.frame")
x <- adonis(otu_table(ps.transf) ~ unlist(sample_data(ps.transf)[,var]), permutations = 999, method = "bray", data=metadata)

#calculate annotation co-ordinates
x.min <- ggplot_build(p.bray)$layout$panel_scales_x[[1]]$range$range[1]
x.max <- ggplot_build(p.bray)$layout$panel_scales_x[[1]]$range$range[2]
x.range <- x.max -x.min
y.min <- ggplot_build(p.bray)$layout$panel_scales_y[[1]]$range$range[1]
y.max <- ggplot_build(p.bray)$layout$panel_scales_y[[1]]$range$range[2]
y.range <- y.max -y.min

#add results to plot
p.bray2 <- p.bray + annotate("label", label = sprintf("italic(R)^2 == %s", round(x$aov.tab$R2[1],3)), x = x.min + x.range/5, y = y.max*1.3, parse = TRUE, label.size = NA, size = 4) + 
  annotate("label", label = sprintf("italic(P) == %s", round(x$aov.tab$`Pr(>F)`[1], 5)), x = x.min + x.range/5, y = y.max*1.1, parse = TRUE, label.size = NA, size = 4)

#plot
p.bray2


####### Figure 3B #######

#set variable, legend name and palate
var <- "Cluster"; pal <- c("red3", "grey25", "chartreuse2");


p.bray <- ggplot(bray.df, aes(x=Axis.1, y=Axis.2, colour=Cluster, fill=Cluster)) + 
  geom_point(size=2, alpha=1) + scale_fill_manual(values = pal) + 
  scale_colour_manual(values = pal) + stat_ellipse(level = 0.8) + theme_classic() + 
  xlab(paste("PCoA 1 [", round(eval.1, 2), "%]", sep = "")) + 
  ylab(paste("PCoA 2 [", round(eval.2, 2), "%]", sep = "")) + 
  theme(legend.position = "right", text = element_text(size=16), axis.ticks.length = unit(0.2, "cm"), 
        axis.text = element_text(size=14), plot.title = element_text(hjust = 0.5), 
        axis.line = element_line(colour = 'black', size = 1.05), 
        axis.ticks = element_line(colour = 'black', size = 1.05)) +
  ggtitle("Bray-Curtis divergence")

#now check PERMANOVA and add resulting values into plot
set.seed(101)
metadata <- as(sample_data(ps.transf), "data.frame")
x <- adonis(otu_table(ps.transf) ~ unlist(sample_data(ps.transf)[,var]), permutations = 999, method = "bray", data=metadata)

#calculate annotation co-ordinates
x.min <- ggplot_build(p.bray)$layout$panel_scales_x[[1]]$range$range[1]
x.max <- ggplot_build(p.bray)$layout$panel_scales_x[[1]]$range$range[2]
x.range <- x.max -x.min
y.min <- ggplot_build(p.bray)$layout$panel_scales_y[[1]]$range$range[1]
y.max <- ggplot_build(p.bray)$layout$panel_scales_y[[1]]$range$range[2]
y.range <- y.max -y.min

#add results to plot
p.bray2 <- p.bray + annotate("label", label = sprintf("italic(R)^2 == %s", round(x$aov.tab$R2[1],3)), x = x.min + x.range/5, y = y.max*1.3, parse = TRUE, label.size = NA, size = 4) + 
  annotate("label", label = sprintf("italic(P) == %s", round(x$aov.tab$`Pr(>F)`[1], 5)), x = x.min + x.range/5, y = y.max*1.1, parse = TRUE, label.size = NA, size = 4)

#plot
p.bray2




####### Figure 3C #######
library(ggrepel)



#set variables
var <- "Cluster"; var1 <- "Cluster"
pal <- c("red3", "grey25", "chartreuse2")


out.donor.bray <- ordinate(ps.transf, method = "MDS", distance = "bray")
bray.df <- cbind(out.donor.bray$vectors[,1:2], sample_data(ps.transf))

eval.1 <- (out.donor.bray$values$Eigenvalues[1]/sum(out.donor.bray$values$Eigenvalues))*100
eval.2 <- (out.donor.bray$values$Eigenvalues[2]/sum(out.donor.bray$values$Eigenvalues))*100
p2 <- ggplot(bray.df, aes(x=Axis.1, y=Axis.2, colour=Cluster, fill=Cluster)) + geom_point(size=2) + scale_fill_manual("Cluster", values = pal) + 
  scale_colour_manual("Cluster", values = pal) + stat_ellipse(level = 0.8) + theme_classic() + 
  xlab(paste("PCoA 1 [", round(eval.1, 2), "%]", sep = "")) + ylab(paste("PCoA 2 [", round(eval.2, 2), "%]", sep = "")) + 
  theme(legend.position = "right", text = element_text(size=16), axis.ticks.length = unit(0.2, "cm"), 
        axis.text = element_text(size=14), plot.title = element_text(hjust = 0.5), 
        axis.line = element_line(colour = 'black', size = 1.05), axis.ticks = element_line(colour = 'black', size = 1.05))

ps2.g.ibd <- tax_glom(ps2.ibd, "Genus", NArm = FALSE)

ps2.g.ibd.prop <- transform_sample_counts(ps2.g.ibd, function(x) x/sum(x))
biplot_otu_table_16S = as(t(otu_table(ps2.g.ibd.prop)), "matrix")
all(rownames(biplot_otu_table_16S) == rownames(tax_table(ps2.g.ibd.prop)))
otu.vec <- c()

for(i in 1:nrow(tax_table(ps2.g.ibd.prop))){
  if(is.na(tax_table(ps2.g.ibd.prop)[i,6]) == FALSE  & grepl("uncultur", ignore.case = TRUE, tax_table(ps2.g.ibd.prop)[i,6]) == FALSE){
    otu.vec <- c(otu.vec, paste("g", tax_table(ps2.g.ibd.prop)[i,6], sep = "_"))
  } else if(is.na(tax_table(ps2.g.ibd.prop)[i,5]) == FALSE  & grepl("uncultur", ignore.case = TRUE, tax_table(ps2.g.ibd.prop)[i,5]) == FALSE){
    otu.vec <- c(otu.vec, paste("f", tax_table(ps2.g.ibd.prop)[i,5], sep = "_"))
  } else if(is.na(tax_table(ps2.g.ibd.prop)[i,4]) == FALSE  & grepl("uncultur", ignore.case = TRUE, tax_table(ps2.g.ibd.prop)[i,4]) == FALSE){
    otu.vec <- c(otu.vec, paste("o", tax_table(ps2.g.ibd.prop)[i,4], sep = "_"))
  } else{
    otu.vec <- c(otu.vec, paste("c", tax_table(ps2.g.ibd.prop)[i,3], sep = "_"))
  }
}

rownames(biplot_otu_table_16S) <- otu.vec

biplot_pcoa_16S = t(out.donor.bray$vectors[,1:2])[,colnames(biplot_otu_table_16S)]

axis1.cor = apply(biplot_otu_table_16S, 1, function(x) cor.test(x, biplot_pcoa_16S[1,], method = "spearman")$estimate)
axis2.cor = apply(biplot_otu_table_16S, 1, function(x) cor.test(x, biplot_pcoa_16S[2,], method = "spearman")$estimate)
axis1.p = apply(biplot_otu_table_16S, 1, function(x) cor.test(x, biplot_pcoa_16S[1,], method = "spearman")$p.value)
axis2.p = apply(biplot_otu_table_16S, 1, function(x) cor.test(x, biplot_pcoa_16S[2,], method = "spearman")$p.value)
axis1.p = p.adjust(axis1.p, method = "fdr")
axis2.p = p.adjust(axis2.p, method = "fdr")
axis1.cor[axis1.p > 0.05] = 0 # zero insignificant correlations
axis2.cor[axis2.p > 0.05] = 0
axis1.cor = axis1.cor[names((sort(abs(axis1.cor), decreasing = TRUE)))] # order by absolute value
axis2.cor = axis2.cor[names((sort(abs(axis2.cor), decreasing = TRUE)))]

top_corr_otu = unique(c(names(axis1.cor)[1:26], names(axis2.cor)[1:26]))
top_corr_otu <- top_corr_otu[!(is.na(top_corr_otu))]
vectors_ax1_end = axis1.cor[top_corr_otu] 
vectors_ax2_end = axis2.cor[top_corr_otu]
axis1.p.top = axis1.p[top_corr_otu] 
axis2.p.top = axis2.p[top_corr_otu]
vectors_df = data.frame(top_corr_otu, vectors_ax1_end, vectors_ax2_end, axis1.p.top, axis2.p.top)
vectors_df$abundance = apply(biplot_otu_table_16S[top_corr_otu,], 1, sum)
vectors_df <- vectors_df[!(vectors_df$vectors_ax1_end == 0 & vectors_df$vectors_ax2_end == 0),]
vectors_df

#plot
p2 + geom_segment(data = vectors_df, inherit.aes = FALSE, show.legend = FALSE, arrow = arrow(angle = 15, type = "closed"),
                  aes(x = 0, y = 0, xend = vectors_ax1_end, yend = vectors_ax2_end)) + 
  geom_text_repel(data = vectors_df, inherit.aes = FALSE, nudge_x = -0.1, nudge_y = 0.1, 
                  aes(vectors_ax1_end, vectors_ax2_end, label = top_corr_otu), size = 5) + 
  ggtitle("Correlation of genera with PCoA axes") + theme(text = element_text(size=18), legend.position = "right")







#################################################################

                        ### DESeq2 ###

#################################################################


library(DESeq2)
library(RColorBrewer)
# 1. Assign to ps.x
ps.x <- ps2
ps.x <- tax_glom(ps.x, "Genus", NArm = FALSE)


# 2. Make a colour vector to keep phylum colour assignments consistent between plots
dd <- unique(tax_table(ps.x)[,2])
dd.col <- c(brewer.pal(8, "Set2"), "#000000", brewer.pal(9, "Set1"))
dd.col <- c(brewer.pal(8, "Set1")[2], brewer.pal(8, "Set1")[4], brewer.pal(8, "Set1")[5], brewer.pal(8, "Set1")[1], "#000000", brewer.pal(9, "Set1"))

names(dd.col)  <- dd
dd.col <- dd.col[1:12]

# 3. Run DESeq2
deseq2 = phyloseq_to_deseq2(ps.x, ~ Neoplasia.st) #manually add this
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
deseq2_geoMeans = apply(counts(deseq2), 1, gm_mean)
deseq2 = estimateSizeFactors(deseq2, geoMeans = deseq2_geoMeans)
deseq2 = DESeq(deseq2, fitType="local")



####### Figure 2F #######
### N0.UC versus Nx.UC

# 4. Select significance threshold for FDR threshold
alpha.var <- 0.2

# 5. Extract relevant results (N0 versus Nx)
res = results(deseq2, contrast=c("Neoplasia.st", "N0.UC", "Nx.UC"))

# 6. Order results at chosen level of (FDR-adjusted) significance
res = res[order(res$padj, na.last=NA), ]
sigtab = res[(res$padj < alpha.var), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.x)[rownames(sigtab), ], "matrix"))
sigtab = sigtab[order(sigtab$log2FoldChange, decreasing = FALSE),]
sigtab$Genus <- as.character(sigtab$Genus)
sigtab$Genus[is.na(sigtab$Genus)] <- as.character(sigtab$Family[is.na(sigtab$Genus)])
sigtab$Genus = factor(sigtab$Genus, levels = unique(sigtab$Genus))

#only keep significant by Wilcoxon test
ps.y <- transform_sample_counts(ps.x, function(x) {x/sum(x)})
ps.y <- subset_taxa(ps.y, rownames(tax_table(ps.x)) %in% rownames(sigtab))
otus <- as.data.frame(as.matrix(otu_table(ps.y)))
otus <- otus[sample_data(ps.y)$Neoplasia.st %in% c("N0.UC", "Nx.UC"),]
samp <- sample_data(ps.y)[sample_data(ps.y)$Neoplasia.st %in% c("N0.UC", "Nx.UC"),]

sig_by_Wilcox <- colnames(otus)[apply(otus, 2, function(x) wilcox.test(x ~ samp$Neoplasia.st.ibd)$p.value) < 0.05]
sigtab <- sigtab[rownames(sigtab) %in% sig_by_Wilcox,]

ggplot(sigtab, aes(y = Genus, x = log2FoldChange, color = Phylum)) + geom_point(size=3) + theme_classic() + 
  theme(text = element_text(size=22), legend.text = element_text(size=14)) + scale_color_manual(values = dd.col) +
  geom_vline(xintercept = 0, colour = "Red", lty = 2) + theme(legend.position="right") + 
  ggtitle(paste("UC: No dysplasia\n", "versus\n", "Dysplasia (any grade)", sep = ""))



### N0.CD versus Nx.CD - no plot

# 4. Select significance threshold for FDR threshold
alpha.var <- 0.2

# 5. Extract relevant results (N0 versus Nx)
res = results(deseq2, contrast=c("Neoplasia.st", "N0.CD", "Nx.CD"))

# 6. Order results at chosen level of (FDR-adjusted) significance
res = res[order(res$padj, na.last=NA), ]
sigtab = res[(res$padj < alpha.var), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.x)[rownames(sigtab), ], "matrix"))
sigtab = sigtab[order(sigtab$log2FoldChange, decreasing = FALSE),]
sigtab$Genus <- as.character(sigtab$Genus)
sigtab$Genus[is.na(sigtab$Genus)] <- as.character(sigtab$Family[is.na(sigtab$Genus)])
sigtab$Genus = factor(sigtab$Genus, levels = unique(sigtab$Genus))

#only keep significant by Wilcoxon test
ps.y <- transform_sample_counts(ps.x, function(x) {x/sum(x)})
ps.y <- subset_taxa(ps.y, rownames(tax_table(ps.x)) %in% rownames(sigtab))
otus <- as.data.frame(as.matrix(otu_table(ps.y)))
otus <- otus[sample_data(ps.y)$Neoplasia.st%in% c("N0.CD", "Nx.CD"),]
samp <- sample_data(ps.y)[sample_data(ps.y)$Neoplasia.st %in% c("N0.CD", "Nx.CD"),]

sig_by_Wilcox <- colnames(otus)[apply(otus, 2, function(x) wilcox.test(x ~ samp$Neoplasia.st)$p.value) < 0.05]
sigtab <- sigtab[rownames(sigtab) %in% sig_by_Wilcox,]

#nothing to plot


####### Figure 2G #######

### N0.Control versus N0.CD

# 4. Select significance threshold for FDR threshold
alpha.var <- 0.2

# 5. Extract relevant results (N0 versus Nx)
res = results(deseq2, contrast=c("Neoplasia.st", "N0.Control", "N0.CD"))

# 6. Order results at chosen level of (FDR-adjusted) significance
res = res[order(res$padj, na.last=NA), ]
sigtab = res[(res$padj < alpha.var), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.x)[rownames(sigtab), ], "matrix"))
sigtab = sigtab[order(sigtab$log2FoldChange, decreasing = FALSE),]
sigtab$Genus <- as.character(sigtab$Genus)
sigtab$Genus[is.na(sigtab$Genus)] <- as.character(sigtab$Family[is.na(sigtab$Genus)])
sigtab$Genus = factor(sigtab$Genus, levels = unique(sigtab$Genus))

#only keep significant by Wilcoxon test
ps.y <- transform_sample_counts(ps.x, function(x) {x/sum(x)})
ps.y <- subset_taxa(ps.y, rownames(tax_table(ps.x)) %in% rownames(sigtab))
otus <- as.data.frame(as.matrix(otu_table(ps.y)))
otus <- otus[sample_data(ps.y)$Neoplasia.st %in% c("N0.Control", "N0.CD"),]
samp <- sample_data(ps.y)[sample_data(ps.y)$Neoplasia.st %in% c("N0.Control", "N0.CD"),]

sig_by_Wilcox <- colnames(otus)[apply(otus, 2, function(x) wilcox.test(x ~ samp$Neoplasia.st)$p.value) < 0.05]
sigtab <- sigtab[rownames(sigtab) %in% sig_by_Wilcox,]
sigtab$Genus <- as.character(sigtab$Genus)
#manually shorten names
sigtab$Genus <- c("Escherichia", "Eubacteri_no")

ggplot(sigtab, aes(y = Genus, x = log2FoldChange, color = Phylum)) + geom_point(size=3) + theme_classic() + 
  theme(text = element_text(size=22), legend.text = element_text(size=14)) + scale_color_manual(values = dd.col) +
  geom_vline(xintercept = 0, colour = "Red", lty = 2) + theme(legend.position="right") + 
  ggtitle(paste("Control\n", "versus\n", "CD (no dysplasia)", sep = ""))




### N0.Control versus N0.UC - no plot

# 4. Select significance threshold for FDR threshold
alpha.var <- 0.2

# 5. Extract relevant results (N0 versus Nx)
res = results(deseq2, contrast=c("Neoplasia.st", "N0.Control", "N0.UC"))

# 6. Order results at chosen level of (FDR-adjusted) significance
res = res[order(res$padj, na.last=NA), ]
sigtab = res[(res$padj < alpha.var), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.x)[rownames(sigtab), ], "matrix"))
sigtab = sigtab[order(sigtab$log2FoldChange, decreasing = FALSE),]
sigtab$Genus <- as.character(sigtab$Genus)
sigtab$Genus[is.na(sigtab$Genus)] <- as.character(sigtab$Family[is.na(sigtab$Genus)])
sigtab$Genus = factor(sigtab$Genus, levels = unique(sigtab$Genus))

#only keep significant by Wilcoxon test
ps.y <- transform_sample_counts(ps.x, function(x) {x/sum(x)})
ps.y <- subset_taxa(ps.y, rownames(tax_table(ps.x)) %in% rownames(sigtab))
otus <- as.data.frame(as.matrix(otu_table(ps.y)))
otus <- otus[sample_data(ps.y)$Neoplasia.st %in% c("N0.Control", "N0.UC"),]
samp <- sample_data(ps.y)[sample_data(ps.y)$Neoplasia.st %in% c("N0.Control", "N0.UC"),]

sig_by_Wilcox <- colnames(otus)[apply(otus, 2, function(x) wilcox.test(x ~ samp$Neoplasia.st)$p.value) < 0.05]
sigtab <- sigtab[rownames(sigtab) %in% sig_by_Wilcox,]

#no results to plot



### N0.Control versus Nx.CD - no plot

# 4. Select significance threshold for FDR threshold
alpha.var <- 0.2

# 5. Extract relevant results (N0 versus Nx)
res = results(deseq2, contrast=c("Neoplasia.st", "N0.Control", "Nx.CD"))

# 6. Order results at chosen level of (FDR-adjusted) significance
res = res[order(res$padj, na.last=NA), ]
sigtab = res[(res$padj < alpha.var), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.x)[rownames(sigtab), ], "matrix"))
sigtab = sigtab[order(sigtab$log2FoldChange, decreasing = FALSE),]
sigtab$Genus <- as.character(sigtab$Genus)
sigtab$Genus[is.na(sigtab$Genus)] <- as.character(sigtab$Family[is.na(sigtab$Genus)])
sigtab$Genus = factor(sigtab$Genus, levels = unique(sigtab$Genus))

#only keep significant by Wilcoxon test
ps.y <- transform_sample_counts(ps.x, function(x) {x/sum(x)})
ps.y <- subset_taxa(ps.y, rownames(tax_table(ps.x)) %in% rownames(sigtab))
otus <- as.data.frame(as.matrix(otu_table(ps.y)))
otus <- otus[sample_data(ps.y)$Neoplasia.st %in% c("N0.Control", "Nx.CD"),]
samp <- sample_data(ps.y)[sample_data(ps.y)$Neoplasia.st %in% c("N0.Control", "Nx.CD"),]

sig_by_Wilcox <- colnames(otus)[apply(otus, 2, function(x) wilcox.test(x ~ samp$Neoplasia.st)$p.value) < 0.05]
sigtab <- sigtab[rownames(sigtab) %in% sig_by_Wilcox,]

#nothing to plot



### N0.Control versus Nx.UC - no plot

# 4. Select significance threshold for FDR threshold
alpha.var <- 0.2

# 5. Extract relevant results (N0 versus Nx)
res = results(deseq2, contrast=c("Neoplasia.st", "N0.Control", "Nx.UC"))

# 6. Order results at chosen level of (FDR-adjusted) significance
res = res[order(res$padj, na.last=NA), ]
sigtab = res[(res$padj < alpha.var), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.x)[rownames(sigtab), ], "matrix"))
sigtab = sigtab[order(sigtab$log2FoldChange, decreasing = FALSE),]
sigtab$Genus <- as.character(sigtab$Genus)
sigtab$Genus[is.na(sigtab$Genus)] <- as.character(sigtab$Family[is.na(sigtab$Genus)])
sigtab$Genus = factor(sigtab$Genus, levels = unique(sigtab$Genus))

#only keep significant by Wilcoxon test
ps.y <- transform_sample_counts(ps.x, function(x) {x/sum(x)})
ps.y <- subset_taxa(ps.y, rownames(tax_table(ps.x)) %in% rownames(sigtab))
otus <- as.data.frame(as.matrix(otu_table(ps.y)))
otus <- otus[sample_data(ps.y)$Neoplasia.st %in% c("N0.Control", "Nx.UC"),]
samp <- sample_data(ps.y)[sample_data(ps.y)$Neoplasia.st %in% c("N0.Control", "Nx.UC"),]
#changed as only single genus (Butyrivibrio) significant by DESeq2
sig_by_Wilcox <- wilcox.test(otus ~ samp$Neoplasia.st)$p.value < 0.05
sigtab <- sigtab[sig_by_Wilcox,]

#nothing to plot


############ Combined IBD subtypes #############

sample_data(ps2)$Neoplasia.st.ibd <- paste(as.character(sample_data(ps2)$Neoplasia2),as.character(sample_data(ps2)$Subtype), sep = ".")
sample_data(ps2)$Neoplasia.st.ibd[sample_data(ps2)$Neoplasia.st.ibd == "N0.UC"] <- "N0.IBD"
sample_data(ps2)$Neoplasia.st.ibd[sample_data(ps2)$Neoplasia.st.ibd == "N0.CD"] <- "N0.IBD"
sample_data(ps2)$Neoplasia.st.ibd[sample_data(ps2)$Neoplasia.st.ibd == "Nx.UC"] <- "Nx.IBD"
sample_data(ps2)$Neoplasia.st.ibd[sample_data(ps2)$Neoplasia.st.ibd == "Nx.CD"] <- "Nx.IBD"
sample_data(ps2)$Neoplasia.st.ibd[sample_data(ps2)$Neoplasia.st.ibd == "A1.UC"] <- "A1.IBD"
sample_data(ps2)$Neoplasia.st.ibd[sample_data(ps2)$Neoplasia.st.ibd == "A1.CD"] <- "A1.IBD"

sample_data(ps2)$Neoplasia.st.ibd <- factor(sample_data(ps2)$Neoplasia.st.ibd, 
                                            levels = c("N0.Control", "A1.Control", "Nx.Control",
                                                       "N0.IBD", "A1.IBD", "Nx.IBD"))




# 1. Assign to ps.x
ps.x <- ps2
ps.x <- tax_glom(ps.x, "Genus", NArm = FALSE)


# 2. Make a colour vector to keep phylum colour assignments consistent between plots
dd <- unique(tax_table(ps.x)[,2])
dd.col <- c(brewer.pal(8, "Set2"), "#000000", brewer.pal(9, "Set1"))
dd.col <- c(brewer.pal(8, "Set1")[2], brewer.pal(8, "Set1")[4], brewer.pal(8, "Set1")[5], brewer.pal(8, "Set1")[1], "#000000", brewer.pal(9, "Set1"))

names(dd.col)  <- dd
dd.col <- dd.col[1:12]



# 3. Run DESeq2
deseq2 = phyloseq_to_deseq2(ps.x, ~ Neoplasia.st.ibd) #manually add this
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
deseq2_geoMeans = apply(counts(deseq2), 1, gm_mean)
deseq2 = estimateSizeFactors(deseq2, geoMeans = deseq2_geoMeans)
deseq2 = DESeq(deseq2, fitType="local")


####### Figure 2E #######

### N0.IBD versus Nx.IBD

# 4. Select significance threshold for FDR threshold
alpha.var <- 0.2

# 5. Extract relevant results (N0 versus Nx)
res = results(deseq2, contrast=c("Neoplasia.st.ibd", "N0.IBD", "Nx.IBD"))

# 6. Order results at chosen level of (FDR-adjusted) significance
res = res[order(res$padj, na.last=NA), ]
sigtab = res[(res$padj < alpha.var), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.x)[rownames(sigtab), ], "matrix"))
sigtab = sigtab[order(sigtab$log2FoldChange, decreasing = FALSE),]
sigtab$Genus <- as.character(sigtab$Genus)
sigtab$Genus[is.na(sigtab$Genus)] <- as.character(sigtab$Family[is.na(sigtab$Genus)])
sigtab$Genus = factor(sigtab$Genus, levels = unique(sigtab$Genus))

#only keep significant by Wilcoxon test
ps.y <- transform_sample_counts(ps.x, function(x) {x/sum(x)})
ps.y <- subset_taxa(ps.y, rownames(tax_table(ps.x)) %in% rownames(sigtab))
otus <- as.data.frame(as.matrix(otu_table(ps.y)))
otus <- otus[sample_data(ps.y)$Neoplasia.st.ibd %in% c("N0.IBD", "Nx.IBD"),]
samp <- sample_data(ps.y)[sample_data(ps.y)$Neoplasia.st.ibd %in% c("N0.IBD", "Nx.IBD"),]

sig_by_Wilcox <- colnames(otus)[apply(otus, 2, function(x) wilcox.test(x ~ samp$Neoplasia.st.ibd)$p.value) < 0.05]
sigtab <- sigtab[rownames(sigtab) %in% sig_by_Wilcox,]

ggplot(sigtab, aes(y = Genus, x = log2FoldChange, color = Phylum)) + geom_point(size=3) + theme_classic() + 
  theme(text = element_text(size=22), legend.text = element_text(size=14)) + scale_color_manual(values = dd.col) +
  geom_vline(xintercept = 0, colour = "Red", lty = 2) + theme(legend.position="right") + 
  ggtitle(paste("IBD: No dysplasia\n", "versus\n", "Dysplasia (any grade)", sep = ""))




### N0.Control versus N0.IBD - no plot

# 4. Select significance threshold for FDR threshold
alpha.var <- 0.2

# 5. Extract relevant results (N0 versus Nx)
res = results(deseq2, contrast=c("Neoplasia.st.ibd", "N0.Control", "N0.IBD"))

# 6. Order results at chosen level of (FDR-adjusted) significance
res = res[order(res$padj, na.last=NA), ]
sigtab = res[(res$padj < alpha.var), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.x)[rownames(sigtab), ], "matrix"))
sigtab = sigtab[order(sigtab$log2FoldChange, decreasing = FALSE),]
sigtab$Genus <- as.character(sigtab$Genus)
sigtab$Genus[is.na(sigtab$Genus)] <- as.character(sigtab$Family[is.na(sigtab$Genus)])
sigtab$Genus[grepl("Ambig", sigtab$Genus)] <- paste0(as.character(sigtab$Family[grepl("Ambig", sigtab$Genus)]),".OTU")
sigtab$Genus = factor(sigtab$Genus, levels = unique(sigtab$Genus))

#only keep significant by Wilcoxon test
ps.y <- transform_sample_counts(ps.x, function(x) {x/sum(x)})
ps.y <- subset_taxa(ps.y, rownames(tax_table(ps.x)) %in% rownames(sigtab))
otus <- as.data.frame(as.matrix(otu_table(ps.y)))
otus <- otus[sample_data(ps.y)$Neoplasia.st.ibd %in% c("N0.Control", "N0.IBD"),]
samp <- sample_data(ps.y)[sample_data(ps.y)$Neoplasia.st.ibd %in% c("N0.Control", "N0.IBD"),]

sig_by_Wilcox <- colnames(otus)[apply(otus, 2, function(x) wilcox.test(x ~ samp$Neoplasia.st.ibd)$p.value) < 0.05]
sigtab <- sigtab[rownames(sigtab) %in% sig_by_Wilcox,]

#nothing to plot



####### Figure 2H #######

### N0.Control versus Nx.IBD

# 4. Select significance threshold for FDR threshold
alpha.var <- 0.2

# 5. Extract relevant results (N0 versus Nx)
res = results(deseq2, contrast=c("Neoplasia.st.ibd", "N0.Control", "Nx.IBD"))

# 6. Order results at chosen level of (FDR-adjusted) significance
res = res[order(res$padj, na.last=NA), ]
sigtab = res[(res$padj < alpha.var), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.x)[rownames(sigtab), ], "matrix"))
sigtab = sigtab[order(sigtab$log2FoldChange, decreasing = FALSE),]
sigtab$Genus <- as.character(sigtab$Genus)
sigtab$Genus[is.na(sigtab$Genus)] <- as.character(sigtab$Family[is.na(sigtab$Genus)])
sigtab$Genus[grepl("Ambig", sigtab$Genus)] <- paste(as.character(sigtab$Family[grepl("Ambig", sigtab$Genus)]),".OTU")
sigtab$Genus = factor(sigtab$Genus, levels = unique(sigtab$Genus))

#only keep significant by Wilcoxon test
ps.y <- transform_sample_counts(ps.x, function(x) {x/sum(x)})
ps.y <- subset_taxa(ps.y, rownames(tax_table(ps.x)) %in% rownames(sigtab))
otus <- as.data.frame(as.matrix(otu_table(ps.y)))
otus <- otus[sample_data(ps.y)$Neoplasia.st.ibd %in% c("N0.Control", "Nx.IBD"),]
samp <- sample_data(ps.y)[sample_data(ps.y)$Neoplasia.st.ibd %in% c("N0.Control", "Nx.IBD"),]

sig_by_Wilcox <- colnames(otus)[apply(otus, 2, function(x) wilcox.test(x ~ samp$Neoplasia.st.ibd)$p.value) < 0.05]
sigtab <- sigtab[rownames(sigtab) %in% sig_by_Wilcox,]

ggplot(sigtab, aes(y = Genus, x = log2FoldChange, color = Phylum)) + geom_point(size=3) + theme_classic() + 
  theme(text = element_text(size=22), legend.text = element_text(size=14)) + scale_color_manual(values = dd.col) +
  geom_vline(xintercept = 0, colour = "Red", lty = 2) + theme(legend.position="right") + 
  ggtitle(paste("Control (N0)\n", "versus\n", "IBD (dysplasia)", sep = ""))



####### Figure 4C #######
#cluster DESeq2

# 1. Assign to ps.x RNAlater_Frozen
ps.x <- ps2.ibd
ps.x <- tax_glom(ps.x, "Genus", NArm = FALSE)
sample_data(ps.x)$Cluster <- factor(sample_data(ps.x)$Cluster, levels = c("C.1", "C.2", "C.3"))
sample_data(ps.x)$Cluster.comb <- as.character(sample_data(ps.x)$Cluster)
sample_data(ps.x)$Cluster.comb[sample_data(ps.x)$Cluster.comb == "C.1"] <- "C.1_2"
sample_data(ps.x)$Cluster.comb[sample_data(ps.x)$Cluster.comb == "C.2"] <- "C.1_2"

sample_data(ps.x)$Cluster.comb <- factor(sample_data(ps.x)$Cluster.comb, levels = c("C.1_2", "C.3"))

# 2. Make a colour vector to keep phylum colour assignments consistent between plots
dd <- unique(tax_table(ps.x)[,2])
dd.col <- c(brewer.pal(8, "Set2"), "#000000", brewer.pal(9, "Set1"))
dd.col <- c(brewer.pal(8, "Set1")[2], brewer.pal(8, "Set1")[4], brewer.pal(8, "Set1")[5], brewer.pal(8, "Set1")[1], "#000000", brewer.pal(9, "Set1"))

names(dd.col)  <- dd


# 3. Run DESeq2 Cluster combined
deseq2 = phyloseq_to_deseq2(ps.x, ~ Cluster.comb) #manually add this
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
deseq2_geoMeans = apply(counts(deseq2), 1, gm_mean)
deseq2 = estimateSizeFactors(deseq2, geoMeans = deseq2_geoMeans)
deseq2 = DESeq(deseq2, fitType="local")


# 4. Select significance threshold for FDR threshold
alpha.var <- 0.2

# 5. Extract relevant results (C.3 versus C.1_2)
res = results(deseq2, contrast=c("Cluster.comb", "C.3", "C.1_2"))

# 6. Order results at chosen level of (FDR-adjusted) significance
res = res[order(res$padj, na.last=NA), ]
sigtab = res[(res$padj < alpha.var), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.x)[rownames(sigtab), ], "matrix"))
sigtab = sigtab[order(sigtab$log2FoldChange, decreasing = FALSE),]
sigtab$Genus = as.character(sigtab$Genus)
sigtab$Genus[grepl("uncult", sigtab$Genus)] <- paste0("Uncultured_", as.character(sigtab$Family[grepl("uncult", sigtab$Genus)]))
sigtab$Genus <- gsub("_group", "", sigtab$Genus)
sigtab$Genus = factor(sigtab$Genus, levels = unique(sigtab$Genus))


#only keep significant by Wilcoxon test
ps.y <- transform_sample_counts(ps.x, function(x) {x/sum(x)})
ps.y <- subset_taxa(ps.y, rownames(tax_table(ps.x)) %in% rownames(sigtab))
otus <- as.data.frame(as.matrix(otu_table(ps.y)))
otus <- otus[sample_data(ps.y)$Cluster.comb %in% c("C.3", "C.1_2"),]
samp <- sample_data(ps.y)[sample_data(ps.y)$Cluster.comb %in% c("C.3", "C.1_2"),]

sig_by_Wilcox <- colnames(otus)[apply(otus, 2, function(x) wilcox.test(x ~ samp$Cluster.comb)$p.value) < 0.05]
sigtab <- sigtab[rownames(sigtab) %in% sig_by_Wilcox,]

pdf(paste0("DESeq2/", paste("C.3vsC.1_2", sep = "_"), ".pdf"), width = 10, height = 6)
ggplot(sigtab, aes(y = Genus, x = log2FoldChange, colour = Phylum)) + geom_point(size=3) + theme_classic() + 
  theme(text = element_text(size=22), legend.text = element_text(size=14)) + scale_color_manual(values = dd.col) +
  geom_vline(xintercept = 0, colour = "Red", lty = 2) + theme(legend.position="right") + 
  ggtitle(paste("Cluster 3\n", "versus\n", "Clusters 1 & 2", sep = ""))
dev.off()




