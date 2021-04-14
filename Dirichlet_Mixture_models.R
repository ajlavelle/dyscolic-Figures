### Determining optimal number of clusters using Dirichlet Multinomial Mixtures (Holmes I, Harris K, Quince C. PlosONE, 2012: doi.org/10.1371/journal.pone.0030126) ###

library(ggrepel)
library(ggplot2)
library(DirichletMultinomial)

# Agglomoration at genus level

ps2.g <- tax_glom(ps2, "Genus", NArm = FALSE)

# 1. Subset IBD

ps2.g.ibd <- subset_samples(ps2.g, sample_data(ps2.g)$Subtype != "Control")

# 2. Make OTU table

otus <- otu_table(ps2.g.ibd)
otus <- as.matrix(as.data.frame(otus))

# 3. Assign taxa names
otu.vec <- c()

for(i in 1:nrow(tax_table(ps2.g.ibd))){
  if(is.na(tax_table(ps2.g.ibd)[i,6]) == FALSE  & grepl("uncultur", ignore.case = TRUE, tax_table(ps2.g.ibd)[i,6]) == FALSE){
    otu.vec <- c(otu.vec, paste("g", tax_table(ps2.g.ibd)[i,6], sep = "_"))
  } else if(is.na(tax_table(ps2.g.ibd)[i,5]) == FALSE  & grepl("uncultur", ignore.case = TRUE, tax_table(ps2.g.ibd)[i,5]) == FALSE){
    otu.vec <- c(otu.vec, paste("f", tax_table(ps2.g.ibd)[i,5], sep = "_"))
  } else if(is.na(tax_table(ps2.g.ibd)[i,4]) == FALSE  & grepl("uncultur", ignore.case = TRUE, tax_table(ps2.g.ibd)[i,4]) == FALSE){
    otu.vec <- c(otu.vec, paste("o", tax_table(ps2.g.ibd)[i,4], sep = "_"))
  } else{
    otu.vec <- c(otu.vec, paste("c", tax_table(ps2.g.ibd)[i,3], sep = "_"))
  }
}

colnames(otus) <- otu.vec


# 4. Cluster - even with set.seed there appears to be some vriability in this
set.seed(121)
fit <- mclapply(1:7, dmn, count=otus, verbose=TRUE)


# 5. Laplace
lplc <- sapply(fit, laplace)
plot(lplc, type="b", xlab="Number of Dirichlet Components",
     ylab="Model Fit")

# 6. Find best

(best <- fit[[which.min(lplc)]])

mixturewt(best)

#confirm match before assigning
all(names(mixture(best, assign = TRUE)) == rownames(sample_data(ps2.ibd)))


sample_data(ps2.ibd)$Cluster <- paste("C", mixture(best, assign = TRUE), sep = ".")




# 7. Identify bacteria

p0 <- fitted(fit[[1]], scale=TRUE)     # scale by theta
p3 <- fitted(best, scale=TRUE)
colnames(p3) <- paste("m", 1:3, sep="")
(meandiff <- colSums(abs(p3 - as.vector(p0))))
sum(meandiff)

diff <- rowSums(abs(p3 - as.vector(p0)))
o <- order(diff, decreasing=TRUE)
cdiff <- cumsum(diff[o]) / sum(diff)
df <- as.data.frame(head(cbind(Mean=p0[o], p3[o,], diff=diff[o], cdiff), 40))



####### Figure 3E #######

# 8. Pie chart

library(scatterpie)
y.axis <- seq(0,0.125,length.out=40)
y.axis <- seq(0,3,length.out=40)

x.axis <- apply(df[,c("m1", "m2", "m3")], 1, function(x) c(1,2,3)[max(x) == x])
df$x.axis <- x.axis
df$y.axis <- y.axis

y.axis <- rownames(df)
x.axis <- apply(df[,c("m1", "m2", "m3")], 1, function(x) c("C.1","C.2","C.3")[max(x) == x])

df[,c("m1", "m2", "m3")] <- t(apply(df[,c("m1", "m2", "m3")], 1, function(x) x*100/sum(x)))
df[,c("m1", "m2", "m3")] <- t(apply(df[,c("m1", "m2", "m3")], 1, function(x) round(x)))

ggplot() + geom_scatterpie(data = df, aes(x=x.axis, y=rev(y.axis), r=sqrt(Mean*100)/15), cols = c("m1", "m2", "m3"), 
                           color = NA, alpha = 0.9) + 
  scale_fill_manual(values = pal) + theme_classic() + 
  geom_text(data = df, inherit.aes = FALSE, aes(2, rev(y.axis), label = rownames(df)), size = 3.5, hjust="outward") + 
  scale_color_manual(values = pal) + geom_scatterpie_legend(sqrt(df$Mean*100)/15, x=2.5, y=1)



#############################################################

                 ### Relative risk ###

#############################################################


# Cross tabulate
ps2.cd <- subset_samples(ps2.ibd, sample_data(ps2.ibd)$Subtype == "CD")
ps2.uc <- subset_samples(ps2.ibd, sample_data(ps2.ibd)$Subtype == "UC")

counts.cd <- table(sample_data(ps2.cd)$Neoplasia, sample_data(ps2.cd)$Cluster)
counts.uc <- table(sample_data(ps2.uc)$Neoplasia, sample_data(ps2.uc)$Cluster)

cd.prop <- sweep(counts.cd, 2, colSums(counts.cd), "/")
uc.prop <- sweep(counts.uc, 2, colSums(counts.uc), "/")



####### Figure 5A and 5b #######
# Make a bar graph
library(tidyr)
pal <- c("blue", "purple", "pink2", "salmon", "red")

cd.prop.long <- gather(as.data.frame(cbind(cd.prop, Neoplasia = rownames(cd.prop))), key = "Cluster", value = "Percentage", C.1:C.3)
uc.prop.long <- gather(as.data.frame(cbind(uc.prop, Neoplasia = rownames(uc.prop))), key = "Cluster", value = "Percentage", C.1:C.3)
cd.prop.long$Percentage <- as.numeric(as.character(cd.prop.long$Percentage))
uc.prop.long$Percentage <- as.numeric(as.character(uc.prop.long$Percentage))

cd.prop.long$Percentage <- cd.prop.long$Percentage*100
uc.prop.long$Percentage <- uc.prop.long$Percentage*100
x <- rbind(cbind(cd.prop.long, Condition = rep("CD", nrow(cd.prop.long))), cbind(uc.prop.long, Condition = rep("UC", nrow(uc.prop.long))))
x$Neoplasia <- factor(x$Neoplasia, levels = c("N0", "A1", "N1", "N2", "N3"))
ggplot(x, aes(x=Neoplasia, y=Percentage, fill=Neoplasia)) + geom_bar(stat = "identity") + 
  facet_grid(Condition~Cluster) + theme_classic() + scale_fill_manual(values = pal) +
  geom_text(aes(label=round(Percentage, 1)), vjust=-0.3, size=3) + ylim(0,100)


# bar graph all neoplasia versus none

cd.dysplasia <- rbind(colSums(counts.cd[1:2,]), colSums(counts.cd[3:5,]))
rownames(cd.dysplasia) <- c("N0/A1", "N1-N3")
uc.dysplasia <- rbind(colSums(counts.uc[1:2,]), colSums(counts.uc[3:5,]))
rownames(uc.dysplasia) <- c("N0/A1", "N1-N3")


cd.prop <- sweep(cd.dysplasia, 2, colSums(cd.dysplasia), "/")
uc.prop <- sweep(uc.dysplasia, 2, colSums(uc.dysplasia), "/")

cd.prop.long <- gather(as.data.frame(cbind(cd.prop, Neoplasia = rownames(cd.prop))), key = "Cluster", value = "Percentage", C.1:C.3)
uc.prop.long <- gather(as.data.frame(cbind(uc.prop, Neoplasia = rownames(uc.prop))), key = "Cluster", value = "Percentage", C.1:C.3)
cd.prop.long$Percentage <- as.numeric(as.character(cd.prop.long$Percentage))
uc.prop.long$Percentage <- as.numeric(as.character(uc.prop.long$Percentage))

cd.prop.long$Percentage <- cd.prop.long$Percentage*100
uc.prop.long$Percentage <- uc.prop.long$Percentage*100
x <- rbind(cbind(cd.prop.long, Condition = rep("CD", nrow(cd.prop.long))), cbind(uc.prop.long, Condition = rep("UC", nrow(uc.prop.long))))

ggplot(x, aes(x=Neoplasia, y=Percentage, fill=Neoplasia)) + geom_bar(stat = "identity") + 
  facet_grid(Condition~Cluster) + theme_classic() + scale_fill_manual(values = c("blue", "red")) +
  geom_text(aes(label=round(Percentage, 1)), vjust=-0.3, size=3) + ylim(0,100)


#ensure correct variable
sample_data(ps2.ibd)$Neoplasia3 <- as.character(sample_data(ps2.ibd)$Neoplasia)
sample_data(ps2.ibd)$Neoplasia3[sample_data(ps2.ibd)$Neoplasia3 == "N0"] <- "N0/A1"
sample_data(ps2.ibd)$Neoplasia3[sample_data(ps2.ibd)$Neoplasia3 == "A1"] <- "N0/A1"
sample_data(ps2.ibd)$Neoplasia3[sample_data(ps2.ibd)$Neoplasia3 == "N1"] <- "N1-N3"
sample_data(ps2.ibd)$Neoplasia3[sample_data(ps2.ibd)$Neoplasia3 == "N2"] <- "N1-N3"
sample_data(ps2.ibd)$Neoplasia3[sample_data(ps2.ibd)$Neoplasia3 == "N3"] <- "N1-N3"



# Odds ratio
ps2.cd <- subset_samples(ps2.ibd, sample_data(ps2.ibd)$Subtype == "CD")
ps2.uc <- subset_samples(ps2.ibd, sample_data(ps2.ibd)$Subtype == "UC")

counts.cd <- table(sample_data(ps2.cd)$Neoplasia, sample_data(ps2.cd)$Cluster)
counts.uc <- table(sample_data(ps2.uc)$Neoplasia, sample_data(ps2.uc)$Cluster)

cd.prop <- sweep(counts.cd, 2, colSums(counts.cd), "/")
uc.prop <- sweep(counts.uc, 2, colSums(counts.uc), "/")



#relative risk for UC - no neoplasia versus neoplasia

#remove adenomas
ps2.uc.noA1 <- subset_samples(ps2.uc, sample_data(ps2.uc)$Neoplasia != "A1")
uc_rr <- table(sample_data(ps2.uc.noA1)$Neoplasia3, sample_data(ps2.uc.noA1)$Cluster)

#set i=1 for C.1 and i=2 for C.2
i=2
j=3


#RR
rel.risk <- (uc_rr[2,i]/(uc_rr[1,i]+uc_rr[2,i]))/(uc_rr[2,j]/(uc_rr[1,j]+uc_rr[2,j]))
rel.risk
#95% CI's
exp(log(rel.risk) - 1.96*sqrt((1/uc_rr[2,i]) + (1/uc_rr[2,j]) - (1/(uc_rr[1,i]+uc_rr[2,i])) - (1/(uc_rr[1,j]+uc_rr[2,j]))))
exp(log(rel.risk) + 1.96*sqrt((1/uc_rr[2,i]) + (1/uc_rr[2,j]) - (1/(uc_rr[1,i]+uc_rr[2,i])) - (1/(uc_rr[1,j]+uc_rr[2,j]))))

#p-value (Sheskin DJ, 2004)
q.val <- log(rel.risk)/sqrt((1/uc_rr[2,i]) + (1/uc_rr[2,j]) - (1/(uc_rr[1,i]+uc_rr[2,i])) - (1/(uc_rr[1,j]+uc_rr[2,j])))

p.val <- (1-pnorm(q.val))*2
p.val

