###### Logistic regression ######
mod.mat <- as.data.frame(as.matrix(sample_data(ps2.bile.ibd))) # using the phyloseq object which includes IBD with bile acid profiles

## merge variables for better levels

#Montreal
mod.mat$Montreal_combined <- as.character(mod.mat$Montreal_combined)
mod.mat$Montreal_combined[mod.mat$Montreal_combined == "L1" & !(is.na(mod.mat$Montreal_combined))] <- "CCD"
mod.mat$Montreal_combined[mod.mat$Montreal_combined == "L2" & !(is.na(mod.mat$Montreal_combined))] <- "CCD"
mod.mat$Montreal_combined[mod.mat$Montreal_combined == "L3" & !(is.na(mod.mat$Montreal_combined))] <- "CCD"
mod.mat$Montreal_combined[mod.mat$Montreal_combined == "E1" & !(is.na(mod.mat$Montreal_combined))] <- "E1_2"
mod.mat$Montreal_combined[mod.mat$Montreal_combined == "E2" & !(is.na(mod.mat$Montreal_combined))] <- "E1_2"
mod.mat <- mod.mat[!(is.na(mod.mat$Montreal_combined)),]

#Age
mod.mat$Age <- as.numeric(as.character(mod.mat$Age))

#Duration
mod.mat$Duration <- as.numeric(as.character(mod.mat$Duration))

#BMI
mod.mat$BMI <- as.numeric(as.character(mod.mat$BMI))

#Age at diagnosis
mod.mat$Age_at_diagnosis <- as.numeric(as.character(mod.mat$Age_at_diagnosis))
mod.mat$Age_at_diagnosis[is.na(mod.mat$Age_at_diagnosis)] <- mod.mat$Age[is.na(mod.mat$Age_at_diagnosis)] - mod.mat$Duration[is.na(mod.mat$Age_at_diagnosis)]

#merge anti-TNFs
mod.mat$antiTNF <- as.character(mod.mat$Infliximab == "Yes" | mod.mat$Adalimumab == "Yes" | mod.mat$Golimumab == "Yes")
mod.mat$Smoke.ever <- as.character(mod.mat$Smoker == "Yes" | mod.mat$Ex.smoker == "Yes")
mod.mat$Rectal_5ASA <- as.character(mod.mat$Rectal_ASA == "Yes" | mod.mat$Enema_ASA == "Yes")

mod.mat$Neoplasia3 <- as.character(mod.mat$Neoplasia3)
mod.mat$Neoplasia3 <- gsub("/","",mod.mat$Neoplasia3)
mod.mat$Neoplasia3 <- gsub("-","",mod.mat$Neoplasia3)
mod.mat$Neoplasia3 <- as.factor(mod.mat$Neoplasia3)

#Select samples only present in both
mod.mat.bile <- mod.mat[rownames(mod.mat) %in% rownames(bile.ibd),]
bile.ibd.mod <- bile.ibd[rownames(bile.ibd) %in% rownames(mod.mat.bile),]
#Normalise the individual bile acids to total
bile.ibd.mod[,1:27] <- sweep(bile.ibd.mod[,1:27],1,bile.ibd.mod$Total_BA,"/")
#Normalise the grouped bile acids to total
bile.ibd.mod[,c(31:37,39:43,51:52)] <- sweep(bile.ibd.mod[,c(31:37,39:43,51:52)],1,bile.ibd.mod$Total_BA,"/")

bile.bin <- bile.ibd.mod[,c(31:53)]

for(i in 1:ncol(bile.bin)){
  temp.numeric.vec <- bile.bin[,i]
  bile.bin[,i] <- as.factor(unlist(lapply(temp.numeric.vec, function(x) if(x > median(temp.numeric.vec)){"High"} else {"Low"})))
  colnames(bile.bin)[i] <- paste(colnames(bile.bin)[i], "bin", sep = ".")
}

#Combine final selected binary bile acids split by median with other variables ensuring sample names match
if(all(rownames(mod.mat.bile) == rownames(bile.bin))){
  mod.mat.bile <- cbind(mod.mat.bile, bile.bin[,c(1:17,19,23)])
} else {
  print("ERROR! Rownames do not match!")
}

#remove missing data and adenomas (A1)
mod.mat.bile <- mod.mat.bile[!(is.na(mod.mat.bile$Smoke.ever)),]
mod.mat.bile <- mod.mat.bile[!(is.na(mod.mat.bile$Appendicectomy)),]

mod.mat.bile <- mod.mat.bile[,((colSums(mod.mat.bile == "No"))/nrow(mod.mat.bile)) <= 0.9 | is.na(colSums(mod.mat.bile == "No"))]
mod.mat.bile <- mod.mat.bile[mod.mat.bile$Neoplasia != "A1",] # remove A1

############ GLM NET ############

library(tidyverse)
library(caret)
library(glmnet)

### variables to select

variables <- c("Age", "Age_at_diagnosis", "Duration", "Gender", "BMI", "Montreal_combined", "PSC_yes_no", 
               "Appendicectomy", "Oral_ASA", "Thiopurine", "UDCA_eng", "Combined_activity.bin", "Neoplasia3",
               "Cluster", "antiTNF", "Smoke.ever", "Rectal_5ASA", "CA.1.bin", "CDCA.1.bin", "DCA.1.bin",                                 
               "LCA.1.bin", "UDCA.1.bin",  "Glyco.BA_Total.bin", "Tauro.BA_Total.bin", "Sulfo.BA_Total.bin", 
               "BA_nonconj.AB_Tot.bin")

sel.data <- mod.mat.bile[,variables]
x <- model.matrix(Neoplasia3~., sel.data)

y <- ifelse(sel.data$Neoplasia3 == "N1N3", 1, 0)

set.seed(1)
cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial")
plot(cv.lasso)


lasso.model <- glmnet(x, y, alpha = 1, family = "binomial",
                      lambda = cv.lasso$lambda.1se)

temp.selected <- lasso.model$beta@Dimnames[[1]][-1][lasso.model$beta@i]

present.vars <- c()
for(i in 1:length(colnames(sel.data))){
  
  present.temp <- temp.selected[grepl(colnames(train.data)[i], temp.selected)]
  if(length(present.temp) > 0){
    pres <- rep(colnames(train.data)[i], length(present.temp))
    present.vars <- c(present.vars, pres)
  }
}
present.vars <- unique(present.vars)

#relevel factors
mod.mat.bile$Cluster <- relevel(mod.mat.bile$Cluster, "C.3")
mod.mat.bile$Tauro.BA_Total.bin <- relevel(mod.mat.bile$Tauro.BA_Total.bin, "Low")
mod.mat.bile$Sulfo.BA_Total.bin <- relevel(mod.mat.bile$Sulfo.BA_Total.bin, "Low")

mod.1 <- glm(Neoplasia3 ~ ., family = "binomial", data = mod.mat.bile[,c("Neoplasia3", present.vars)])

#forward-backward stepwise logistic regression
xx <- step(mod.1, direction = "both")
xx


# Coefficient plots

library(broom)
library(ggstance)
model1 <- xx

m1.coef <- tidy(model1, conf.int = TRUE)

m1.coef <- dplyr::filter(m1.coef, term!="(Intercept)")

#Make variable names more readable
m1.coef$term <- c("Disease extent: E1/E2[UC]", "Disease extent: E3[UC]", "Thiopurine: Yes", 
                  "Clinical activity: Active", "DMM cluster: C.1", "DMM cluster: C.2", "Rectal 5-ASA: Yes", 
                  "Taurine-conjugated BAs: High", "Sulfo-conjugated BAs: High")


####### Figure 6B #######
ggplot(m1.coef, aes(x=estimate, y=factor(term, levels = ))) + geom_point(colour="dodgerblue", size=3) + 
  geom_errorbarh(data=m1.coef,aes(xmin=conf.low, xmax=conf.high), colour="dodgerblue", height=0.2) + theme_classic()  + 
  geom_vline(xintercept = 0, colour="red", lty=2) + ggtitle("Coefficient plot: logistic regression") + xlab("Beta coefficient") +
  scale_y_discrete(limits=rev(levels(m1.coef$term))) + ylab("")



