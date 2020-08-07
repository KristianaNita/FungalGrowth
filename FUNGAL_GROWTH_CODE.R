#Load libraries
library(data.table)
library(stringr)
library(ggplot2)
library(gridExtra)
library(factoextra) #for fviz visualisation
library(cluster) #for clustering
library(seqinr) #fasta reader
library(glmnet) #lasso
#library(tidyverse)

#Load data
cryptoWU.samples <- fread("CryptoWakeUpSampleSheetPlusTags.txt") #36
gene.expression.dt <- fread("CW-kallisto-abundance-foldchange-long-bygene.txt") #293472
promoter.dt <- fread("H99_allorfs_promoter500nt_5mercounts.txt") #6879

#Order samples
cryptoWU.samples <- cryptoWU.samples[order(cryptoWU.samples$Code)]

#Number of genes
length(unique(gene.expression.dt$Gene)) #8152

#Merge gene expression and samplesheets
gene.expression.dt <- merge(cryptoWU.samples, gene.expression.dt, by.x="Code", by.y="Code",all.y=TRUE) #293472
gene.expression.dt <- gene.expression.dt[,!c(7,8)] #erase unnecessary columns 

#Extract specific letters from a string
gene.expression.dt$identifier <- str_sub(gene.expression.dt$Gene, 6, -5)
#Have only coding genes
gene.expression.dt <- gene.expression.dt[identifier == "0"]
gene.expression.dt <- gene.expression.dt[,!11] #erase unnecessary columns 
#250380 rows

#Number of genes
length(unique(gene.expression.dt$Gene)) #6955

#Round Estimated counts into 0 decimals
gene.expression.dt$EstCounts <- as.integer(round(gene.expression.dt$EstCounts,0))

#Find the genes that have 0 counts and remove them 
names.dt <- unique(gene.expression.dt$Gene[gene.expression.dt$EstCounts == 0])
nonzero.gene.expression  <- gene.expression.dt[! gene.expression.dt$Gene %in% names.dt,]

length(unique(nonzero.gene.expression$Gene))

#Normalise variables
nonzero.gene.expression$FoldChange <- log2(nonzero.gene.expression$FoldChange)
nonzero.gene.expression$TPM <- log2(1+nonzero.gene.expression$TPM)

#Create a data table with the genes and foldchanges
gene.foldchanges <- data.frame(matrix(nrow=36, ncol=0))
for(i in 1:length(unique(nonzero.gene.expression$Gene))){
  gene.foldchanges[,i]  = nonzero.gene.expression[which(nonzero.gene.expression$Gene == unique(nonzero.gene.expression$Gene)[i]),10]
}

colnames(gene.foldchanges) <- unique(nonzero.gene.expression$Gene)
rownames(gene.foldchanges) <- unique(nonzero.gene.expression$Code)
design.matrix <- t(gene.foldchanges)

#Create a dataframe
design.matrix <- as.data.frame(design.matrix)

# Plot time 0 and 180 
#0 cold vs 10 cold
w <- as.data.frame(cbind(design.matrix$YC0A, design.matrix$YC1A))
setnames(w, "V1", "YCOA")
setnames(w, "V2", "YC1A")
w.melt <- melt(w) 

p0 <- ggplot(w.melt ) + geom_density(aes(x = value, colour = variable)) + theme_minimal() +
  scale_color_manual(values=c( "#E69F00", "deeppink4")) +
  ggtitle("Densities of YC0A and YC1A") + xlab("Foldchange") + ylab("Density") +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) 

#180 medium vs 120 hot and cold
w1 <- as.data.frame(cbind(design.matrix$YM9A, design.matrix$YC4A, design.matrix$YH4A))
setnames(w1, "V1", "YM9A")
setnames(w1, "V2", "YC4A")
setnames(w1, "V3", "YH4A")
w1.melt <- melt(w1)

p180 <- ggplot(w1.melt ) + geom_density(aes(x = value, colour = variable)) + theme_minimal() +
  scale_color_manual(values=c( "#E69F00", "deeppink4","grey38")) +
  ggtitle("Densities of YM9A,YC4A and YH4A") + xlab("Foldchange") + ylab("Density") +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) 

grid.arrange(p0,p180,nrow=1)

#Remove tha samples with Time 0 and 180
idx_0_180 <- which(nonzero.gene.expression$Time %in% c("0","180"))#19124
nonzero.gene.expression <- nonzero.gene.expression[!idx_0_180,]#152992
design.matrix <- design.matrix[,-c(17,18,35,36)]


#CNAG_03012 largest TPM value
CNAG_03012.RPMI <- nonzero.gene.expression[Gene == "CNAG_03012" & Rep == "A" & Medium == "RPMI+"]
CNAG_03012.YPD <- nonzero.gene.expression[Gene == "CNAG_03012" & Rep == "A" & Medium == "YPD"]

p1 <- ggplot(data=CNAG_03012.RPMI, aes(x=Time, y=TPM, group = Temp, colour = Temp)) +
  theme_minimal() + geom_line() +
  scale_color_manual(values=c('deeppink4','#E69F00'))+
  geom_point( size=3, shape=21, fill="white") +
  ggtitle("RPMI+ and Replicate A") + xlab("Time points") + ylab("TPM") 

p2 <- ggplot(data=CNAG_03012.YPD, aes(x=Time, y=TPM, group = Temp, colour = Temp)) +
  theme_minimal() + geom_line() +
  scale_color_manual(values=c('deeppink4','#E69F00'))+
  geom_point( size=3, shape=21, fill="white") +
  ggtitle("YPD and Replicate A") + xlab("Time points") + ylab("TPM") 



#Outlier CNAG_06917
CNAG_06917.RPMI <- nonzero.gene.expression[Gene == "CNAG_06917" & Rep == "A" & Medium == "RPMI+"]
CNAG_06917.YPD <- nonzero.gene.expression[Gene == "CNAG_06917" & Rep == "A" & Medium == "YPD"]

p3 <- ggplot(data=CNAG_06917.RPMI, aes(x=Time, y=TPM, group = Temp, colour = Temp)) +
  theme_minimal() + geom_line() +
  scale_color_manual(values=c('deeppink4','#E69F00'))+
  geom_point( size=3, shape=21, fill="white") +
  ggtitle("RPMI+ and Replicate A") + xlab("Time points") + ylab("TPM") 

p4 <- ggplot(data=CNAG_06917.YPD, aes(x=Time, y=TPM, group = Temp, colour = Temp)) +
  theme_minimal() + geom_line() +
  scale_color_manual(values=c('deeppink4','#E69F00'))+
  geom_point( size=3, shape=21, fill="white") +
  ggtitle("YPD and Replicate A") + xlab("Time points") + ylab("TPM") 

grid.arrange(arrangeGrob(p1, p2, top="CNAG_03012"), arrangeGrob(p3, p4, top="CNAG_06917"), ncol=2)



#Density of TPM
with(nonzero.gene.expression, hist(FoldChange))
with(nonzero.gene.expression,hist(TPM))


#Histograms of each sample
for (i in 1:32){
  hist(design.matrix[,i])
}


#PCA for factors
factor.pca <- prcomp(gene.foldchanges, scale = TRUE)
fviz_pca_ind(factor.pca, geom=c("point","text"))+
  scale_color_gradient2(low="blue", mid="red",
                        high="green", midpoint=4,)

#Correlation between conditions
design.matrix.dt <- as.data.table(design.matrix)
corr <- cor(design.matrix.dt[,1:32][, .SD, .SDcols = sapply(design.matrix.dt[,1:32], is.numeric)], use = "pairwise.complete")
library(corrplot)
#Plot the correlation matrix
corrplot(corr, order="hclust", method = "color", hclust.method = "ward.D",tl.cex = 0.3,
         diag=FALSE, tl.col="black",title="Correlation matrix (ordered by hierarchical clustering)",
         mar=c(0,1,2,0))



NoReplicates <- as.data.frame(matrix(0,nrow(design.matrix.dt), ncol=16))
NoReplicates[,1] = apply(design.matrix.dt[ ,c(2*1-1, 2*1)],1, mean)
NoReplicates[,2] = apply(design.matrix.dt[ ,c(2*2-1, 2*2)],1, mean)
NoReplicates[,3] = apply(design.matrix.dt[ ,c(2*3-1, 2*3)],1, mean)
NoReplicates[,4] = apply(design.matrix.dt[ ,c(2*4-1, 2*4)],1, mean)
NoReplicates[,5] = apply(design.matrix.dt[ ,c(2*5-1, 2*5)],1, mean)
NoReplicates[,6] = apply(design.matrix.dt[ ,c(2*6-1, 2*6)],1, mean)
NoReplicates[,7] = apply(design.matrix.dt[ ,c(2*7-1, 2*7)],1, mean)
NoReplicates[,8] = apply(design.matrix.dt[ ,c(2*8-1, 2*8)],1, mean)
NoReplicates[,9] = apply(design.matrix.dt[ ,c(2*9-1, 2*9)],1, mean)
NoReplicates[,10] = apply(design.matrix.dt[ ,c(2*10-1, 2*10)],1, mean)
NoReplicates[,11] = apply(design.matrix.dt[ ,c(2*11-1, 2*11)],1, mean)
NoReplicates[,12] = apply(design.matrix.dt[ ,c(2*12-1, 2*12)],1, mean)
NoReplicates[,13] = apply(design.matrix.dt[ ,c(2*13-1, 2*13)],1, mean)
NoReplicates[,14] = apply(design.matrix.dt[ ,c(2*14-1, 2*14)],1, mean)
NoReplicates[,15] = apply(design.matrix.dt[ ,c(2*15-1, 2*15)],1, mean)
NoReplicates[,16] = apply(design.matrix.dt[ ,c(2*16-1, 2*16)],1, mean)


colnames(NoReplicates) <- c("RC1","RC2","RC3","RC4","RH1","RH2","RH3","RH4",
                            "YC1","YC2","YC3","YC4","YH1","YH2","YH3","YH4")
row.names(NoReplicates)<- rownames(design.matrix)

#Standardize the variables
NoReplicates.sc <- scale(NoReplicates)
#Kmeans clustering
# Elbow method
fviz_nbclust(NoReplicates.sc, kmeans, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)+
  labs(subtitle = "Elbow method")


#CLUSTER ANALYSIS
set.seed(2)
GE.km <- kmeans(NoReplicates.sc, centers = 3, nstart=100)
t(GE.km$centers)# Cluster centres
print(GE.km)


# Within-cluster variability and size 
GE.km$withinss # Cluster homogeneity

GE.km$size # Cluster sizes
knitr::kable(t(GE.km$centers),"latex")
apply(t(GE.km$centers), 2, sd)
apply(t(GE.km$centers), 2, mean)

#Plot clusters
fviz_cluster(GE.km, NoReplicates.sc, labelsize = 6, palette = "jco", ggtheme = theme_minimal())


#Cluster allocations
clusters <- GE.km$cluster
outputs <- as.data.table(clusters, keep.rownames = TRUE)
setnames(outputs, "rn", "Genes")


# Append groupinf information
design.matrix.clusters <- as.data.frame(cbind(clusters,NoReplicates.sc)) 

################################### QUESTION 2 #######################################
#PCA
summary(NoReplicates)
apply(NoReplicates,2,sd) # Standard deviations


dm.pca <- prcomp(NoReplicates.sc) # PCA on the correlation matrix
summary(dm.pca)
#Compute the variance explained by the first two components
perc.expl <- dm.pca$sdev^2 / sum(dm.pca$sdev^2)
cat("The variance explained by the first two components is:",
    signif(sum(perc.expl[1:2]),4)*100,"%","\n")
cat("The variance explained by the first three components is:",
    signif(sum(perc.expl[1:3]),4)*100,"%","\n")

dm.pca$rotation #Loadings
predict(dm.pca) # PC scores for the first genes
fviz_eig(dm.pca) # Scree plot

# First 3 PCs scores
comp <- data.frame(dm.pca$x[,1:3])

#3 PCs plot
library(RColorBrewer)
library(scales)
palette(alpha(brewer.pal(5,'Paired'), 0.6))
plot(comp, col=GE.km$clust, pch=16)


#Variables
fviz_pca_var(dm.pca, geom = c("point", "text")) 

#Bibplot
fviz_pca_biplot(dm.pca,label ="var", col.ind=design.matrix.clusters$clusters) +
  scale_color_gradient2(low="blue", mid="green",
                        high="green", midpoint=4)

#Associations of variables with PCs
var <- get_pca_var(dm.pca)
corrplot(var$cos2, is.corr = FALSE)


###################################QUESTION 3#############################
################# Perform lasso regression for PCs
motifs <- promoter.dt[,!1]
motifs <- as.matrix(motifs)

FiveMers <- matrix(0,nrow(motifs),ncol(motifs))
colnames(FiveMers) <- colnames(motifs)

for(i in 1:ncol(motifs)){
  new_name = c()
  nucl =rev(unlist(strsplit(colnames(motifs)[i], "")))
  
  if(sum(which(nucl == "A")>0)>0){
    new_name[which(nucl == "A")] = "T"
  }
  if(sum(which(nucl == "T")>0)>0){
    new_name[which(nucl == "T")] = "A"
  }
  if(sum(which(nucl == "C")>0)>0){
    new_name[which(nucl == "C")] = "G"
  }
  if(sum(which(nucl == "G")>0)>0){
    new_name[which(nucl == "G")] = "C"
  }
  
  FiveMers[,i] <- apply(cbind(motifs[,which(colnames(motifs) == paste(new_name, collapse = ""))],
                              motifs[,which(colnames(motifs) == paste(nucl, collapse = ""))]), 1, mean)
  colnames(FiveMers)[i] <- paste(nucl, collapse = "")
}

FiveMers <- as.data.frame(FiveMers)
FiveMers <- FiveMers[,1:512]
FiveMers$Gene <- promoter.dt[,1]


#Create data frame on PCs
PCs <- as.data.frame(predict(dm.pca)[,1:3])
setDT(PCs,keep.rownames = TRUE)
setnames(PCs, "rn", "Gene")

#Merge promoter kmers dataset and the PCs
fivemers.pc.dt <- merge(PCs, FiveMers, by.x="Gene", by.y="Gene",all.x=TRUE)
fivemers.pc.dt <- fivemers.pc.dt[!which(is.na(fivemers.pc.dt$AAAAA))]

#Create data partition
set.seed(100) 
index = sample(1:nrow(fivemers.pc.dt), 0.7*nrow(fivemers.pc.dt)) 
train = fivemers.pc.dt[index,] # Create the training data 
test = fivemers.pc.dt[-index,] # Create the test data
dim(train)
dim(test)

#Training set
x_train <- as.matrix(train[,5:516])
y_train <- as.matrix(train[,2:4])
#Test set
x_test <- as.matrix(test[,5:516])
y_test <- as.matrix(test[,2:4])



lambda_seq <- 10^seq(3,-3, by=-.1)
#Grouped lasso
lasso1 <- cv.glmnet(x = x_train, y = y_train,
                   lambda = lambda_seq, family = "mgaussian")
#Fit the model
fit_lasso1 <- glmnet(x = x_train, y = y_train,
                        lambda = lasso1$lambda.min, family = "mgaussian")
plot(lasso1)
lasso1$lambda.min
lasso1$lambda.1se

#Coefficients of the model
coefficients<- coef(fit_lasso1)
i.length <- length(coefficients[["PC1"]]@i) #194

#Create a matrix with the coefficients of the PCs
motif.coef <- data.frame(matrix(nrow=i.length, ncol=3))
motif.coef[,1] <- coefficients[["PC1"]]@x
motif.coef[,2] <- coefficients[["PC2"]]@x
motif.coef[,3] <- coefficients[["PC3"]]@x

colnames(motif.coef) <- c("PC1","PC2","PC3")
row.names(motif.coef)<- coefficients[["PC1"]]@Dimnames[[1]][which(coefficients[["PC1"]]@x != 0)]
motif.coef <- motif.coef[-1,]

#Calculate the L2-norm
motif.coef <- cbind(motif.coef, L2_norm = apply(X = motif.coef, MARGIN = 1, FUN = norm, '2'))
motif.coef <- motif.coef[order(-motif.coef$L2_norm),]
head(motif.coef,15)
knitr::kable(head(motif.coef,15),"latex")
#Calculate the sum of squares
mean.coef <- apply(motif.coef[1:3], 2, mean)
sumofsq <- function(x)(sum((x[1]-mean.coef[1])^2,(x[2]-mean.coef[2])^2,(x[3]-mean.coef[3])^2))
motif.coef <- cbind(motif.coef, SumOfSquares = apply(motif.coef, 1, sumofsq) )
motif.coef <- motif.coef[order(-motif.coef$SumOfSquares),]
head(motif.coef[,-4],15)
knitr::kable(head(motif.coef[,-4],15),"latex")



###### QUESTION 3B #####
#Creating 6-mers
six_mers <- fread("promoter_6mer_counts.txt")


motifs6 <- six_mers[,!1]
motifs6 <- as.matrix(motifs6)

SixMers <- matrix(0,nrow(motifs6),ncol(motifs6))
colnames(SixMers) <- colnames(motifs6)

for(i in 1:ncol(motifs6)){
  new_name = c()
  nucl =rev(unlist(strsplit(colnames(motifs6)[i], "")))
  
  if(sum(which(nucl == "A")>0)>0){
    new_name[which(nucl == "A")] = "T"
  }
  if(sum(which(nucl == "T")>0)>0){
    new_name[which(nucl == "T")] = "A"
  }
  if(sum(which(nucl == "C")>0)>0){
    new_name[which(nucl == "C")] = "G"
  }
  if(sum(which(nucl == "G")>0)>0){
    new_name[which(nucl == "G")] = "C"
  }
  
  SixMers[,i] <- apply(cbind(motifs6[,which(colnames(motifs6) == paste(new_name, collapse = ""))],
                              motifs6[,which(colnames(motifs6) == paste(nucl, collapse = ""))]), 1, mean)
  colnames(SixMers)[i] <- paste(nucl, collapse = "")
}

SixMers <- as.data.frame(SixMers)
SixMers <- SixMers[,1:(ncol(motifs6)/2)]
SixMers$Gene <- six_mers[,1]

SixMers.pc.dt <- merge(PCs, SixMers, by.x="Gene", by.y="Gene",all.x=TRUE)
SixMers.pc.dt <- SixMers.pc.dt[!which(is.na(SixMers.pc.dt$AAAAAA))]

set.seed(100) 
index6 = sample(1:nrow(SixMers.pc.dt), 0.7*nrow(SixMers.pc.dt)) 
train6 = SixMers.pc.dt[index6,] # Create the training data 
test6 = SixMers.pc.dt[-index6,] # Create the test data
dim(train6)
dim(test6)


lambda_seq <- 10^seq(3,-3, by=-.1)
#Training set
x1_train<- as.matrix(train6[,5:2052])
y1_train <- as.matrix(train6[,2:4])
#Test set
x1_test<- as.matrix(test6[,5:2052])
y1_test <- as.matrix(test6[,2:4])


lasso11 <- cv.glmnet(x = x1_train, y = y1_train,
                    lambda = lambda_seq, family = "mgaussian")

fit_lasso11 <- glmnet(x = x1_train, y = y1_train,
                     lambda = lasso11$lambda.min, family = "mgaussian")

coef6 <- coef(fit_lasso11)
i.length6 <- length(coef6[["PC1"]]@i)

#Create a matrix with the coefficients of the PCs
motif.coef6 <- data.frame(matrix(nrow=i.length6, ncol=3))
motif.coef6[,1] <- coef6[["PC1"]]@x
motif.coef6[,2] <- coef6[["PC2"]]@x
motif.coef6[,3] <- coef6[["PC3"]]@x

colnames(motif.coef6) <- c("PC1","PC2","PC3")
row.names(motif.coef6)<- coef6[["PC1"]]@Dimnames[[1]][which(coef6[["PC1"]]@x != 0)]
motif.coef6 <- motif.coef6[-1,]

motif.coef6 <- cbind(motif.coef6, L2_norm = apply(X = motif.coef6, MARGIN = 1, FUN = norm, '2'))
motif.coef6 <- motif.coef6[order(-motif.coef6$L2_norm),]
head(motif.coef6,20)
knitr::kable(head(motif.coef6,20),"latex")


#4-mers
four_mers <- fread("promoter_4mer_counts.txt")


motifs4 <- four_mers[,!1]
motifs4 <- as.matrix(motifs4)

FourMers <- matrix(0,nrow(motifs4),ncol(motifs4))
colnames(FourMers) <- colnames(motifs4)

for(i in 1:ncol(motifs4)){
  new_name = c()
  nucl =rev(unlist(strsplit(colnames(motifs4)[i], "")))
  
  if(sum(which(nucl == "A")>0)>0){
    new_name[which(nucl == "A")] = "T"
  }
  if(sum(which(nucl == "T")>0)>0){
    new_name[which(nucl == "T")] = "A"
  }
  if(sum(which(nucl == "C")>0)>0){
    new_name[which(nucl == "C")] = "G"
  }
  if(sum(which(nucl == "G")>0)>0){
    new_name[which(nucl == "G")] = "C"
  }
  
  FourMers[,i] <- apply(cbind(motifs4[,which(colnames(motifs4) == paste(new_name, collapse = ""))],
                             motifs4[,which(colnames(motifs4) == paste(nucl, collapse = ""))]), 1, mean)
  colnames(FourMers)[i] <- paste(nucl, collapse = "")
}

FourMers <- as.data.frame(FourMers)
FourMers <- FourMers[,1:(ncol(motifs4)/2)]
FourMers$Gene <- four_mers[,1]

FourMers.pc.dt <- merge(PCs, FourMers, by.x="Gene", by.y="Gene",all.x=TRUE)
FourMers.pc.dt <- FourMers.pc.dt[!which(is.na(FourMers.pc.dt$AAAA))]

set.seed(100) 
index4 = sample(1:nrow(FourMers.pc.dt), 0.7*nrow(FourMers.pc.dt)) 
train4 = FourMers.pc.dt[index4,] # Create the training data 
test4 = FourMers.pc.dt[-index4,] # Create the test data
dim(train4)
dim(test4)


lambda_seq <- 10^seq(3,-3, by=-.1)
#Training set
x2_train<- as.matrix(train4[,5:132])
y2_train <- as.matrix(train4[,2:4])

#Test set
x2_test<- as.matrix(test4[,5:132])
y2_test <- as.matrix(test4[,2:4])


lasso2 <- cv.glmnet(x = x2_train, y = y2_train,
                     lambda = lambda_seq, family = "mgaussian")

fit_lasso2 <- glmnet(x = x2_train, y = y2_train,
                      lambda = lasso2$lambda.min, family = "mgaussian")

coef4 <- coef(fit_lasso2)
i.length4 <- length(coef4[["PC1"]]@i)

#Create a matrix with the coefficients of the PCs
motif.coef4 <- data.frame(matrix(nrow=i.length4, ncol=3))
motif.coef4[,1] <- coef4[["PC1"]]@x
motif.coef4[,2] <- coef4[["PC2"]]@x
motif.coef4[,3] <- coef4[["PC3"]]@x

colnames(motif.coef4) <- c("PC1","PC2","PC3")
row.names(motif.coef4)<- coef4[["PC1"]]@Dimnames[[1]][which(coef4[["PC1"]]@x != 0)]
motif.coef4 <- motif.coef4[-1,]

motif.coef4 <- cbind(motif.coef4, L2_norm = apply(X = motif.coef4, MARGIN = 1, FUN = norm, '2'))
motif.coef4 <- motif.coef4[order(-motif.coef4$L2_norm),]
head(motif.coef4,10)
knitr::kable(head(motif.coef4,10),"latex")













