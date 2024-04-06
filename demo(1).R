###################################################
#This is for final project
#FPG¿Õ¸¹ÑªÌÇ£¬PPG·¹ºóÑªÌÇ
rm(list=ls())
library(ggplot2)
library(GGally)
library(factoextra)
library(gridExtra)
library(corrplot)
protein <- read.csv('protein/protein.csv')
head(protein)
str(protein)
#######~~Descriptive analysis~~~~~~###################~~~~~~~
summary(protein)
nrow(protein)

#Removing the missing values
protein <- na.omit(protein)

set.seed(20191126)

#calulate the correlation coefficient
cor <- cor(protein[,-1])
round(cor, 2)

#correlation plot
corrplot(cor, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

#density plot
p1<-ggplot(data=protein,aes(RedMeat,..density..))+
  geom_histogram(color='white',fill='gray60',binwidth = 3)+ 
  geom_line(stat='density',colour="red",size=0.7)
p2<-ggplot(data=protein,aes(WhiteMeat,..density..))+
  geom_histogram(color='white',fill='gray60',binwidth = 3)+
  geom_line(stat='density',colour="red",size=0.7)
p3<-ggplot(data=protein,aes(Eggs,..density..))+
  geom_histogram(color='white',fill='gray60',binwidth = 3)+
  geom_line(stat='density',colour="red",size=0.7)
p4<-ggplot(data=protein,aes(Milk,..density..))+
  geom_histogram(color='white',fill='gray60',binwidth = 3)+
  geom_line(stat='density',colour="red",size=0.7)
p5<-ggplot(data=protein,aes(Fish,..density..))+
  geom_histogram(color='white',fill='gray60',binwidth = 3)+
  geom_line(stat='density',colour="red",size=0.7)
p6<-ggplot(data=protein,aes(Cereals,..density..))+
  geom_histogram(color='white',fill='gray60',binwidth = 3)+
  geom_line(stat='density',colour="red",size=0.7)
p7<-ggplot(data=protein,aes(Starch,..density..))+
  geom_histogram(color='white',fill='gray60',binwidth = 3)+
  geom_line(stat='density',colour="red",size=0.7)
p8<-ggplot(data=protein,aes(Nuts,..density..))+
  geom_histogram(color='white',fill='gray60',binwidth = 3)+
  geom_line(stat='density',colour="red",size=0.7)
p9<-ggplot(data=protein,aes(Fr.Veg,..density..))+
  geom_histogram(color='white',fill='gray60',binwidth = 3)+
  geom_line(stat='density',colour="red",size=0.7)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,nrow=3)

## Order the countries with respect to protein consumption

p1 <- ggplot(protein, aes(x = RedMeat, y = reorder(Country, RedMeat))) + 
  geom_point(size = 6)+ylab("Country")
p2 <- ggplot(protein, aes(x = WhiteMeat, y = reorder(Country, WhiteMeat))) + 
  geom_point(size = 6)+ylab("Country")
p3 <- ggplot(protein, aes(x = Eggs, y = reorder(Country, Eggs))) + 
  geom_point(size = 3)+ylab("Country")
p4 <- ggplot(protein, aes(x = Milk, y = reorder(Country, Milk))) + 
  geom_point(size = 6)+ylab("Country")
p5 <- ggplot(protein, aes(x = Fish, y = reorder(Country, Fish))) + 
  geom_point(size = 6)+ylab("Country")
p6 <- ggplot(protein, aes(x = Cereals, y = reorder(Country, Cereals))) + 
  geom_point(size = 6)+ylab("Country")
p7 <- ggplot(protein, aes(x = Starch, y = reorder(Country, Starch))) + 
  geom_point(size = 3)+ylab("Country")
p8 <- ggplot(protein, aes(x = Nuts, y = reorder(Country, Nuts))) + 
  geom_point(size = 3)+ylab("Country")
p9 <- ggplot(protein, aes(x = Fr.Veg, y = reorder(Country, Fr.Veg))) + 
  geom_point(size = 3)+ylab("Country")
grid.arrange(p3,p7,p8,p9,nrow=2)


########~~~~~~~Principal Component Analysis~~~~~~~~~##################
protein_pca <- prcomp(protein[, -1], scale = TRUE)
round(protein_pca$rotation,3)
summary(protein_pca)

#Pick the number of principal components 
par(mfrow=c(1,1))
plot(protein_pca,main="")
mtext(side=1,"European Protein Principal Components",line=1,font=2)

#Graph the principal components
v_12 <- fviz_pca_var(protein_pca,axes = c(1,2))
v_34 <- fviz_pca_var(protein_pca,axes = c(3,4))
grid.arrange(v_12,v_34,nrow=1)

#Bioplot the individuals and principal components
b_12 <- fviz_pca_biplot(protein_pca,axes=c(1,2))+
  geom_text(label=protein$Country, colour="red")
b_34 <- fviz_pca_biplot(protein_pca,axes=c(3,4))+
  geom_text(label=protein$Country, colour="red")
grid.arrange(b_12,b_34,nrow=1)

########~~~~~~~~~~ Hierarchical clustering~~~~~~~~~~~~~~~~~~~########
# scale the data
protein_new <- scale(protein[,-1])
rownames(protein_new) <- protein$Country
hc <- hclust(dist(protein_new), "single")
plot(hc)
hc2 <- hclust(dist(protein_new), "complete")
plot(hc2)
hc3 <- hclust(dist(protein_new), "average")
plot(hc3)
hc4 <- hclust(dist(protein_new), "ward.D")
plot(hc4)
rect.hclust(hc4,5,border = "red")


#cut the tree into 5 cluster and extract the members and get the mean of each cluster

categ<-cutree(hc4,5)
cluster_1 <- apply(protein_new[which(categ==1),],2,mean)
cluster_2 <- apply(protein_new[which(categ==2),],2,mean)
cluster_3 <- apply(protein_new[which(categ==3),],2,mean)
cluster_4 <- apply(protein_new[which(categ==4),],2,mean)
cluster_5 <- apply(protein_new[which(categ==5),],2,mean)
hc_means <- data.frame(rbind(cluster_1,cluster_2,cluster_3,cluster_4,cluster_5))
rownames(hc_means) <- c(1:5)
round(hc_means,3)

##visualize the clustering by projecting the data onto principals

protein_pca_subset <- data.frame(protein_pca$x[,1:4],cluster=as.factor(categ))
rownames(protein_pca_subset) <- protein$Country
p1 <- ggplot(protein_pca_subset,aes(x=PC1,y=PC2,colour=cluster,label=rownames(protein_pca_subset)))+
  geom_point()+geom_text(size=6.5)+ theme_bw()+theme(panel.grid=element_blank())
p2 <- ggplot(protein_pca_subset,aes(x=PC3,y=PC4,color=cluster,label=rownames(protein_pca_subset)))+
  geom_point()+geom_text(size=6.5)+theme_bw()+theme(panel.grid=element_blank())
grid.arrange(p1,p2,nrow=1)



#plot the ggpairs using PCA 

ggpairs(protein_pca_subset,aes(color=cluster, alpha=0.4))



##########~~~~~k-means clustering~~~~########################################3#
library(cluster)
library(fpc)
#Let us pick the best k for K-means clustering using "wss"
K<-10
wss<-rep(0,K)
for (k in 1:K){
  wss[k] <- sum(kmeans(protein_new,k)$withinss)
}
plot(1:K,wss,typ="b",ylab="Total within cluster sum of squares",xlab="Number of clusters (k)",col="blue",lwd=2)




# use silhouette measure to pick optimal k
#run kmeansrun() from k=1:10 and get the average silhouette value for each k
kmeans_asw <- kmeansruns(protein_new, krange=1:10, criterion="asw")

kmeans_asw$bestk # Average silhoutte maximized at cluster 3

# extract the criteron values.
df.sil <- data.frame(ave.sil=kmeans_asw$crit,k=1:10)
plot(df.sil$k,df.sil$ave.sil,type="b",xlab="k", ylab="Average silhoutte width",col="blue", lwd=2)
abline(v=3,lty=2,col="red",lwd=2)

# we use K=3 to do k-means clustering
protein_kmeans<-kmeans(protein_new,3,nstart = 100,iter.max = 100)

# list of cluster assignments
ind=order(protein_kmeans$cluster)
data.frame(protein$Country[ind],protein_kmeans$cluster[ind])



#Show the centers of k-means
round(protein_kmeans$centers,3)

#Plot the silhouette
sil <- silhouette(protein_kmeans$cluster, dist(protein_new))
rownames(sil) <- protein$Country
summary(sil)
plot(sil,col=2:4,border=NA)

# visualize the clusters by projecting the data on the principal components.
p1 <- fviz_cluster(protein_kmeans, data = protein_new, axes=c(1,2),ellipse.type = "norm") + theme_minimal() + ggtitle("k = 3") 
p2 <- fviz_cluster(protein_kmeans, data = protein_new, axes=c(3,4),ellipse.type = "norm") + theme_minimal() + ggtitle("k = 3") 
grid.arrange(p1,p2,nrow=1)


###########~~mixture model-based clustering~~~~~#########
library(mclust)
mod1 = Mclust(protein_new)
summary(mod1,parameters=TRUE)

#Extract the clusters
protein[which(mod1$classification==1),"Country"]
protein[which(mod1$classification==2),"Country"]
protein[which(mod1$classification==3),"Country"]
protein[which(mod1$classification==4),"Country"]

# Model selection with BIC
fviz_mclust_bic(mod1)


#Project the clustering onto the principal components
p1 <- fviz_mclust(mod1, "classification",geom="point",axes=c(1,2))+geom_text(aes(color=factor(mod1$classification)),label=protein$Country)
p2 <-  fviz_mclust(mod1, "classification",geom="point",axes=c(3,4))+geom_text(aes(color=factor(mod1$classification)),label=protein$Country)
grid.arrange(p1,p2,nrow=1)

#Extract the cluster mean
cluster_mean <- summary(mod1,parameters=TRUE)["mean"]
cluster_mean

