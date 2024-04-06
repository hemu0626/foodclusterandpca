###################################################
#This is for final project

rm(list=ls())
library(ggplot2)
library(GGally)
library(factoextra)
library(gridExtra)
library(corrplot)
library(MASS)
food<- read.csv('D:/FYP/food.csv')
head(food)
str(food)
#######~~Descriptive analysis~~~~~~###################~~~~~~~
summary(food)
nrow(food)

#Removing the missing values
food <- na.omit(food)
set.seed(20191126)

#calulate the correlation coefficient
cor <- cor(food[,-1])
round(cor, 2)

#correlation plot
corrplot(cor, type = "upper", order = "AOE", 
         tl.col = "black", tl.srt = 45,addCoef.col = "gray")
#find the variables with strong correlation
symnum(cor)

#density plot
p1<-ggplot(data=food,aes(Fried,..density..))+
  geom_histogram(color='white',fill='gray60',binwidth = 0.01)+ 
  geom_line(stat='density',colour="red",size=0.5)
p2<-ggplot(data=food,aes(Grilled,..density..))+
  geom_histogram(color='white',fill='gray60',binwidth = 0.01)+
  geom_line(stat='density',colour="red",size=0.5)
p3<-ggplot(data=food,aes(Spicy,..density..))+
  geom_histogram(color='white',fill='gray60',binwidth = 0.01)+
  geom_line(stat='density',colour="red",size=0.5)
p4<-ggplot(data=food,aes(Tongue.numbing ,..density..))+
  geom_histogram(color='white',fill='gray60',binwidth = 0.01)+
  geom_line(stat='density',colour="red",size=0.5)
p5<-ggplot(data=food,aes(Instant.boiled,..density..))+
  geom_histogram(color='white',fill='gray60',binwidth = 0.01)+
  geom_line(stat='density',colour="red",size=0.5)
p6<-ggplot(data=food,aes(Sweet,..density..))+
  geom_histogram(color='white',fill='gray60',binwidth = 0.01)+
  geom_line(stat='density',colour="red",size=0.5)
p7<-ggplot(data=food,aes(Diabetes,..density..))+
  geom_histogram(color='white',fill='gray60',binwidth = 0.01)+
  geom_line(stat='density',colour="red",size=0.7)
p8<-ggplot(data=food,aes(Hypertension,..density..))+
  geom_histogram(color='white',fill='gray60',binwidth = 0.1)+
  geom_line(stat='density',colour="red",size=0.7)
p9<-ggplot(data=food,aes(FPG,..density..))+
  geom_histogram(color='white',fill='gray60',binwidth = 0.5)+
  geom_line(stat='density',colour="red",size=0.7)
p10<-ggplot(data=food,aes(PPG,..density..))+
  geom_histogram(color='white',fill='gray60',binwidth = 0.5)+
  geom_line(stat='density',colour="red",size=0.7)
p11<-ggplot(data=food,aes(BMI,..density..))+
  geom_histogram(color='white',fill='gray60',binwidth = 0.5)+
  geom_line(stat='density',colour="red",size=0.7)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,nrow=4)

#find the maximum and minimum value in each column 
apply(food[,-1],2,min)
apply(food[,-1],2,max)
## Order the provinces with respect to dietary preference

p1 <- ggplot(food, aes(x = Fried, y = reorder(Provinces, Fried))) + 
  geom_point(size = 6)+ylab("Provinces")
p2 <- ggplot(food, aes(x = Grilled, y = reorder(Provinces, Grilled))) + 
  geom_point(size = 6)+ylab("Provinces")
p3 <- ggplot(food, aes(x = Spicy, y = reorder(Provinces, Spicy))) + 
  geom_point(size = 3)+ylab("Provinces")
p4 <- ggplot(food, aes(x = Tongue.numbing, y = reorder(Provinces, Tongue.numbing))) + 
  geom_point(size = 6)+ylab("Provinces")
p5 <- ggplot(food, aes(x = Instant.boiled, y = reorder(Provinces, Instant.boiled))) + 
  geom_point(size = 6)+ylab("Provinces")
p6 <- ggplot(food, aes(x = Sweet, y = reorder(Provinces, Sweet))) + 
  geom_point(size = 6)+ylab("Provinces")
p7 <- ggplot(food, aes(x = Diabetes, y = reorder(Provinces, Diabetes))) + 
  geom_point(size = 3)+ylab("Provinces")
p8 <- ggplot(food, aes(x = Hypertension, y = reorder(Provinces, Hypertension))) + 
  geom_point(size = 3)+ylab("Provinces")

grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,nrow=2)

########~~~~~~~Principal Component Analysis~~~~~~~~~##################
library(FactoMineR)
food.pca2<-PCA(food[, -1],scale.unit=T,ncp=5,graph=T)
print(food.pca2)
summary(food.pca2)

#Pick the number of principal components 
get_eigenvalue(food.pca2)
fviz_eig(food.pca2,addlabels = TRUE)

#visualize cos2 
corrplot(get_pca_var(food.pca2)$cos2,is.corr = FALSE)
fviz_cos2(food.pca2,choice = "var",axes = 1:2)
##Graph the principal components
fviz_pca_var(food.pca2,axes = c(1,2),
             col.var="cos2",
             gradient.cols=c("#00AFBB","#E7B800","#FC4E07"),
             repel=TRUE)# Avoid text overlapping

#Bioplot the individuals and principal components
#visualze sample scores
get_pca_ind(food.pca2)
food.pca2$ind
b_12 <-fviz_pca_ind(food.pca2)
b_34 <- fviz_pca_biplot(food.pca2,axes=c(1,2))+
  geom_text(label=food$Provinces, colour="red")
grid.arrange(b_12,b_34,nrow=1)

########~~~~~~~~~~ Hierarchical clustering~~~~~~~~~~~~~~~~~~~########
# scale the data
food_new <- scale(food[,-1])
rownames(food_new) <- food$Provinces
hc <- hclust(dist(food_new), "single")
plot(hc)
hc2 <- hclust(dist(food_new), "complete")
plot(hc2)
hc3 <- hclust(dist(food_new), "average")
plot(hc3)
hc4 <- hclust(dist(food_new), "ward.D")
plot(hc4)
rect.hclust(hc4,5,border = "red")

#cut the tree into 5 cluster and extract the members and get the mean of each cluster

categ<-cutree(hc4,5)
cluster_1 <- apply(food_new[which(categ==1),],2,mean)
cluster_2 <- apply(food_new[which(categ==2),],2,mean)
cluster_3 <- apply(food_new[which(categ==3),],2,mean)
cluster_4 <- apply(food_new[which(categ==4),],2,mean)
cluster_5 <- apply(food_new[which(categ==5),],2,mean)
hc_means <- data.frame(rbind(cluster_1,cluster_2,cluster_3,cluster_4,cluster_5))
rownames(hc_means) <- c(1:5)
round(hc_means,3)

##visualize the clustering by projecting the data onto principals
food_pca <- prcomp(food[, -1], scale = TRUE)
round(food_pca$rotation,3)
food_pca_subset <- data.frame(food_pca$x[,1:4],cluster=as.factor(categ))
rownames(food_pca_subset) <- food$Provinces
p1 <- ggplot(food_pca_subset,aes(x=PC1,y=PC2,colour=cluster,label=rownames(food_pca_subset)))+
  geom_point()+geom_text(size=6.5)+ theme_bw()+theme(panel.grid=element_blank())
grid.arrange(p1,nrow=1)

#plot the ggpairs using PCA 

ggpairs(food_pca_subset,aes(color=cluster, alpha=0.4))

##########~~~~~k-means clustering~~~~########################################3#
library(cluster)
library(fpc)
#Let us pick the best k for K-means clustering using "wss"
K<-10
wss<-rep(0,K)
for (k in 1:K){
  wss[k] <- sum(kmeans(food_new,k)$withinss)
}
plot(1:K,wss,typ="b",ylab="Total within cluster sum of squares",xlab="Number of clusters (k)",col="blue",lwd=2)

# use silhouette measure to pick optimal k
#run kmeansrun() from k=1:10 and get the average silhouette value for each k
kmeans_asw <- kmeansruns(food_new, krange=1:10, criterion="asw")

kmeans_asw$bestk # Average silhoutte maximized at cluster 3

# extract the criteron values.
df.sil <- data.frame(ave.sil=kmeans_asw$crit,k=1:10)
plot(df.sil$k,df.sil$ave.sil,type="b",xlab="k", ylab="Average silhoutte width",col="blue", lwd=2)
abline(v=3,lty=2,col="red",lwd=2)

# we use K=3 to do k-means clustering
food_kmeans<-kmeans(food_new,3,nstart = 100,iter.max = 100)

# list of cluster assignments
ind=order(food_kmeans$cluster)
data.frame(food$Provinces[ind],food_kmeans$cluster[ind])

#Show the centers of k-means
round(food_kmeans$centers,3)

#Plot the silhouette
sil <- silhouette(food_kmeans$cluster, dist(food_new))
rownames(sil) <- food$Country
summary(sil)
plot(sil,col=2:4,border=NA)

# visualize the clusters by projecting the data on the principal components.
p1 <- fviz_cluster(food_kmeans,data = food_new, axes=c(1,2),ellipse.type = "norm") + theme_minimal() + ggtitle("k = 3") 
grid.arrange(p1,nrow=1)

###########~~mixture model-based clustering~~~~~#########
library(mclust)
mod1 = Mclust(food_new,G=3)
summary(mod1,parameters=TRUE)
?Mclust
#Extract the clusters
food[which(mod1$classification==1),"Provinces"]
food[which(mod1$classification==2),"Provinces"]
food[which(mod1$classification==3),"Provinces"]
food[which(mod1$classification==4),"Provinces"]

# Model selection with BIC
fviz_mclust_bic(mod1)

#Project the clustering onto the principal components
p1 <- fviz_mclust(mod1, "classification",geom="point",axes=c(1,2))+geom_text(aes(color=factor(mod1$classification)),label=food$Provinces)
p2 <- fviz_mclust(mod1, "classification",geom="point",axes=c(3,4))+geom_text(aes(color=factor(mod1$classification)),label=food$Provinces)
grid.arrange(p1,p2,nrow=1)

#Extract the cluster mean
cluster_mean <- summary(mod1,parameters=TRUE)["mean"]
cluster_mean

### fuzzy-cmeans~~~~~#########
library(e1071)
library(factoextra)
set.seed(123)
df <- scale(food[,-1])
# Compute fuzzy clustering
cm <- cmeans(df, 4)
cm
# Visualize using corrplot
library(corrplot)
colnames(cm$membership)<-c('PC1','PC2','PC3','PC4')
rownames(cm$membership)<-food$Provinces
corrplot(cm$membership,is.corr = FALSE)
cm$membership
# Observation groups/clusters
cm$cluster
#visualize clusters
fviz_cluster(list(data = food_new, cluster=cm$cluster), 
             ellipse.type = "norm",
             ellipse.level = 0.68,
             palette = "jco",
             ggtheme = theme_minimal())

