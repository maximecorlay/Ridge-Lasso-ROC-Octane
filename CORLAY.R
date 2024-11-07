## STA203 - TP noté - Maxime Corlay (seul)
rm(list=objects());graphics.off()

###################
### 2- analyse exploratoire
###################
# Q2.1 
df_train = read.table("gasolineTrain.txt",header=TRUE)
df_test = read.table("gasolineTest.txt",header=TRUE)

xtrain = as.matrix(df_train[,-1])
xtest = as.matrix(df_test[,-1])
ytrain = df_train$octane
ytest = df_test$octane

par(mfrow=c(1,1))
boxplot(xtrain)

matplot(t(xtrain),type="l")

library(corrplot)
corrplot(cor(df_train[,-1]), method = "color",tl.pos='n')

x=colnames(df_train[,-1])
G1=x[which(x=="X900.nm"):which(x=="X1134.nm")]
G2=c(x[which(x=="X1134.nm"):which(x=="X1152.nm")],x[which(x=="X1642.nm"):which(x=="X1700.nm")])
G3=x[which(x=="X1154.nm"):which(x=="X1640.nm")]

# Combinaison des groupes de colonnes en un seul ensemble de colonnes
colonnes_a_calculer <- c(G1, G2, G3)

# Calcul de la moyenne des valeurs sur les lignes des colonnes spécifiées
moyenne_lignes <- apply(df_train[colonnes_a_calculer], 1, mean)

c1=apply(df_train[,G1],1,mean)
c2=apply(df_train[,G2],1,mean)
c3=apply(df_train[,G3],1,mean)

df_moy<-data.frame(
  G1=c1,
  G2=c2,
  G3=c3
)
corrplot(cor(df_moy),method="circle")

# Q2.2
library(FactoMineR)
library(ggplot2)

res.pca=PCA(df_train[,-1],ncp=6)
res.pca$ind$octane <- df_train[,1]

# représenter l'éboulis des valeurs propres, et les pourcentages
df_valeurs_propres<-data.frame(
  valeurs_propres=(res.pca$eig)[,1],
  pourcentage_inertie=(res.pca$eig)[,2]
)
pourcentage_inertie<-ggplot(df_valeurs_propres,aes(1:35,pourcentage_inertie))
pourcentage_inertie<-pourcentage_inertie+geom_col(color="green",fill=3)
approx<-round(df_valeurs_propres$pourcentage_inertie,1)
pourcentage_inertie<-pourcentage_inertie+geom_label(aes(label=round(approx,3)))
pourcentage_inertie

# représenter le nuage sur les six premiers axes principaux
df_nuage1<-data.frame(
  dim1=(res.pca$ind$coord)[,1],
  dim2=(res.pca$ind$coord)[,2],
  octane=res.pca$ind$octane
)
df_nuage2<-data.frame(
  dim3=(res.pca$ind$coord)[,3],
  dim4=(res.pca$ind$coord)[,4],
  octane=res.pca$ind$octane
)
df_nuage3<-data.frame(
  dim5=(res.pca$ind$coord)[,5],
  dim6=(res.pca$ind$coord)[,6],
  octane=res.pca$ind$octane
)

g1 <- ggplot(df_nuage1,aes(dim1,dim2,color=octane))
g1 <- g1+geom_point(shape=15,size=3)  # nuage de points sans éparpillement aléatoire
g1 <- g1+scale_color_gradientn(colors = c("blue", "red"), na.value = "grey50") 
g1 <- g1+geom_hline(aes(yintercept=0))
g1 <- g1+geom_vline(aes(xintercept=0))
g1
g2 <- ggplot(df_nuage2,aes(dim3,dim4,color=octane))
g2 <- g2+geom_point(shape=15,size=3)  # nuage de points sans éparpillement aléatoire
g2 <- g2+scale_color_gradientn(colors = c("blue", "red"), na.value = "grey50") 
g2 <- g2+geom_hline(aes(yintercept=0))
g2 <- g2+geom_vline(aes(xintercept=0))
g2
g3 <- ggplot(df_nuage3,aes(dim5,dim6,color=octane))
g3 <- g3+geom_point(shape=15,size=3)  # nuage de points sans éparpillement aléatoire
g3 <- g3+scale_color_gradientn(colors = c("blue", "red"), na.value = "grey50") 
g3 <- g3+geom_hline(aes(yintercept=0))
g3 <- g3+geom_vline(aes(xintercept=0))
g3

library(patchwork)
g<-g1+g2+g3
g

###################
### 3- régression pénalisée
###################
# Q3.1
library(glmnet)

# Estimer le modèle de régression ridge avec la fonction glmnet
grid=10^seq(6,-10,length=100)
x = as.matrix(df_train[,-1])
y = df_train[,1]
ridge.fit = glmnet(x,y,alpha=0,lambda=grid)

# Le recalculer en fonction des estimées des autres paramètres suivant la formule de la section 1.2
intercept.perso<-function(X, Y, k_star) {
  n <- nrow(X)
  p <- ncol(X)
  Id_n <- diag(n)
  Id_p <- diag(p)
  A <- solve(t(X) %*% X + k_star * Id_p)
  theta_0_star <- solve(n-t(rep(1, n)) %*% X %*% A %*% t(X) %*% rep(1, n),t(rep(1, n)) %*% (Id_n - X %*% A %*% t(X)) %*% Y)
  return(theta_0_star)
}
theta_results <- numeric(100)
for (i in c(1:100)) {
  k_star <- grid[i]
  theta_results[i] <- intercept.perso(x, y, k_star)
}
df_grid_intercept<-data.frame(
  grid=seq(6,-10,length=100),
  intercept=coef(ridge.fit)["(Intercept)",],
  theta_star=theta_results
)
g <- ggplot(df_grid_intercept,aes(grid,intercept))
g <- g+geom_point(size=1,color="red")  
g <- g+geom_point(aes(grid,theta_star),color="blue")
g <- g+geom_hline(aes(yintercept=0))
g <- g+geom_vline(aes(xintercept=0))
g <- g+labs(title = "Intercept en fonction de lambda",
            x = "log(Lambda)",
            y = "Intercept")
g <- g + annotate(
  "text",
  x = Inf, y = Inf,
  label = "Calculé avec glmnet",
  hjust = 1, vjust = 4,
  color = "red"
)
g <- g + annotate(
  "text",
  x = Inf, y = Inf,
  label = "Calculé matriciellement avec 1.2",
  hjust = 1, vjust = 2,
  color = "blue"
)
g


# centrons seulement xtrain=X
xc=scale(xtrain,scale=FALSE)
ridge.fit = glmnet(xc,ytrain,alpha=0,lambda=grid)
df_grid_intercept<-data.frame(
  grid=seq(6,-10,length=100),
  intercept=coef(ridge.fit)["(Intercept)",]
)
g <- ggplot(df_grid_intercept,aes(grid,intercept))
g <- g+geom_point(size=1)  # nuage de points sans éparpillement aléatoire
g <- g+scale_color_gradientn(colors = c("blue", "red"), na.value = "grey50") 
g <- g+geom_hline(aes(yintercept=0))
g <- g+geom_hline(aes(yintercept=87))
g <- g+geom_vline(aes(xintercept=0))
g <- g+labs(title = "Intercept en fonction de lambda",
            x = "log(Lambda)",
            y = "Intercept")
g

mean(ytrain)

# centrons aussi ytrain=Y
yc=scale(ytrain,scale=FALSE)
ridge.fit = glmnet(xc,yc,alpha=0,lambda=grid)
df_grid_intercept<-data.frame(
  grid=seq(6,-10,length=100),
  intercept=coef(ridge.fit)["(Intercept)",]
)
g <- ggplot(df_grid_intercept,aes(grid,intercept))
g <- g+geom_point(size=1)  # nuage de points sans éparpillement aléatoire 
g <- g+geom_hline(aes(yintercept=0))
g <- g+geom_vline(aes(xintercept=0))
g <- g+labs(title = "Intercept en fonction de lambda",
            x = "log(Lambda)",
            y = "Intercept")
g

# réduisons xtrain et ytrain
xcr=scale(xtrain)
ycr=scale(ytrain)
res.svd <- svd(xcr)   # décomposition en valeurs singulières
D=1/(res.svd$d)[-36]
D[36]=0
U=res.svd$u
V=res.svd$v
DD=diag(D)
A0=V%*%DD%*%t(U) # matrice A0

# vérification que A0 est une pseudo-inverse de x
M=A0%*%xcr%*%A0-A0
N=xcr%*%A0%*%xcr-xcr
all(M<1e-10)  # true
all(N<1e-10)  # true

# coefficients limites
theta_hat_lim=A0 %*% ycr
print(theta_hat_lim)

# Q3.2 
# dans cette question, on suppose les données standardisées
# on va utiliser xcr et ycr (xtrain et ytrain centrées réduites)
# mais lm.ridge exige un dataframe en entrée donc on standardise aussi df_train
dfcr=scale(df_train)

# bibliothèques et grille de départ
library(MASS)
library(glmnet)
grid=10^seq(6,-10,length=100)

# comparaison de lm.ridge et de glmnet
ridge.fit.lm = lm.ridge(octane~.,data=as.data.frame(dfcr),lambda=36*grid)
ridge.fit.glmnet = glmnet(xcr,ycr,alpha=0,lambda=grid)
ridge.fit.lm.coef=ridge.fit.lm$coef
ridge.fit.glmnet.coef=as.matrix(ridge.fit.glmnet$beta)
err=abs(ridge.fit.lm.coef - ridge.fit.glmnet.coef) / abs(ridge.fit.lm.coef)   # erreur relative
dferr <- expand.grid(lignes = 1:nrow(err), colonnes = 1:ncol(err))
dferr$val <- as.vector(err)
g<-ggplot(dferr, aes(x = lignes, y = colonnes, fill = cut(val, breaks = c(0,1,10,100,1000,Inf))))
g<-g+geom_tile() 
g<-g+scale_fill_manual(values = c("grey","blue", "green", "orange", "red")  , name = "Err. relative", 
                       labels = c("<1","1<...<10", "10<...<100", "100<...<1000", ">1000"))
g<-g+theme_minimal() 
g<-g+labs(title = "Erreurs relatives des coefficients entre glmnet et lm.ridge")
g
# je ne sais pas si c'est le meilleur moyen de comparer deux matrices de coefficients lorsque
# les coefficients sont aussi faibles, mais en tout cas, on constate que les résultats 
# des deux fonctions sont très proches. Juste pour s'en convaincre :
(ridge.fit.lm.coef)[1,1]  # -1.968386e-07
(ridge.fit.glmnet.coef)[1,1] # -1.969142e-07

# maintenant on essaie de calculer matriciellement l'optimum (en utilisant la Q1.2)
perso<-function(x,y,lambda){
  n<-nrow(x)
  p<-ncol(x)
  Id_n<-diag(n)
  Id_p<-diag(p)
  return(solve(t(x)%*%x+lambda*Id_p,t(x)%*%y))
}
ridge.fit.perso.coef<-matrix(0, nrow = 401, ncol = 100)
for (i in c(1:100)) {
  lambda<-grid[i]
  ridge.fit.perso.coef[,i]<-perso(xcr,ycr,36*lambda)  # multiplier par 36 pour comparer à glmnet
}
(ridge.fit.perso.coef)[1,1]  # -1.940875e-07
# ridge.fit.perso.coef est cohérent avec ridge.fit.lm.coef et ridge.fit.glmnet.coef
# Pour s'en convaincre de manière plus large, je propose de comparer l'erreur relative avec glmnet
err2=abs(ridge.fit.perso.coef - ridge.fit.glmnet.coef) / abs(ridge.fit.perso.coef)   # erreur relative
dferr2 <- expand.grid(lignes = 1:nrow(err2), colonnes = 1:ncol(err2))
dferr2$val <- as.vector(err2)
g<-ggplot(dferr2, aes(x = lignes, y = colonnes, fill = cut(val, breaks = c(0,1,10,100,1000,Inf))))
g<-g+geom_tile() 
g<-g+scale_fill_manual(values = c("grey","blue", "green", "orange", "red"), name = "Err. relative", 
                       labels = c("<1","1<...<10", "10<...<100", "100<...<1000", ">1000"))
g<-g+theme_minimal() 
g<-g+labs(title = "Erreurs relatives des coefficients entre glmnet et Q1.2")
g

# Enfin, pour conclure cette question, je propose de comparer la valeur de la fonction objectif
M=ridge.fit.lm.coef
N=ridge.fit.glmnet.coef
O=ridge.fit.perso.coef
sum((ycr-xcr%*%M[,1])^2)+36*grid[1]*sum(M^2) # 1532657938
sum((ycr-xcr%*%N[,1])^2)+36*grid[1]*sum(N^2) # 1706351910
sum((ycr-xcr%*%O[,1])^2)+36*grid[1]*sum(O^2) # 1574033159
objM=rep(1,100)
objN=rep(1,100)
objO=rep(1,100)
for(i in c(1:100)){
  objM[i]=log(sum((ycr-xcr%*%M[,i])^2)+36*grid[i]*sum(M^2))
  objN[i]=log(sum((ycr-xcr%*%N[,i])^2)+36*grid[i]*sum(N^2))
  objO[i]=log(sum((ycr-xcr%*%O[,i])^2)+36*grid[i]*sum(O^2))
}
dfobj<-data.frame(
  x = 1:100,
  y = c(objM, objN, objO),
  suite = rep(c("lm", "glmnet", "perso"), each = 100)
)
g<-ggplot(data=dfobj,aes(x = x, y = y, color = suite))
g<-g+geom_line() 
g<-g+labs(title = "Valeurs du logarithme népérien des fonctions objectif",x = "Index de la grid",y = "Log de la valeur objectif",color="Méthode")
g<-g+theme_minimal()
g

# Q3.3
# sélection des indices de validation croisée
library(pls)
set.seed(87)  # initialisation de la graine
cvseg=cvsegments(nrow(df_train),k=4,type="random")  # création des quatre folds

# grille adaptée à glmnet 
gridglmnet=seq(0,1,0.001)

# validation croisée avec cv.glmnet
res.cv.glmnet=cv.glmnet(x=xcr,y=ycr,lambda=gridglmnet,type.measure="deviance",folds=4)
res.cv.glmnet # par défaut type.measure est RSS
# Le lambda optimal est 0.011

# MSE en fonction de lambda
plot(res.cv.glmnet)

# maintenant, nous allons réobtenir ce résultat de manière indépendante 
# nous utilisons glmnet, mais pas cv.glmnet
errcross=matrix(0, nrow = length(gridglmnet), ncol = 4)
for (i in 1:length(gridglmnet)) {
  for (j in 1:4) {
    val=dfcr[unlist(cvseg[j]),]  # ens de validation (de longueur 9)
    app=dfcr[-unlist(cvseg[j]),]  # ens d'apprentissage (de longueur 27)
    ridge.fit.glmnet2=glmnet(app[,-1],app[,1],lambda=gridglmnet[i])
    predictions=predict(ridge.fit.glmnet2,newx = val[,-1])
    errcross[i, j]<-sum((val[,1]-predictions)^2)/9
  }
}
errmoy<-apply(errcross,1, mean)
gridglmnet[which.min(errmoy)] # Lambda=0.014

# MSE en fonction de lambda
sd_values=apply(errcross,1,sd)
plot(errmoy, type = "l", ylim = c(min(errmoy-1.96*sd_values/sqrt(1001)),max(errmoy+1.96*sd_values/sqrt(1001))),xlab="Index de la grille",ylab="MSE",main="MSE=f(lambda)")
lines(errmoy+1.96*sd_values/sqrt(1001), col = "red")
lines(errmoy-1.96*sd_values/sqrt(1001), col = "red")

sd_values=apply(errcross,1,sd)
dfmse=data.frame(
  lambda=gridglmnet,
  min=errmoy-1.96*sd_values/sqrt(1001),
  moy=errmoy,
  max=errmoy+1.96*sd_values/sqrt(1001)
)
g<-ggplot(data=dfmse)
g<-g+geom_errorbar(aes(x=lambda,ymin=min,ymax=max),color="grey")
g<-g+geom_line(aes(x=lambda,y=moy),color="black")
g<-g+geom_vline(aes(xintercept=0.013))
g<-g+labs(title="MSE=f(lambda)")
g

# réajustement sur tout xtrain
ridge.fit.glmnet3=glmnet(xcr,ycr,lambda=0.014)

# calcul de l'erreur de généralisation
xtcr=scale(xtest)
ytcr=scale(ytest)
predictions3=predict(ridge.fit.glmnet3,newx=xtest)
mean((ytcr-predictions3)^2) # MSE de généralisation = 0.9610059

# bonus
# Comparaison des régressions Ridge et Lasso
ridge.fit.glmnet = glmnet(xcr,ycr,alpha=0,lambda=seq(0,1,0.01)/36)
plot(ridge.fit.glmnet)
lasso.fit.glmnet = glmnet(xcr,ycr,alpha=1,lambda=seq(0,1,0.01)/36)
plot(lasso.fit.glmnet)

###################
### 4- régression logistique pénalisée
###################

# Q4.1
# on utilise la fonction sigmoide, fct réciroque de la fonction logit
z=1/(1+exp(88-ytrain))
ztest=1/(1+exp(88-ytest))

# Q4.2
library(glmnet)
set.seed(42)  
foldid <- sample(rep(1:4, length.out = 36))
cv.glmnet(xcr,z, alpha = 0, family = "gaussian", foldid = foldid) # Lambda=2.077
cv.glmnet(xcr,z, alpha = 1, family = "gaussian", foldid = foldid) # Lambda=0.006055

# réajustement sur tout xtrain
ridge.fit.glmnet4=glmnet(xcr,z,alpha=0,lambda=2.077)
lasso.fit.glmnet4=glmnet(xcr,z,alpha=1,lambda=0.006055)

# calcul de l'erreur de généralisation
ridge.predictions4=predict(ridge.fit.glmnet4,newx=xtcr)
mean((ztest-ridge.predictions4)^2) # MSE de généralisation ridge = 0.0123163
lasso.predictions4=predict(lasso.fit.glmnet4,newx=xtcr)
mean((ztest-lasso.predictions4)^2) # MSE de généralisation ridge = 0.008750137

# Les méthodes Ridge et Lasso donnent des erreurs de généralisation équivalentes.

