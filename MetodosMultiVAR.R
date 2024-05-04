datos <- as.data.frame(state.x77)

datos

#cambiar nombres de variables con espacios:
colnames(datos)[4] <- "Life.Exp"
colnames(datos)[6] <- "HS.Grand"

pca_lic <- prcomp(datos, scale = TRUE)

pca_lic

#Dimensiones n y p
dim(datos)
n <- dim(datos)[1]
p <- dim(datos)[2]

#Gráfico de dispersión 
pairs(datos, col="brown", pch=19, main="Datos State.x77")

#Vector de medias y matriz de covarianzas
vmedias <- colMeans(datos)
vmedias

mcov <- cov(datos)
mcov

#Obtencion de los componentes principales
  #autovalores y matriz de autovalores
AutoV <- eigen(mcov)
AutoV

  #para guardarlos
eigenval <- AutoV$values #autovalores
eigenval
eigenvec <- AutoV$vectors #autovectores
eigenvec

  #proporcion de variabilidad:
prop.var <- eigenval / sum(eigenval) 
prop.var
  #proporcion de variabilidad acumulada 
prop.var.acum <- cumsum(eigenval) /sum(eigenval)
prop.var.acum


# obtener los componetes a traves de la matriz de correlacion 

mcorr <- cor(datos)
eigencorr <- eigen(mcorr)
eigencorr

eigen.val <- eigencorr$values 
eigen.vec <- eigencorr$vectors

  #proporcion de variabilidad
prop.vari <- eigen.val / sum(eigen.val)
prop.vari.acum <- cumsum(eigen.val) / sum(eigen.val)
prop.vari.acum

  #media de los eigen
mean(eigen.val)

#obtener los scores, las nuevas variables

ones <- matrix(rep(1,n), nrow = n, ncol=1)  #vector de unos
datos.cen  <- as.matrix(datos) - ones %*% vmedias
datos.cen

Dx <- diag(diag(mcov))
Dx

y <- datos.cen %*% solve(Dx)^(1/2)
y #datos estandarizados

scores <- y%*% eigen.vec
colnames(scores) <- c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8")
scores

  #graficados
pairs(scores, main="scores", col="darkblue",pch=19)

#grafica screeplot
screeplot(princomp(datos, cor = T), main="screeplot", col="darkgreen", type = "lines",pch=19)

#graficar los componentes

plot(scores[,1:2], xlab = "Primera componente", ylab = "Segunda componente", col="darkblue"
       , pch=19, main="primera y segunda componente principal")
text(scores[,1:2],labels = rownames(datos), pos=1, col= "darkblue")

plot(scores[,c(1,3)], xlab = "Primera componente", ylab = "tercera componente", col="darkblue"
       , pch=19, main="primera y tercera componente principal")
text(scores[,c(1,3)],labels = rownames(datos), pos=1, col= "darkblue")


plot(scores[,2:3], xlab = "Segundaa componente", ylab = "tercera componente", col="darkblue"
       , pch=19, main="segunda y tercera componente principal")
text(scores[,2:3],labels = rownames(datos), pos=1, col= "darkblue")

#correlacion de los componentes que elegimos y las variables

corr <- diag(eigen.val[1:3]^(1/2)) %*% t(eigen.vec[,1:3])
corr


#metodo lic___________________________________________________________________________

#install.packages("ISLR")
library(ISLR)

names(NCI60)

datos_nci <- NCI60$data
dim(datos_nci)

head(datos_nci)[, 1:6]
# Tipos de cáncer distintos en el set de datos
unique(NCI60$labs)
# Número de muestras por tipo de cáncer
table(NCI60$labs)

# Media de la expresión de cada gen (muestra de los 10 primeros). 
# (MARGIN = 2 para que se aplique la función a las columnas)
apply(X = datos_nci, MARGIN = 2, FUN = mean)[1:10]

# Varianza de la expresión de cada gen (muestra de los 10 primeros)
apply(X = datos_nci, MARGIN = 2, FUN = var)[1:10]

pca_nci <- prcomp(datos_nci, scale = TRUE)

pca_nci

names(pca_nci)

# Muestra de los primeros 6 elementos del vector de loadings de los 5 primeros componentes
head(pca_nci$rotation)[, 1:5]
dim(pca_nci$rotation)

head(pca_nci$x)[,1:5]

pca_nci$sdev

# Varianza explicada por cada componente
pca_nci$sdev^2

summary(pca_nci)


#con la funcion PCA_______________________________________________________________________

#install.packages("FactoMineR")
library(FactoMineR)

pca2.nci <- PCA(X = datos_nci, scale.unit = TRUE, ncp = 64, graph = FALSE)

print(pca2.nci)

head(pca2.nci$eig)

#representacion grafica_________________________________________________________________
#install.packages("factoextra")
library(factoextra)

fviz_pca_ind(pca_nci, geom.ind = "point", 
             col.ind = "#FC4E07", 
             axes = c(1, 2), 
             pointsize = 1.5) 

colores <- function(vec){
  # la función rainbow() devuelve un vector que contiene el número de colores distintos
  col <- rainbow(length(unique(vec)))
  return(col[as.numeric(as.factor(vec))])
}

par(mfrow = c(1,2))
# Observaciones sobre PC1 y PC2
plot(pca_nci$x[,1:2], col = colores(NCI60$labs), 
     pch = 19, 
     xlab = "Z1", 
     ylab = "Z2")
# Observaciones sobre PC1 y PC3
plot(pca_nci$x[,c(1, 3)], col = colores(NCI60$labs), 
     pch = 19, 
     xlab = "Z1", 
     ylab = "Z3")

fviz_pca_var(pca_nci, col.var = "cos2", 
             geom.var = "arrow", 
             labelsize = 2, 
             repel = FALSE)
var <- get_pca_var(pca_nci)
var

# Datos de expresion de 5 genes
datos_nci9 <- NCI60$data[, 1:9]

# PCA
pca_nci9 <- prcomp(datos_nci9, scale = TRUE)

biplot(pca_nci9, scale = 0, cex = 0.5, col = c("dodgerblue3", "deeppink3"))


# Normalizamos los scores al rango [0, 1]
pca_nci9$x[, c(1:2)] <- pca_nci9$x[, c(1:2)]/sqrt(nrow(datos_nci9)-1)

biplot(pca_nci9, scale = 0, cex = 0.5, col = c("dodgerblue3", "deeppink3"))

#grafico de codos que mas bonito es
fviz_screeplot(pca_nci, addlabels = TRUE, ylim = c(0, 20), barcolor="lightblue")
# Top 10 variables que más contribuyen a PC1
fviz_contrib(pca_nci, choice = "var", axes = 1, top = 10)


#var explicada en porcentaje
PVE <- 100*pca_nci$sdev^2/sum(pca_nci$sdev^2)
PVE


par(mfrow = c(1,2))

plot(PVE, type = "o", 
     ylab = "PVE", 
     xlab = "Componente principal", 
     col = "darkblue")
plot(cumsum(PVE), type = "o", 
     ylab = "PVE acumulada", 
     xlab = "Componente principal", 
     col = "brown3")

#    ANALISIS FACTORIAL________________________________________________________________

#install.packages("polycor")
#install.packages("ggcorrplot")
library(psych)
library(polycor)
library(ggcorrplot)

bfi_s <- bfi[1:200,1:25] # subconjunto de datos
mat_cor <- hetcor(bfi_s)$correlations #matriz de correlación policorica
ggcorrplot(mat_cor,type="lower",hc.order = T)

#prueba para saber si la matriz es factorizable 
cortest.bartlett(mat_cor)->p_esf
p_esf$p #se rechaza h0 sí estan correlacionados

KMO(mat_cor)


### prueba de dos modelos con tres factores
modelo1<-fa(mat_cor,
            nfactors = 3,
            rotate = "none",
            fm="mle") # modelo máxima verosimilitud

modelo2<-fa(mat_cor,
            nfactors = 3,
            rotate = "none",
            fm="minres") # modelo minimo residuo
######comparando las comunalidades
sort(modelo1$communality,decreasing = T)->c1
sort(modelo2$communality,decreasing = T)->c2
head(cbind(c1,c2))

sort(modelo1$uniquenesses,decreasing = T)->u1
sort(modelo2$uniquenesses,decreasing = T)->u2
head(cbind(u1,u2))

scree(mat_cor)

fa.parallel(mat_cor,n.obs=200,fa="fa",fm="minres")

#Rotaciones
library(GPArotation)
rot<-c("none", "varimax", "quartimax","Promax")
bi_mod<-function(tipo){
  biplot.psych(fa(bfi_s,nfactors = 2,fm="minres",rotate = tipo),main = paste("Biplot con rotación ",tipo),col=c(2,3,4),pch = c(21,18),group = bfi[,"gender"])  
}
sapply(rot,bi_mod)


#interpretacion
modelo_varimax<-fa(mat_cor,nfactors = 5,rotate = "varimax",
                   fa="minres")
fa.diagram(modelo_varimax)

print(modelo_varimax$loadings,cut=0) 


#ESCALADO MULTIDIMENSIONAL _________________________________________________________________
#distancias entre ciudades

data.dist <- eurodist
data.dist

#cantidad de ciudades
n<- nrow(data.dist)

#escalado multidimensional clasico

mds.cities <- cmdscale(data.dist, eig = T)
#acá están los autovalores mds.cities$eig

#grafico de autovalores
plot(mds.cities$eig, pch=19, col="darkblue", xlab="Number", ylab="eigen values", type="o")
abline(a=0,b=0, col="darkred")
#hay autovalores negativos, consideramos ahí tambien la camtidad de coordenadas principales

#medida de presicion 

mp <- sum(abs(mds.cities$eig[1:2])) / sum(abs(mds.cities$eig))
mp

#obteniendo las coordenadas con k=2
mds.cities <- cmdscale(data.dist, eig=TRUE, k=2)

x1 <- mds.cities$points[,1]
x2 <- mds.cities$points[,2]

plot(x1,x2,pch=19,col="darkblue",xlim=range(x1)+c(0,600))
text(mds.cities$points,pos=4,labels=rownames(mds.cities$points),col="red")

x2 <- -x2

plot(x1,x2,pch=19,col="blue",xlim=range(x1)+c(0,600))
text(x1,x2,pos=4,labels=rownames(data.dist),col="red")
  

#escalado no metrico

rm(list=ls()) #borramos los resultados del ej anterior

library(MASS)
#install.packages("HSAUR2")
library(HSAUR2)


data("voting", package = "HSAUR2")

voting.mds <- isoMDS(voting) #escalado no metrico
voting.mds$points

plot(voting.mds$points[,1],voting.mds$points[,2], main="escalamiento no metrico",
     pch=19, col="darkblue", xlab = "primera coordenada", ylab = "segunda coordenada")
text(voting.mds$points, labels = rownames(voting.mds$points), pos=1, col="darkred")


# cluster k medias_______________________________________________________











