#VISUALIZACION DATOS INICIALES

d = 15 #primeros dias a mirar
gd = datos_gd[,1:d]

library(reshape2)
library(ggplot2)
df = data.frame(cbind(as.matrix(1:d),t(gd[1:5,]))) #indicar los indices de las observaciones a graficar
meltdf = melt(df,id="V1")
ggplot(meltdf,aes(x=V1,y=value,colour=variable,group=variable))+ ggtitle(paste("GD en primeros ",d," dias")) + geom_line() + facet_wrap(~ variable, scales = 'free_y', ncol = 1)


#===============================================================================================
#SUAVIZAMIENTO DE LOS DATOS

n = dim(datos_gd)[1]
d = 15 #dias a mirar de las observaciones
gd = datos_gd[,1:d] #datos_gd es una matriz que tiene cada cultivo como fila y GD como columnas
gd = t(gd)

dayrng = c(1,d)
daytime = 1:d

knots = seq(1,d,2)
norder = 8
nbasis = length(knots) + norder - 2
bbasis = create.bspline.basis(dayrng,nbasis,norder,knots) #crear base funcional tipo B-Splines


#escoger parametro de suavizacion
curv.Lfd = int2Lfd(2)
curv.fdPar = fdPar(bbasis,curv.Lfd,10)
lambdas = 10^seq(-4,4,by=0.5)
mean.gcv = rep(0,length(lambdas))

for(i in 1:length(lambdas)){
  curv.fdPari = curv.fdPar
  curv.fdPari$lambda = lambdas[i]
  gdSmoothi = smooth.basis(daytime,gd,curv.fdPari)
  mean.gcv[i] = mean(gdSmoothi$gcv)
}
plot(lambdas,mean.gcv,type='b',log='x') 

lambda = lambdas[which.min(mean.gcv)]
curv.fdPar$lambda = lambda
gdSmooth = smooth.basis(daytime,gd,curv.fdPar)
gdfd = gdSmooth$fd

plot(gdfd,lty =1) #graficar todos los datos funcionales
plotfit.fd(gd,daytime,gdfd) #graficar cada suavizacion 

#===============================================================================================
#REGRESIONES
#===============================================================================================
dias = c(15,21,25,35,41,49,59,65)
mse_FR = 1:length(dias)
mse_PCR = 1:length(dias)
mse_MVR = 1:length(dias)

#===============================================================================================
#REGRESION FUNCIONAL POR BASES FUNCIONALES Y PENALIZACION

library(caret)
flds = createFolds(1:n, k = 10, list = TRUE, returnTrain = FALSE) #para d=59,65 tomar n = n-1 y quitar observacion 203

d = 15
dayrng = c(1,d)
gdfd = gdfd15

daytime = 1:d
knots = seq(1,d,2)
norder = 4
nbasis = length(knots) + norder - 2
betabasis = create.bspline.basis(dayrng,nbasis,norder,knots)
lambdas = c(10,50,100,250,500,750,1000,1500,2000,5000,10000)
curv.Lfd = int2Lfd(2)
mse_lambdas = 1:length(lambdas)

#calibracion del parametro de suavizacion
for (i in 1:length(lambdas)){ #va variando el lambda
  bwtlist = list(len=2)
  cbasis = create.constant.basis(dayrng)
  bwtlist[[1]] = fdPar(cbasis)
  beta.fdPar  = fdPar(betabasis,curv.Lfd,lambdas[i])
  bwtlist[[2]] = beta.fdPar
  
  mses = 1:10
  for (k in 1:10){
    test = flds[[k]]
    train = which(!(1:n %in% test))
    xlist = list(len=2)
    xlist[[1]] = rep(1,length(train))
    xlist[[2]] = gdfd[train]
    
    modelk = fRegress(corte.bloque[train],xlist,bwtlist)
    beta = modelk$betaestlist[[2]]$fd
    intercept = modelk$betaestlist[[1]]$fd$coefs
    xmean = mean(gdfd[train])
    ymean = mean(corte.bloque[train])
    pred = replicate(length(test),ymean-inprod(xmean,beta))+inprod(gdfd[test],beta)
    mses[k]=mean((pred-corte.bloque[test])^2)
  }
  mse_lambdas[i]=mean(mses)
}
plot(lambdas,mse_lambdas,type='b',col="red",xlab = "lambda penalizacion",ylab="MSE por 10-folds") #graficar parametro de suavizacion y MSE

bwtlist[[2]]$lambda = lambdas[which.min(mse_lambdas)]
xlist[[1]] = rep(1,n)
xlist[[2]] = gdfd
model = fRegress(corte.bloque,xlist,bwtlist)

plot(model$betaestlist[[2]], main = paste("Funcion Beta Estimada para",d,"dias") , col = "blue") #graficar funcion de regresion estimada

sstot = sum((corte.bloque.trim-mean(corte.bloque.trim))^2)
ssres = sum((corte.bloque.trim-model$yhatfdobj)^2)
1 - ssres/sstot #calculo de R^2
min(mse_lambdas) #MSE por 10-folds

mse_FR[which(dias == d)] = min(mse_lambdas)
#===============================================================================================
#REGRESION FUNCIONAL POR COMPONENTES PRINCIPALES

d=15
gdfd = gdfd15

dayrng = c(1,d)
daytime = 1:d

#definir base funcional para representar los componentes principales
knots = seq(1,d,2)
norder = 4
nbasis = length(knots) + norder - 2
pcbasis = create.bspline.basis(dayrng,nbasis,norder,knots)
curv.Lfd = int2Lfd(2)

#escoger cantidad de componentes principales para la regresion
mse_pc = 2:10
for (p in 2:10){ #varia numero de CP
  gdPCA = pca.fd(gdfd,nharm=p)
  
  pca.fdPar = fdPar(pcbasis,curv.Lfd,1) #se suavizan las componentes principales funcionales 
  gdPCA = pca.fd(gdfd,nharm=p,harmfdPar=pca.fdPar)
  Xmat = gdPCA$scores
  mses = 1:10
  for (k in 1:10){ #varia fold
    test = flds[[k]]
    train = which(!(1:332 %in% test))
    lmodel = lm(corte.bloque[train]~.,data=data.frame(Xmat[train,]))
    mses[k]=mean((predict.lm(lmodel,newdata=data.frame(Xmat)[test,])-corte.bloque[test])^2)
  }
  mse_pc[p-1] = mean(mses)
}
plot(2:10,mse_pc,type='b',col="red",xlab = "componente principal",ylab="MSE por 10-folds") #grafica de calibracion de cantidad de CP

#modelo de regresion por CPF
p = which.min(mse_pc)+1
pca.fdPar = fdPar(pcbasis,curv.Lfd,1)
gdPCA = pca.fd(gdfd,nharm=p,harmfdPar=pca.fdPar)

plot(gdPCA$varprop,type='b', ylab = "% de varianza expicada", xlab = "componente principal",col="red") #graficar varianza explicada por CPF
plot(gdPCA$harmonics,lty=1) #graficar CPF

Xmat = gdPCA$score
lmodel = lm(corte.bloque~Xmat)
summary(lmodel)
Mcoefs = lmodel$coef
PCs = gdPCA$harmonics
beta = 0
for(n in 1:p){
  beta = beta + Mcoefs[n+1]*PCs[n] 
}
plot(beta,ylab = "values",main = paste("Beta estimada para",d,"dias por CP"),col="forestgreen") #reconstruccion de la funcion de regresion estimada beta

mse_PCR[which(dias == d)] = min(mse_pc)
#===============================================================================================
#REGRESION LINEAL MULTIVARIADA

d = 15
gd = datos_gd[,1:d]


mses = 1:10
for (l in 1:10){
  test = flds[[l]]
  train = which(!(1:332 %in% test))
  modeli = lm(corte.bloque[train] ~ ., data = data.frame(gd[train,]))
  mses[l] = mean((predict(modeli,data.frame(gd)[test,])-corte.bloque[test])^2)
}
mean(mses) #MSE por 10-folds
mse_MVR[which(dias==d)]=mean(mses) #MSE por 10-folds

model = lm(corte.bloque ~ ., data = data.frame(gd))
summary(model)$r.squared #R^2 del model de regresion multivariada

mse_MVR[which(dias == d)] = mean(mses)
#===============================================================================================
#COMPARACION REGRESIONES

library(ggplot2)

df1 = data.frame(group = "FR", mse_FR)
df2 = data.frame(group = "PCR", mse_PCR)
df3 = data.frame(group = "MVR", mse_MVR)
names(df1) = c("group","MSE")
names(df2) = c("group","MSE")
names(df3) = c("group","MSE")
df = rbind(df1,df3,df2)
df = cbind(c(dias,dias,dias),df)
names(df)[1]="Dias"
df$Dias = as.factor(df$Dias)
ggplot(data=df, aes(x=Dias, y=MSE, fill=group)) + geom_bar(stat="identity", position=position_dodge())+ scale_fill_brewer(palette="Paired")+theme_minimal()
