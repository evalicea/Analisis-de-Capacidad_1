##############################################################################################################
#####********************UNIVERSIDAD POLITECNICA DE PENJAMO**********************************#################
#####************************CONTROL ESTADISTICO*********************************************#################
#####*********************INGENIERIA EN BIOTECNOLOGIA ***************************************#################
#####***************************MAYO-AGOSTO 2024*********************************************#################
##############################################################################################################
###############SCRIPT PARA EL ANALISIS DE LA CAPACIDAD DE PRODUCCION DEL CONTENIDO DE PASTAS##################
##############################################################################################################

#Desde un proyecto de R y asegurándose que en la misma carpeta donde se encuentra el proyecto esté la base de datos 
#linea1.csv, abrir el script y ejecutarlo.
#Es necesario limpiar todos los objetos en memoria, por ello se escribe la siguiente instrucción
rm(list=ls()) 
#Instalar las librerías
install.packages("SixSigma")
library(SixSigma)
install.packages("modeest")
library(modeest)
#Guardar en el objeto linea1, la informacion de los pesos de paquetes de pastas
linea1 <- read.table("linea1.csv", sep = ",", header = TRUE)
linea1
names(linea1)
vlinea1 = c(linea1$x1, linea1$x2, linea1$x3, linea1$x4, linea1$x5)
vlinea1
n_linea1= length(vlinea1)
n_linea1
#Estadísticos de tendencia central
media_linea1 = mean(vlinea1)
mediana_linea1 = median(vlinea1)
moda_linea1 = mlv(vlinea1)
#Estadísticos de dispersión
desviacion_est_linea1 = sd(vlinea1)
varianza_linea1 = var(vlinea1)
max_linea1 = max(vlinea1)
min_linea1 = min(vlinea1)
rango_linea1 = max_linea1 - min_linea1
#Generar un data frame con la información de los estadísticos
tabla_estadisticos_linea1= data.frame(n_linea1, 
                                      media_linea1, 
                                      mediana_linea1, 
                                      moda_linea1, 
                                      desviacion_est_linea1, 
                                      varianza_linea1, 
                                      rango_linea1, 
                                      max_linea1, 
                                      min_linea1)
write.csv(tabla_estadisticos_linea1, file = "tabla_estadisticos_linea1.csv")
tabla_estadisticos_linea1
hist(vlinea1, col="darkblue", breaks="Sturges", main= "Pesos de paquetes línea 1", 
     xlab= "Pesos (g)", ylab="Frecuencia")
hist(vlinea1, col= c("lightsteelblue", "lightsteelblue1", "lightsteelblue2", 
                     "lightsteelblue3", "lightsteelblue4"), 
     breaks="Sturges", main= "Pesos de paquetes línea 1", 
     xlab= "Pesos (g)", ylab="Frecuencia")
1
#Prueba de hipótesis
alpha= 0.05
mu0= 200
zTest = (media_linea1 - mu0)/ (desviacion_est_linea1/sqrt(n_linea1)) 
zTest
zalphaiz <- qnorm(alpha/2,lower.tail = TRUE)
zalphaiz
zalphader <- qnorm(alpha/2, lower.tail = FALSE)
zalphader
x = pretty(-3:3,n = 50) 
y = dnorm(x)
plot(x,y, ylab = "Densidad", type = "l", col = "black", 
     main="Prueba de hipotesis bilateral")
xi <- c(-3,seq(-3,zalphaiz,by = 0.001),zalphaiz)
yi <- c(0,dnorm(seq(-3,zalphaiz,by = 0.001)),0)
polygon(xi,yi,col = "red",lty = 1,lwd=1)
cor.x = c(zalphaiz, seq(zalphaiz,zalphader, by=0.001),zalphader)
cor.y = c(0,dnorm(seq(zalphaiz, zalphader, by=0.001)),0)
polygon(cor.x,cor.y, col = "lightgrey",lty = 1,lwd=2)
xd = c(zalphader,seq(zalphader,3, by=0.001),3)
yd = c(0,dnorm(seq(zalphader,3,by=0.001)),0)
polygon(xd,yd, angle = 45,col = "red",lty = 1,lwd=1)
text(zTest, 0.02, "ZTest", cex= NULL, pos=4, col="black")
abline(v = zTest, col = "red",lty=3)
text(0, 0.2 , "Zona de aceptacion H0", cex= NULL, pos=1, col="black")
#Estudio de capacidad de proceso con índices
#Determinar si los datos siguien una distribución normal
qqnorm(vlinea1) 
qqline (vlinea1, col = "red")
#Shapiro.test se usa para contrastar si un conjunto de datos 
#sigue una distribución normal o no
shapiro.test(vlinea1)
#H0: los datos provienen de una distribución normal, 
#Ha: los datos no provienen de una distribución normal
#sí p-value ≤ 0.05 se rechaza la hipótesis nula 
#si p-value ≥ 0.05 no se rechaza la hipótesis nula.
#Estimación de μx y σx
#n=5 por lo que 𝑑2=2.326 (Revisar apéndice)
row.sums <- apply(linea1, 1, mean)
mini <- apply(linea1, 1 , min)
maxi <- apply(linea1, 1 , max)
ES <- 208
EI <- 192
mu <- mean(row.sums)
N <- 200
rango = maxi - mini
d2 = 2.326
rdesvest = mean(rango)/d2
Cp = (ES-EI)/(6*rdesvest) 
Cr = (6*rdesvest)/(ES-EI) 
Cpi = (mu-EI)/(3*rdesvest) 
Cps = (ES-mu)/(3*rdesvest) 
Cpk = min(Cpi,Cps) 
#Índice de Taguchi 
tau = sqrt((rdesvest^2)+((mu-N)^2)) 
Cpm = (ES-EI)/(6*tau)
#Cálculo de índice Z
zi = 3*Cpi 
Pzi = 1-pnorm(zi) 
PPMFEi = Pzi*1000000 
zs = 3*Cps 
Pzs = 1-pnorm(zs) 
PPMFEs = Pzs*1000000
zbench = qnorm(1-(Pzi+Pzs))
Pzbench = 1-pnorm(zbench) 
PPMFE = Pzbench*1000000
#Para generar un resumen de todos los índices de capacidad se utiliza la función data.frame
Resumen_indices_linea1 = data.frame(Cp,Cpi,Cps,Cpk,zi,Pzi,PPMFEi,zs,Pzs,PPMFEs,zbench,Pzbench,PPMFE,Cpm)
Resumen_indices_linea1
write.csv(Resumen_indices_linea1, file = "Resumen_indices_linea1.csv")
#Representación gráfica de los resultados, guardamos el histograma en un objeto 
#llamado h, posteriormente le colocaremos las especificaciones y el ajuste de la 
#línea normal
h = hist(vlinea1, breaks= "sturges", freq = TRUE, right = TRUE, 
         col = "lightblue", border = "black", 
         main = "Peso de paquetes de pasta", 
         xlim = c (160,230), xlab = "Peso", 
         ylab = "Frecuencia", las=1)
abline(v = ES, col = "red",lty=3) 
text(EI, 25, "EI", cex= NULL, pos=2, col="black") 
abline(v = EI, col = "red",lty=3)
text(ES, 25, "ES", cex= NULL, pos=4, col="black") 
abline(v = N , col = "red",lty=1) 
text(N, 30, "Objetivo", cex= NULL, pos=2, col="black")
myx = seq(min(vlinea1), 230, length.out= 50) 
mymean = mean(vlinea1) 
mysd = sd(vlinea1) 
normal = dnorm(x = myx, mean = mymean, sd = mysd) 
multiplier = h$counts / h$density 
lines (myx, normal * multiplier[1], col = "blue", lwd = 2)

#Información faltante
boxplot(vlinea1)
boxplot(vlinea1, col="lightblue", ylab= "gramos", main="Gráfico de caja del peso de paquetes")
#Límites reales
LRS= (media_linea1 + (3*desviacion_est_linea1))
LRS
LRI = (media_linea1 - (3*desviacion_est_linea1))
LRI
zTest

boxplot(vlinea1)

#Analisis de capacidad con la librería SixSigma
library(SixSigma)
ss.study.ca(vlinea1, LSL= 192, USL= 208,
            Target = 200, alpha=0.5,
            f.main= "Analisis de Capacidad línea 1",
            f.colours= c("#1b859d", "#356e84", "#50576b", "#6a3f52", "#842839"),"midnightblue")


