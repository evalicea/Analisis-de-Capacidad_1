#Guardar en el objeto linea3, la informacion de los pesos de paquetes de pastas
linea3 <- read.table("linea3.csv", sep = ",", header = TRUE)
linea3
names(linea3)
vlinea3 = c(linea3$x1, linea3$x2, linea3$x3, linea3$x4, linea3$x5)
vlinea3
n_linea3= length(vlinea3)
n_linea3
#Estadísticos de tendencia central
media_linea3 = mean(vlinea3)
mediana_linea3 = median(vlinea3)
moda_linea3 = mlv(vlinea3)
#Estadísticos de dispersión
desviacion_est_linea3 = sd(vlinea3)
varianza_linea3 = var(vlinea3)
max_linea3 = max(vlinea3)
min_linea3 = min(vlinea3)
rango_linea3 = max_linea3 - min_linea3
#Generar un data frame con la información de los estadísticos
tabla_estadisticos_linea3= data.frame(n_linea3, 
                                      media_linea3, 
                                      mediana_linea3, 
                                      moda_linea3, 
                                      desviacion_est_linea3, 
                                      varianza_linea3, 
                                      rango_linea3, 
                                      max_linea3, 
                                      min_linea3)
write.csv(tabla_estadisticos_linea3, file = "tabla_estadisticos_linea3.csv")
tabla_estadisticos_linea3
hist(vlinea3, col="darkblue", breaks="Sturges", main= "Pesos de paquetes línea 1", 
     xlab= "Pesos (g)", ylab="Frecuencia")
hist(vlinea3, col= c("lightsteelblue", "lightsteelblue1", "lightsteelblue2", 
                     "lightsteelblue3", "lightsteelblue4"), 
     breaks="Sturges", main= "Pesos de paquetes línea 1", 
     xlab= "Pesos (g)", ylab="Frecuencia")

#Prueba de hipótesis
alpha= 0.05
mu0= 200
zTest = (media_linea3 - mu0)/ (desviacion_est_linea3/sqrt(n_linea3)) 
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
polygon(cor.x,cor.y, col = "grey",lty = 1,lwd=2)
xd = c(zalphader,seq(zalphader,3, by=0.001),3)
yd = c(0,dnorm(seq(zalphader,3,by=0.001)),0)
polygon(xd,yd, angle = 45,col = "red",lty = 1,lwd=1)
text(zTest, 0.02, "ZTest", cex= NULL, pos=4, col="black")
abline(v = zTest, col = "black",lty=3)
text(0, 0.2 , "Zona de aceptacion H0", cex= NULL, pos=1, col="black")
#Estudio de capacidad de proceso con índices
#Determinar si los datos siguien una distribución normal
qqnorm(vlinea3) 
qqline (vlinea3, col = "red")
#Shapiro.test se usa para contrastar si un conjunto de datos 
#sigue una distribución normal o no
shapiro.test(vlinea3)
#H0: los datos provienen de una distribución normal, 
#Ha: los datos no provienen de una distribución normal
#sí p-value ≤ 0.05 se rechaza la hipótesis nula 
#si p-value ≥ 0.05 no se rechaza la hipótesis nula.
#Estimación de μx y σx
#n=5 por lo que 𝑑2=2.326 (Revisar apéndice)
row.sums <- apply(linea3, 1, mean)
mini <- apply(linea3, 1 , min)
maxi <- apply(linea3, 1 , max)
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
Resumen_indices_linea3 = data.frame(Cp,Cpi,Cps,Cpk,zi,Pzi,PPMFEi,zs,Pzs,PPMFEs,zbench,Pzbench,PPMFE,Cpm)
Resumen_indices_linea3
write.csv(Resumen_indices_linea3, file = "Resumen_indices_linea3.csv")
#Representación gráfica de los resultados, guardamos el histograma en un objeto 
#llamado h, posteriormente le colocaremos las especificaciones y el ajuste de la 
#línea normal
h = hist(vlinea3, breaks= "sturges", freq = TRUE, right = TRUE, 
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
myx = seq(min(vlinea3), 230, length.out= 50) 
mymean = mean(vlinea3) 
mysd = sd(vlinea3) 
normal = dnorm(x = myx, mean = mymean, sd = mysd) 
multiplier = h$counts / h$density 
lines (myx, normal * multiplier[1], col = "blue", lwd = 2)

#Información faltante
boxplot(vlinea3)
boxplot(vlinea3, col="lightblue", ylab= "gramos", main="Gráfico de caja del peso de paquetes")
#Límites reales
LRS= (media_linea3 + (3*desviacion_est_linea3))
LRS
LRI = (media_linea3 - (3*desviacion_est_linea3))
LRI
zTest

#Modificación histograma

h = hist(vlinea3, breaks= "sturges", freq = TRUE, right = TRUE, col = "darkslategray", 
         border = "black", main = "Peso de paquetes de pasta en linea 3", 
         xlim = c (160,230), xlab = "Peso", ylab = "Frecuencia", 
         las=1, ylim= c(0, 35))
abline(v = ES, col = "firebrick4",lty=3) 
text(EI, 26, "EI", cex= NULL, pos=2, col="black") 
abline(v = EI, col = "red",lty=3)
text(ES, 26, "ES", cex= NULL, pos=4, col="black") 
abline(v = N , col = "red",lty=3) 
text(N, 30, "Objetivo", cex= NULL, pos=2, col="black")
myx = seq(160, 230, length.out= 50) 
mymean = mean(vlinea3) 
mysd = sd(vlinea3) 
normal = dnorm(x = myx, mean = mymean, sd = mysd) 
multiplier = h$counts / h$density 
lines (myx, normal * multiplier[1], col = "red", lwd = 2)

#Analisis de capacidad con la librería SixSigma
library(SixSigma)
ss.study.ca(vlinea3, LSL= 192, USL= 208,
            Target = 200, alpha=0.5,
            f.main= "Analisis de Capacidad línea 3",
            f.colours= c("midnightblue", "mediumvioletred", "mediumseagreen", "mediumslateblue","mediumspringgreen", "mediumturquoise"))


