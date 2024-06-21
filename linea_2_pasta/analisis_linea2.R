#Guardar en el objeto linea2, la informacion de los pesos de paquetes de pastas
linea2 <- read.table("linea2.csv", sep = ",", header = TRUE)
linea2
names(linea2)
vlinea2 = c(linea2$x1, linea2$x2, linea2$x3, linea2$x4, linea2$x5)
vlinea2
n_linea2= length(vlinea2)
n_linea2
#Estadísticos de tendencia central
media_linea2 = mean(vlinea2)
mediana_linea2 = median(vlinea2)
library(modeest)
moda_linea2 = mlv(vlinea2)
#Estadísticos de dispersión
desviacion_est_linea2 = sd(vlinea2)
varianza_linea2 = var(vlinea2)
max_linea2 = max(vlinea2)
min_linea2 = min(vlinea2)
rango_linea2 = max_linea2 - min_linea2
#Generar un data frame con la información de los estadísticos
tabla_estadisticos_linea2= data.frame(n_linea2, 
                                      media_linea2, 
                                      mediana_linea2, 
                                      moda_linea2, 
                                      desviacion_est_linea2, 
                                      varianza_linea2, 
                                      rango_linea2, 
                                      max_linea2, 
                                      min_linea2)
write.csv(tabla_estadisticos_linea2, file = "tabla_estadisticos_linea2.csv")
tabla_estadisticos_linea2
hist(vlinea2, col="darkblue", breaks="Sturges", main= "Pesos de paquetes línea 1", 
     xlab= "Pesos (g)", ylab="Frecuencia")
hist(vlinea2, col= c("lightsteelblue", "lightsteelblue1", "lightsteelblue2", 
                     "lightsteelblue3", "lightsteelblue4"), 
     breaks="Sturges", main= "Pesos de paquetes línea 1", 
     xlab= "Pesos (g)", ylab="Frecuencia")

#Prueba de hipótesis
alpha= 0.05
mu0= 200
zTest = (media_linea2 - mu0)/ (desviacion_est_linea2/sqrt(n_linea2)) 
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
qqnorm(vlinea2) 
qqline (vlinea2, col = "red")
#Shapiro.test se usa para contrastar si un conjunto de datos 
#sigue una distribución normal o no
shapiro.test(vlinea2)
#H0: los datos provienen de una distribución normal, 
#Ha: los datos no provienen de una distribución normal
#sí p-value ≤ 0.05 se rechaza la hipótesis nula 
#si p-value ≥ 0.05 no se rechaza la hipótesis nula.
#Estimación de μx y σx
#n=5 por lo que 𝑑2=2.326 (Revisar apéndice)
row.sums <- apply(linea2, 1, mean)
mini <- apply(linea2, 1 , min)
maxi <- apply(linea2, 1 , max)
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
Resumen_indices_linea2 = data.frame(Cp,Cpi,Cps,Cpk,zi,Pzi,PPMFEi,zs,Pzs,PPMFEs,zbench,Pzbench,PPMFE,Cpm)
Resumen_indices_linea2
write.csv(Resumen_indices_linea2, file = "Resumen_indices_linea2.csv")
#Representación gráfica de los resultados, guardamos el histograma en un objeto 
#llamado h, posteriormente le colocaremos las especificaciones y el ajuste de la 
#línea normal
h = hist(vlinea2, breaks= "sturges", freq = TRUE, right = TRUE, 
         col = "lightblue", border = "black", 
         main = "Peso de paquetes de pasta linea 2", 
         xlim = c (170,230), xlab = "Peso", 
         ylab = "Frecuencia", las=1,
         ylim = c(0, 35))
abline(v = ES, col = "red",lty=3) 
text(EI, 26, "EI", cex= NULL, pos=2, col="black") 
abline(v = EI, col = "red",lty=3)
text(ES, 26, "ES", cex= NULL, pos=4, col="black") 
abline(v = N , col = "red",lty=1) 
text(N, 30, "Objetivo", cex= NULL, pos=2, col="black")
myx = seq(min(vlinea2), 230, length.out= 50) 
mymean = mean(vlinea2) 
mymean_2 = mean(vlinea2)
mysd = sd(vlinea2) 
mysd_2 = sd(vlinea2)
normal = dnorm(x = myx, mean = mymean, sd = mysd) 
multiplier = h$counts / h$density 
lines (myx, normal * multiplier[1], col = "blue", lwd = 2)

#Información faltante
boxplot(vlinea2)
boxplot(vlinea2, col="lightblue", ylab= "gramos", main="Gráfico de caja del peso de paquetes")
#Límites reales
LRS= (media_linea2 + (3*desviacion_est_linea2))
LRS
LRI = (media_linea2 - (3*desviacion_est_linea2))
LRI
zTest
#Analisis de capacidad con la librería SixSigma
library(SixSigma)
ss.study.ca(vlinea2, LSL= 192, USL= 208,
            Target = 200, alpha=0.5,
            f.main= "Analisis de Capacidad línea 2",
            f.colours= c("mediumspringgreen", "mediumturquoise", "mediumvioletred","midnightblue", "mediumseagreen", "mediumslateblue"))

h_2 = hist(vlinea2, breaks= "sturges", freq = TRUE, right = TRUE, col = "lightblue", border = "black", main = "Peso de paquetes de pasta linea 2", xlim = c (170,230), ylim= c(0, 35), xlab = "Peso", ylab = "Frecuencia", las=1)
abline(v = ES, col = "red",lty=3)
text(EI, 26, "EI", cex= NULL, pos=2, col="black")
abline(v = EI, col = "red",lty=3)
text(ES, 26, "ES", cex= NULL, pos=4, col="black")
abline(v = N , col = "red",lty=3)
text(N, 30, "Objetivo", cex= NULL, pos=2, col="black")
myx = seq(170, 230, length.out= 50)
mymean = mean(vlinea2)
mysd = sd(vlinea2)
normal = dnorm(x = myx, mean = mymean, sd = mysd)
multiplier = h$counts / h$density
lines (myx, normal * multiplier[1], col = "blue", lwd = 2)
EI_2 = 192
ES_2 = 208
N_2 = 200



h_2 = hist(vlinea2, breaks= "sturges", freq = TRUE, right = TRUE, col = "lightblue", border = "black", main = "Peso de paquetes de pasta", xlim = c (170,230), ylim= c(0, 35), xlab = "Peso", ylab = "Frecuencia", las=1)
abline(v = ES_2, col = "red",lty=3)
text(EI_2, 26, "EI_2", cex= NULL, pos=2, col="black")
abline(v = EI_2, col = "red",lty=3)
text(ES_2, 26, "ES_2", cex= NULL, pos=4, col="black")
abline(v = N_2 , col = "red",lty=3)
text(N_2, 30, "Objetivo", cex= NULL, pos=2, col="black")
myx = seq(170, 230, length.out= 50)
mymean_2 = mean(vlinea2)
mysd_2 = sd(vlinea2)
normal_2 = dnorm(x = myx, mean = mymean_2, sd = mysd_2)
multiplier = h$counts / h$density
lines (myx, normal * multiplier[1], col = "blue", lwd = 2)
