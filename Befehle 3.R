# Bayes Statistik
function(){
  xachse <- seq(3,6,0.001)
  plot(xachse, dnorm(xachse, mean=mu0, sd=sigma0), type="l", main="A-Priori-Verteilung", xlab="x", ylab="Dichte")
  
  mu1 <- ((n/sigma^2)*xquer + (1/sigma0^2)*mu0)/(n/sigma^2+1/sigma0^2)
  sigma1 <- sqrt(1/(n/sigma^2+1/sigma0^2))
  
  plot(xachse,dnorm(xachse, mean=mu1, sd=sigma1), type="l", main="A-Posteriori-Verteilung", xlab="x", ylab="Dichte")
  
  #---
  # Likelihood
  likelihood<-function(mu) prod(dnorm(x1,mean=mu,sd=sigma))
  L<-sapply(xachse,likelihood)
  
  plot(xachse,L, type="l", main="Likelihood", xlab="x", ylab="Likelihood")

}

# Log-likelihood-Methode
function(){

  l <- function(xs, mu, sigma) {
    sum(log(dnorm(x=xs,mean=mu,sd=sigma)))
  }
  l(xs,0,1)
  
  l.max = -100000000
  mu.min = NA;
  sig.min = NA;
  for (mu in seq(0.5,0.6,0.001)) {
    for (sigma in seq(0.6,0.75,0.001)) {
      lik <- l(xs, mu, sigma)
      if (lik > l.max){
        l.max <- lik
        mu.min <- mu;
        sig.min <- sigma;
      }
    }
  }
  l.max;mu.min;sig.min

}

# Bootstrap Methode
function(){

# jedesmal zu beginn!
library(boot)

# Funktion f??r den Sch??tzer
#"eine Bootstrap Methode f??r ...(die Funktion die impementiert werden muss)"
# ind => durchlauf der Daten // data => die Daten selber
schaetzer <- function(data,ind) 1/mean(data[ind])

# Die Bootstrap Methode
boot(Daten, Funktion mit 2 Variablen!, Widerholungen)
boot.erg<-boot(data=x, statistic=schaetzer, R=500)

str(boot.erg)
hist(boot.erg$t) # hist(bootbefehl$t)  => $t ist wichtig!
var(boot.erg$t)  #gesch??tzte Varianz
mean(boot.erg$t)-boot.erg$t0  #gesch??tzter Bias

# Bootstrap-Konfidenzintervall
boot.ci(boot.erg,conf=0.95,type=c("perc", "bca"))   #bca-Intervall ist besser

}

# Gleichungen aus Skript
function(){
Yi = alpha + beta*xi + Ei   (mit Ei unab ~ N(0,sigma^2))

#---
# von Hand
x = unique2010$ATM
y = unique2010$Pax

b1 = sum( (x - mean(x)) * (y - mean(y)) ) / sum((x - mean(x))^2)
b0 = mean(y) - b1 * mean(x)

# Residuen
y_hat = b0 - b1*x
y - y_hat

#R^2
1 - sum((y - yhat)^2) / sum((y - mean(y))^2)
 
# Standartabweichung der Residuen
sig.eps <- sqrt(sum((y - yhat)^2)/(n-2))

# Standarbweichung der Sch??tzer
# b1
sig.beta1 <- sig.eps * sqrt(1/sum((x-mean(x))^2))
# b0
sig.beta0 <- sig.eps * sqrt(1/n + mean(x)^2/sum((x-mean(x))^2))

# t-Werte
t0 <- b0 / sig.beta0
t1 <- b1 / sig.beta1

#p-Werte
t <- c(t0, t1)
2*(1-pt(abs(t), df=n-2))

}

# Residuen Check und Konfidenz Intervall
function(){
par(mfrow=c(2,3))
plot(lm.erg, which=1:6)

# confidenz Intervalle usw.
x<-c(42.5, 47.5, 52.5, 57.5, 62.5, 67.5)
y<-c(1.00, 0.94, 1.12, 1.11, 1.11, 1.08)

plot(x,y,ylim=c(0.5,1.5))
abline(lm.erg)

x2<-data.frame(x=seq(40,70,0.5))
conf<-predict(lm.erg,x2,interval="confidence")
lines(x2$x,conf[,2])
lines(x2$x,conf[,3])
pred<-predict(lm.erg,x2,interval="prediction")
lines(x2$x,pred[,2],lty=2)
lines(x2$x,pred[,3],lty=2)

}