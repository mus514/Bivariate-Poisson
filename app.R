## Diapo 2 Poison Biv

## Loi Poisson Bivariee
poisBiv <- function(x1, x2, lam1, lam2, alpha)
{
  lmin <- min(x1, x2)
  sum(dpois(0:lmin, alpha) * exp(-(lam1-alpha))*(lam1-alpha)^(x1 -(0:lmin))/factorial(x1 - (0:lmin))
      * exp(-(lam2-alpha))*(lam2-alpha)^(x2 -(0:lmin))/factorial(x2 - (0:lmin)))
}

m <- 0:100
lambda1 <- 10
lambda2 <- 5
a <- 1

# Consctuction des valeurs de F_M
fxxp <- sapply(m, function(x) sapply(m, function(y) poisBiv(x, y, lambda1, lambda2, a)))
# Construction de fs
fssp <- fss <-sapply(m+1, function(x) sum(sapply(1:x, function(y) fxxp[y, x - y + 1])))

sum(fxxp)
sum(fssp)
vxx <- 0:(length(m) - 1)
pix <- round(sapply(c(0, 5, 10, 15, seq(20, 50, 10)) ,function(x) sum(fssp*pmax(vxx - x, 0))), 5)
corre <- (sum(sapply(m, function(x) sapply(m, function(y)   sum(x*y*fxxp[x+1, y+1])))) - 50) /sqrt(50)


## FFT pour les valeurs de fs
nfft <- 2^20

fc <- c(0, 1,rep(0, nfft-2))
ffc <- fft(fc)
fs <- Re(fft(exp((lambda1-a)*(ffc-1))*exp((lambda2-a)*(ffc -1))*exp(a*(ffc^2 - 1)), inverse = T)/nfft)
sum(fs)
vs <- 0:(nfft-1)
pix1 <- round(sapply(c(0, 5, 10, 15, seq(20, 50, 10)) ,function(x) sum(fs*pmax(vs - x, 0))), 5)
ex <- sum(vs*fs)


fm1 <- dpois(m,lambda1)
fm2 <- dpois(m,lambda2)
c(sum(fm1), sum(fm2))
Fm1 <- cumsum(fm1)
Fm2 <- cumsum(fm2)

## W(x,x)


Fxx <- t(sapply(m+1, function(x) sapply(m+1, function(y) pmax(Fm1[x] + Fm2[y] - 1, 0))))
fxx <- matrix(0, nrow = length(m), ncol = length(m))
fxx[1, 1] <- Fxx[1,1]
fxx[2:101, 1] <- Fxx[2:101, 1] - Fxx[1:100, 1] 
fxx[1, 2:101] <- Fxx[1, 2:101] - Fxx[1,1:100] 
fxx[2:101, 2:101] <- Fxx[2:101, 2:101] - Fxx[1:100, 2:101]- Fxx[2:101, 1:100] + Fxx[1:100, 1:100]
sum(fxx)
corre1 <- (sum(sapply(m, function(x) sapply(m, function(y)   sum(x*y*fxx[x+1, y+1])))) - lambda2*lambda1) /sqrt(lambda2*lambda1)


##M(x,x )

Fxxm <- t(sapply(m+1, function(x) sapply(m+1, function(y) min(Fm1[x], Fm2[y]))))
fxxm <- matrix(0, nrow = length(m), ncol = length(m))
fxxm[1, 1] <- Fxxm[1,1]
fxxm[2:101, 1] <- Fxxm[2:101, 1] - Fxxm[1:100, 1] 
fxxm[1, 2:101] <- Fxxm[1, 2:101] - Fxxm[1,1:100] 
fxxm[2:101, 2:101] <- Fxxm[2:101, 2:101] - Fxxm[1:100, 2:101]- Fxxm[2:101, 1:100] + Fxxm[1:100, 1:100]
sum(fxxm)
corre2 <- (sum(sapply(m, function(x) sapply(m, function(y)   sum(x*y*fxxm[x+1, y+1])))) - lambda2*lambda1) /sqrt(lambda2*lambda1)


##Fs W

fss <-sapply(m+1, function(x) sum(sapply(1:x, function(y) fxx[y, x - y + 1])))
sum(fss)
vss <- 0:(length(fss) -1)
pis <- round(sapply(c(0, 5, 10, 15, seq(20, 50, 10)),function(x) sum(fss*pmax(vss - x, 0))), 5)
es <- sum(fss*vss)

##Fs M
fssp <-sapply(m+1, function(x) sum(sapply(1:x, function(y) fxxm[y, x - y + 1])))
sum(fssp)
piss <- round(sapply(c(0, 5, 10, 15, seq(20, 50, 10)),function(x) sum(fssp*pmax(vss - x, 0))), 5)
ess <- sum(fssp*vss)

# Compraison
cbind(pis, pix, pix1, piss)
round(cbind(corre, corre1, corre2), 4)
cbind(ex, es, ess, lambda2+lambda1)

