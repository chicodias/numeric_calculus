## LABORATORIO DE SIMULAÇÃO E COMPUTAÇÃO
## Francisco Rosa Dias de Miranda
## NUSP: 4402962 RG: 50688552
## Curso: Bac. em Estatística
## 1º Semestre de 2017

library(nloptr)

## Entrada: Dados  
t1 <- c(0.01, 0.19, 0.51, 0.57, 0.70, 0.73, 0.75, 0.75, 1.11, 1.16, 1.21, 1.22, 1.24, 1.48, 1.54, 1.59, 1.61, 1.61, 1.62, 1.62, 1.71, 1.75, 1.77, 1.79, 1.88, 1.90, 1.93, 2.01, 2.16, 2.18, 2.30, 2.30, 2.41, 2.44, 2.57, 2.61, 2.62, 2.72, 2.76, 2.84, 2.96, 2.98, 3.19, 3.25, 3.31) ##45 q queimaram
t2 <- c(1.19, 3.50, 3.50, 3.50,3.50) ##5 caras q nao queimaram


min_func <- function(alpha, beta, gamma,t1,t2){
  
  s1 <- sum( (log(beta) + (beta - 1)*log(t1 + alpha)) -( beta*log(gamma)) - (((t1 + alpha)/gamma)**beta) + (alpha/gamma)**beta )
  
  s2 <- sum((alpha/gamma)**beta - ((t2 + alpha)/gamma)**beta)
  
  return(-1*(s1 + s2))
}

# weibull a posteriori

weibull <- function(alpha,beta,gamma, t1, t2){
  
  p1 <- prod(((beta*((t1 + alpha)**(beta - 1)))/(gamma**beta))* exp((-(((t1 + alpha)/gamma)**beta)))/exp(-((alpha/gamma)**beta)))
  
  p2 <- prod(exp((-(((t2 + alpha)/gamma)**beta)))/exp(-((alpha/gamma)**beta)))
  
  return (p1*p2)
}

fn <- function(x) min_func(x[0],x[1],x[2],t1, t2)
x0 <- c(1.0, 1.0, 1.0)

#otimizacao da funcao
res <- slsqp(x0,fn)

mi = res$par[2]*gamma(1 + 1.0/res$par[1])
p_est = res$par[0]/mi

# valor maximo de f* = f(theta*)
maximo <- 1 - weibull(res$par[0],res$par[1],res$par[2], t1, t2)
