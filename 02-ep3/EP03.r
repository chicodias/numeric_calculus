## LABORATORIO DE SIMULAÇÃO E COMPUTAÇÃO
## Francisco Rosa Dias de Miranda
## NUSP: 4402962 RG: 50688552
## Curso: Bac. em Estatística
## 1º Semestre de 2017

  f <- quote(exp((-0.2506*(1+t)))*0.5*(1+cos(0.6440*t))) #expressao que queremos saber a integral
  ft <- function(t){eval(f)}  

  n <- 10000 #numero de pontos
  n0 <- 1000 #numero de pontos para burn-in

  x <- 0 ##valor inicial de x

##burn-in com n0 iteradas para garantir a convergencia da cadeia
  for (i in 1:n0){
  
    u <- qnorm((0.5+0.4*runif(1)),0,12)  ##valores com distr. normal entre 0 e 15
  
    alfa <- min (1,(pnorm(x)*ft(u))/(ft(x)*pnorm(u))) ##Probabilidade de aceitação
  
    if (runif(1) <= alfa)
      x <- u
  }

  ##criando um vetor y[] com o numero de amostras n
  y <- numeric()  
  
  while(length(y)< n){
  
    u <- qnorm((0.5+0.4*runif(1)),0,12) ##gero valores na distribuicao normal entre 0 e 15
    alfa <- min (1,(pnorm(x)*ft(u))/(ft(x)*pnorm(u))) ##probabilidade de aceitacao
    
    if (runif(1) <= alfa){
      x <- u
      y <- c(x,y)
    }
  }

  I <-  (1/n)*sum(ft(y)) ##vlr da integral
  
