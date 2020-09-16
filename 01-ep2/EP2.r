## LABORATORIO DE SIMULAÇÃO E COMPUTAÇÃO
## Francisco Rosa Dias de Miranda
## NUSP: 4402962 RG: 50688552
## Curso: Bac. em Estatística
## 1º Semestre de 2017

  prec <- 0.0009 ##precisão da integral
  f <- quote(exp((-0.2506*(1+t))*0.5*(1+cos(0.6440*t)))) #expressao que queremos saber a integral
  ft <- function(t){eval(f)}  
  a <- 0 ##intervalo de integração
  b <- 1

#Monte Carlo Simples
  nMCS <- 10^3 ##Numero de ptos inicial
  SigMCS <- 1 ##Erro = 1 para primeira iteração

  while(SigMCS>prec){
    nMCS <-nMCS+10^3 ##Aumentando a amostragem em 1000 pontos a cada iteração
  
    U <- ft(runif(nMCS,a,b)) ##Calculando f nos n pontos de uma distribuição uniforme em [a,b]
  
    SigMCS <<- sqrt((sum(U^2) / nMCS) - var(U)) / sqrt(nMCS) ##Erro de Monte Carlo Simples
  
    Imcs <<- mean(U) ##Valor da Integral 
    }

#Hit or Miss
 
  h <<- 1      ##Altura do retangulo  
  nHOM <- 10^3 ##Numero de ptos inicial
  SigHOM <- 1 ##Erro = 1 para primeira iteração 
    
  while(SigHOM>prec){
    nHOM <- nHOM+10^3
  
    Ut = runif(nHOM,a,b)   ##Uniformes para (Ti,Yi)
    Uy = runif(nHOM,a,h)
  
    ns <- Uy <= ft(Ut) ##Avaliando em quais pontos há um 'hit'
  
    Ihom <<- (b-a)*h*mean(ns) ##Valor da Integral por Hit or Miss
    SigHOM <<- sd(ns)/sqrt(nHOM) ##Incerteza por Hit or Miss
  }
  

##Importance sampling
  
  nIS <- 10^3 ##Numero de ptos inicial
  SigIS <- 1 ##Erro para a primeira iteração

  while(SigIS>prec){
    nIS <- nIS+10^3

    B <- rbeta(nIS,1,1.5) ##Usando uma distribuição Beta com parâmetros 1 e 1.5

    Vm <- ft(B) / dbeta(B,1,1.5)  ##Valor da função dividido pela densidade
    
    Iis <- mean(Vm) ##Valor da Integral Importance Sampling
    SigIS <- sd(Vm)/sqrt(nIS) ##Incerteza Importance Sampling
  }
  
##Função Quadratica como Variável de Controle
  
  nCV <- 10^3 ##Numero de ptos inicial
  SigCV <- 1 ##Erro para primeira iteracao
  
  g <- quote(0.05770457*t^2 -0.19983003*t + 0.77882512) ##função quadratica utilizada
  gt <- function(t){eval(g)}
  
  while(SigCV>prec){
    nCV <- nCV+10^3
    
    U <- runif(nCV,a,b)  ##Distribuição Uniforme em [a,b]
    Mxi <- ft(U)-gt(U)
    
    Icv <- integrate(gt,a,b)$value+mean(Mxi) ##Valor da Integral por Variavel de Controle
    SigCV <- sd(Mxi)/sqrt(nCV)  #Erro amostral
  }
  
##Imprimindo os resultados
  
  R <- matrix(c(nMCS,Imcs,SigMCS,nHOM,Ihom,SigHOM,nIS,Iis,SigIS,nCV,Icv,SigCV),nrow=4,ncol=3,byrow = TRUE)
  rownames(R) <- c("Crude","Hit or Miss","Importance Sampl","Var de Controle")
  colnames(R) <- c("n:","Vlr Numerico:","Precisão:")
  print(R)
    
  