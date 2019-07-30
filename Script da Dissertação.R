#########################################################################################
###################    Mestrado em Matemática e Aplicações - UA    ######################
###################                 2018/2019                      ######################
###################                                                ######################
###################           Dados Composicionais                 ######################
###################                                                ######################
###################              Letícia Leite                     ######################
#########################################################################################

# Leitura dos dados ####
library(XLConnectJars)
library(XLConnect)
ficheiro2005 = loadWorkbook("MatrizBase2005_ler.xlsx")

#Abrir os packages que serão usados no script
library(compositions)
library(robCompositions)
library(robustbase)
library(mvoutlier)
library(gplots)
library(MASS)
library(rgdal) #Biblioteca de abstração de geoespacial
library(stringi)

###################****************************************###################
###################           Matriz IdadeI 2005           ###################
###################****************************************###################

# Leitura da matriz ####
Y1inicial = readWorksheet(ficheiro2005, sheet="IdadeI2005")
Y1 = Y1inicial
GrandTotal = Y1[-309, 7] #Número de habitantes internos que entraram, até à data de 2011, em cada município
Y1 = Y1[-309,-c(1,7)] #Retirar a variável 'GrandTotal'
colnames(Y1) = c("Idade0-14","Idade15-24","Idade25-39","Idade40-64","Idade65+")

dadosY1 = Y1/GrandTotal #Distribuição das idades dos habitantes internos que entraram em cada município
colnames(dadosY1) = colnames(Y1)

#Variável que guarda só os nomes dos municípios
Y1names = Y1inicial[-309,-7] #Adicionar a coluna dos nomes e retirar a coluna do 'GrandTotal'
colnames(Y1names) = c("Municípios","Idade0-14","Idade15-24","Idade25-39","Idade40-64","Idade65+")

#------------------------------------------------------------------------------#
# Metodologias NUMÉRICAS ####

### PBS segundo Filzmoser ####
dadosY1ilr = pivotCoord(dadosY1) #Coordenadas ilr-transformadas 
names(dadosY1ilr) = paste0("z", 1:ncol(dadosY1ilr))
Base = orthbasis(ncol(dadosY1)) #Base ortonormal
matrizContr = Base$V #Matriz de contrastes (colunas ortonormais)
PBS = Base$basisv #Tabela da PBS

#********************************************************************************#

### Estimador MCD ####
robustezY1 = covMcd(dadosY1ilr, alpha = 0.75, nsamp = "deterministic") #Algoritmo Fast-MCD
#NOTA: obtém-se as estimativas robustas da média e da matriz de covariância
summary(robustezY1)

dY1 = determinant(robustezY1$cov,logarithm = TRUE)$modulus[1] #-13.00288
ddY1 = determinant(robustezY1$cov,logarithm = FALSE)$modulus[1] #2.253831e-06
#NOTA: ln(2.253831e-06)=-13.00288

#Outliers de 'covMcd'
outCovY1 = which(robustezY1$mcd.wt==0) #34 outliers
Y1names[outCovY1,1] #Nomes dos municípios outliers

#Gráfico de 'covMcd'
plot(robustezY1, which = "distance", classic = TRUE)

#********************************************************************************#

### Distância de Mahalanobis ####
#Adjusted Quantile Plot
#Função aq.plot original alterada
aq.plot = function (x, delta = qchisq(0.975, df = ncol(x)), quan = 1/2,alpha = 0.05) 
{
  if (is.vector(x) == TRUE || ncol(x) == 1) {
    stop("x must be at least two-dimensional")
  }
  covr <- covMcd(x, alpha = quan)
  mcd <- determinant(covr$cov, logarithm = TRUE)$modulus[1] #Acrescentei esta linha
  dist <- mahalanobis(x, center = covr$center, cov = covr$cov)
  s <- sort(dist, index = TRUE)
  z <- x
  if (ncol(x) > 2) {
    p <- princomp(x, covmat = covr)
    z <- p$scores[, 1:2]
    sdprop <- (p$sd[1] + p$sd[2])/sum(p$sd)
    cat("Projection to the first and second robust principal components.\\n")
    cat("Proportion of total variation (explained variance): ")
    cat(sdprop)
    cat("\\n")
  }
  par(mfrow = c(2, 2), mai = c(0.8, 0.6, 0.2, 0.2), mgp = c(2.4,1, 0))
  plot(z, col = 3, type = "n", xlab = "", ylab = "")
  text(z, dimnames(as.data.frame(z))[[1]], col = 3, cex = 0.8)
  plot(s$x, (1:length(dist))/length(dist), col = 3, xlab = "Ordered squared robust distance", 
       ylab = "Cumulative probability", type = "n")
  text(s$x, (1:length(dist))/length(dist), as.character(s$ix),col = 3, cex = 0.8)
  t <- seq(0, max(dist), by = 0.01)
  lines(t, pchisq(t, df = ncol(x)), col = 6)
  abline(v = delta, col = 5)
  text(x = delta, y = 0.4, paste(100 * (pchisq(delta, df = ncol(x))), 
                                 "% Quantile", sep = ""), col = 5, pos = 2, srt = 90,cex = 0.8)
  xarw <- arw(x, covr$center, covr$cov, alpha = alpha)
  if (xarw$cn < Inf) {
    abline(v = xarw$cn, col = 4)
    text(x = xarw$cn, y = 0.4, "Adjusted Quantile", col = 4,pos = 4, srt = 90, cex = 0.8)
  }
  plot(z, col = 3, type = "n",
       main = paste("Outliers based on ",100 * (pchisq(delta, df = ncol(x))), "% quantile", sep = ""), 
       xlab = "", ylab = "")
  if (any(dist > delta)) {
    text(z[dist > delta, 1], z[dist > delta, 2],
         dimnames(as.data.frame(x)[dist > delta, ])[[1]], col = 2, cex = 0.8)
  }
  if (any(dist <= delta)) {
    text(z[dist <= delta, 1], z[dist <= delta, 2],
         dimnames(as.data.frame(x)[dist <= delta, ])[[1]], col = 3, cex = 0.8)
  }
  plot(z, col = 3, type = "n", main = "Outliers based on adjusted quantile", 
       xlab = "", ylab = "")
  if (xarw$cn < Inf) {
    text(z[dist > xarw$cn, 1], z[dist > xarw$cn, 2],
         dimnames(as.data.frame(x)[dist > xarw$cn, ])[[1]], col = 2, cex = 0.8)
  }
  text(z[dist <= xarw$cn, 1], z[dist <= xarw$cn, 2],
       dimnames(as.data.frame(x)[dist <= xarw$cn, ])[[1]], col = 3, cex = 0.8)
  
  o <- (sqrt(dist) > max(sqrt(xarw$cn), sqrt(qchisq(0.975, dim(x)[2]))))
  aqplot <- list(outliers = o, MCD1 = mcd, MatrizCov = covr$cov, media = covr$center)
  #Acrescentei os dois últimos parâmetros
  return(aqplot)
}

#Usar a função aqplot alterada
adjquanY1 = aq.plot(dadosY1ilr, delta = qchisq(0.975, df = ncol(dadosY1ilr)),
                    quan = 0.75, alpha = 0.975)

#Valor do estimador MCD
adjquanY1$MCD1

#Outliers
outAQplotY1 = which(adjquanY1$outliers == TRUE) #Obtém-se as observações que são outliers
length(outAQplotY1)
Y1names[outAQplotY1,1] #Nomes dos municípios outliers

#SOLUÇÃO: guardar os outliers associados ao menor valor de MCD
out1Y1Final = c(3,11,34,63,67,70,85,86,100,103,104,121,132,138,141,144,145,148,156,159,188,
                200,212,233,267,270,284,294,296,307)
# 30 outliers com MCD=-13.01813

#Nomes dos municípios outliers para a solução
Y1names[c(3,11,34,63,67,70,85,86,100,103,104,121,132,138,141,144,145,148,156,159,188,
          200,212,233,267,270,284,294,296,307),1]

#................................................................................#

#Distance-Distance Plot
#Função ddplot original alterada
dd.plot = function (x, quan = 1/2, alpha = 0.025, ...) 
{
  if (!is.matrix(x) && !is.data.frame(x)) 
    stop("x must be matrix or data.frame")
  rob <- covMcd(x, alpha = quan)
  mcd <- determinant(rob$cov, logarithm = TRUE)$modulus[1] #Acrescentei esta linha
  xarw <- arw(x, rob$center, rob$cov, alpha = alpha)
  distcla <- sqrt(mahalanobis(x, center = apply(x, 2, mean), 
                              cov = cov(x)))
  distrob <- sqrt(mahalanobis(x, center = rob$center, cov = rob$cov))
  plot(distcla, distrob, main = "Distance-Distance Plot", xlab = "Mahalanobis Distance", 
       ylab = "Robust Distance", type = "n", ...)
  if (xarw$cn != Inf) {
    alpha <- sqrt(c(xarw$cn, qchisq(c(0.75, 0.5, 0.25), ncol(x))))
  }
  else {
    alpha <- sqrt(qchisq(c(0.975, 0.75, 0.5, 0.25), ncol(x)))
  }
  abline(h = alpha[1])
  abline(v = alpha[1])
  abline(a = 0, b = 1)
  lpch <- c(3, 3, 16, 1, 1)
  lcex <- c(1.5, 1, 0.5, 1, 1.5)
  lalpha <- length(alpha)
  xs <- scale(x) - min(scale(x))
  eucl <- sqrt(apply(xs^2, 1, sum))
  rbcol <- rev(rainbow(nrow(x), start = 0,
                       end = 0.7))[as.integer(cut(eucl, nrow(x), labels = 1:nrow(x)))]
  rd <- distrob
  for (j in 1:lalpha) {
    if (j == 1) {
      points(distcla[rd >= alpha[j]], distrob[rd >= alpha[j]], 
             pch = lpch[j], cex = lcex[j], col = rbcol[rd >= alpha[j]])
    }
    if (j > 1 & j < lalpha) 
      points(distcla[rd < alpha[j - 1] & rd >= alpha[j]], 
             distrob[rd < alpha[j - 1] & rd >= alpha[j]], 
             cex = lcex[j], pch = lpch[j], col = rbcol[rd < alpha[j - 1] & rd >= alpha[j]])
    if (j == lalpha) {
      points(distcla[rd < alpha[j - 1] & rd >= alpha[j]], 
             distrob[rd < alpha[j - 1] & rd >= alpha[j]], 
             cex = lcex[j], pch = lpch[j], col = rbcol[rd < alpha[j - 1] & rd >= alpha[j]])
      points(distcla[rd < alpha[j]], distrob[rd < alpha[j]], 
             pch = lpch[j + 1], cex = lcex[j + 1], col = rbcol[rd < alpha[j]])
    }
  }
  o <- (rd > min(sqrt(xarw$cn), sqrt(qchisq(0.975, dim(x)[2]))))
  ddplot <- list(outliers = o, md.cla = distcla, md.rob = distrob, MCD2 = mcd) #Acrescentei o último parâmetro
  return(ddplot)
}

#Usar a função ddplot alterada
DDplotY1 = dd.plot(dadosY1ilr, quan = 0.75, alpha = 0.975)

#Valor do estimador MCD
DDplotY1$MCD2

#Outliers
outDDplotY1 = which(DDplotY1$outliers == TRUE)
length(outDDplotY1)
Y1names[outDDplotY1,1]

#SOLUÇÃO: guardar os outliers associados ao menor valor de MCD
out2Y1Final = c(3,11,34,63,67,70,81,85,86,100,103,104,108,121,132,138,141,144,145,148,156,
                159,188,200,204,212,233,262,267,270,284,294,296,307)
# 34 outliers com MCD=-13.01813

#Nomes dos municípios outliers para a solução
Y1names[c(3,11,34,63,67,70,81,85,86,100,103,104,108,121,132,138,141,144,145,148,156,159,188,
          200,204,212,233,262,267,270,284,294,296,307),1]

#................................................................................#

#Outlier CoDa
#Função outCoDa original alterada
outCoDa = function (x, quantile = 0.975, method = "robust", h = 1/2, coda = TRUE) 
{
  if (dim(x)[2] < 2) 
    stop("need data with at least 2 variables")
  covEst <- function(x, type) {
    standard <- function(x) {
      list(mean = colMeans(x, na.rm = TRUE), varmat = cov(x))
    }
    robust <- function(x) {
      v <- robustbase::covMcd(x)
      mcd <- base::determinant(v$cov, logarithm = TRUE)$modulus[1] #Acrescentei esta linha
      list(mean = v$center, varmat = v$cov, MCD = mcd) #Acrescentei o último parâmetro
    }
    switch(type, standard = standard(x), robust = robust(x))
  }
  if (!is.logical(coda) & !is.function(coda)) {
    stop("coda must be logical or function")
  }
  if (!is.logical(coda)) {
    x <- coda(x)
  }
  else if (coda) {
    x <- pivotCoord(x)
  }
  cv <- covEst(x, "robust")
  cvc <- covEst(x, "standard")
  dM <- sqrt(mahalanobis(x, center = cv$mean, cov = cv$varmat))
  dMc <- sqrt(mahalanobis(x, center = cvc$mean, cov = cvc$varmat))
  limit <- sqrt(qchisq(p = quantile, df = ncol(x) - 1))
  res <- list(mahalDist = dM, limit = limit, outlierIndex = dM > limit, method = method,
              om2 = dMc > limit, m2 = dMc, coda = coda, MCD3 = cv$MCD) #Acrescentei o último parâmetro
  class(res) <- "outCoDa"
  invisible(res)
}

#Usar a função outCoDa alterada
outRobY1 = outCoDa(dadosY1, quantile = 0.975, method = "robust", h = 231)
str(outRobY1)

#Valor do estimador MCD
outRobY1$MCD3

#Outliers
outROBY1 = which(outRobY1$outlierIndex == TRUE)
Y1names[outROBY1,1]

#SOLUÇÃO: guardar os outliers associados ao menor valor de MCD
out3Y1Final = c(3,11,18,20,32,34,59,63,67,70,81,85,86,100,101,102,103,104,108,121,132,138,
                141,144,145,148,153,156,159,163,188,200,204,212,222,223,233,248,262,267,270,
                277,279,284,288,294,296,307)
# 48 outliers com MCD=-13.16091

#Nomes dos municípios outliers para a solução
Y1names[c(3,11,18,20,32,34,59,63,67,70,81,85,86,100,101,102,103,104,108,121,132,138,
          141,144,145,148,153,156,159,163,188,200,204,212,222,223,233,248,262,267,270,
          277,279,284,288,294,296,307),1]

#------------------------------------------------------------------------------#
# Metodologias GRÁFICAS ####
#Função mvoutlier.CoDa original alterada
mvoutlier.CoDa = function (x, quan = 0.75, alpha = 0.025,
                           col.quantile = c(0,0.05, 0.1, 0.5, 0.9, 0.95, 1),
                           symb.pch = c(3, 3, 16, 1,1), symb.cex = c(1.5, 1, 0.5, 1, 1.5),
                           adaptive = TRUE)
{
  if (!is.matrix(x) && !is.data.frame(x)) 
    stop("x must be matrix or data.frame")
  if (ncol(x) < 3) 
    stop("x must have at least 3 compositional parts")
  Z <- pivotCoord(x)
  V <- (orthbasis(ncol(x))$V)
  Zj <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  for (j in 1:ncol(x)) {
    Zj[, j] <- pivotCoord(cbind(x[, j], x[, -j]))[, 1]
  }
  dimnames(Zj)[[2]] <- names(x)
  rob <- covMcd(Z, alpha = quan)
  mcd <- determinant(rob$cov, logarithm = TRUE)$modulus[1] #Acrescentei esta linha
  if (adaptive) {
    Zarw <- arw(Z, rob$center, rob$cov, alpha = alpha)
    if (Zarw$cn != Inf) {
      alpha1 <- sqrt(c(Zarw$cn, qchisq(c(0.75, 0.5, 0.25), ncol(Z))))
    }
    else {
      alpha1 <- sqrt(qchisq(c(0.975, 0.75, 0.5, 0.25), ncol(Z)))
    }
    rd2 <- mahalanobis(Z, center = Zarw$m, cov = Zarw$c)
  }
  else {
    cutoff <- qchisq(1 - alpha, ncol(Z))
    rd2 <- mahalanobis(Z, center = rob$center, cov = rob$cov)
    Zarw <- list(m = rob$center, c = rob$cov, cn = cutoff, w = as.numeric(rd2 < cutoff))
    alpha1 <- sqrt(qchisq(c(0.975, 0.75, 0.5, 0.25), ncol(Z)))
  }
  rd <- sqrt(rd2)
  covobj <- list(center = Zarw$m, cov = Zarw$c, n.obs = length(rd), mah = rd)
  Z.pca <- suppressWarnings(princomp(Z, covmat = covobj, cor = FALSE))
  pcaclr <- Z.pca
  eval <- eigen(Zarw$c)$values
  pcaclr$sdev <- sqrt(eval)
  pcaclr$scores <- Z.pca$scores
  pcaclr$loadings <- V %*% Z.pca$loadings
  dimnames(pcaclr$loadings)[[1]] <- names(x)
  pcaobj <- list(method = "robust", eigenvalues = eval, princompOutputClr = pcaclr)
  class(pcaobj) <- "pcaCoDa"
  Zcent <- scale(Zj, center = apply(Zj, 2, median), scale = FALSE)
  eucl <- apply(abs(Zcent), 1, median)
  out <- (!Zarw$w)
  lq <- length(col.quantile)
  colcol <- rev(rainbow(lq - 1, start = 0,
                        end = 0.7))[as.integer(cut(eucl, quantile(eucl, col.quantile, labels = 1:(lq - 1))))]
  colbw <- rev(gray(seq(from = 0.1, to = 0.9,
                        length = lq - 1)))[as.integer(cut(eucl, quantile(eucl, col.quantile, labels = 1:(lq - 1))))]
  pchvec <- rep(symb.pch[1], nrow(Zj))
  cexvec <- rep(symb.cex[1], nrow(Zj))
  if (length(symb.pch) == 5 & length(symb.cex) == 5) {
    lalpha <- length(alpha1)
    for (j in 1:(lalpha)) {
      pchvec[rd < alpha1[j]] <- symb.pch[j + 1]
      cexvec[rd < alpha1[j]] <- symb.cex[j + 1]
    }
  }
  mvoutlierCoDa <- list(ilrvariables = Zj, outliers = out, pcaobj = pcaobj, colcol = colcol,
                        colbw = colbw, pchvec = pchvec, cexvec = cexvec,
                        media = rob$center, cov = rob$cov, MCD4 = mcd)
  #Acrescentei os três últimos parâmetros
  class(mvoutlierCoDa) <- "mvoutlierCoDa"
  return(mvoutlierCoDa)
}

#Usar a função mvoutlierCoDa alterada
resY1 = mvoutlier.CoDa(dadosY1, quan = 0.75, alpha = 0.975, adaptive = TRUE)
str(resY1)

#Valor do estimador MCD
resY1$MCD4

#Outliers
outResY1 = which(resY1$outliers == TRUE)
length(outResY1)
Y1names[outResY1,1]

#Matriz de covariância
matrizCovY1 = resY1$cov
colnames(matrizCovY1) = c("z1","z2","z3","z4")
row.names(matrizCovY1) = c("z1","z2","z3","z4")
matrizCovY1

#Vetor da média aritmética
resY1$media

#................................................................................#

#SOLUÇÃO: guardar os outliers associados ao menor valor de MCD
out4Y1Final = c(3,11,34,63,67,70,85,86,100,103,104,121,132,138,141,144,145,148,156,159,188,200,
                212,233,262,267,270,284,294,296,307)
# 31 outliers com MCD=-13.01813

#Nomes dos municípios outliers para a solução
Y1names[c(3,11,34,63,67,70,85,86,100,103,104,121,132,138,141,144,145,148,156,159,188,200,
          212,233,262,267,270,284,294,296,307),1]

#Matriz de covariância associada ao menor valor de MCD
matrizCovY1Final = read.table("MatrizCovarianciaY1.txt", sep = "", dec = ".", header = TRUE)
# para MCD=-13.01813

#Média aritmética associada ao menor valor de MCD
mediaY1Final = c(-1.5537044,-0.5070277,0.4915172,0.8187603)
# para MCD=-13.01813

#Novamente a função anterior, mas alterou-se o valor de MCD e as estimativas robustas dos parâmetros \mu e \Sigma
mvoutlier.CoDa = function (x, quan = 0.75, alpha = 0.025,
                           col.quantile = c(0,0.05, 0.1, 0.5, 0.9, 0.95, 1),
                           symb.pch = c(3, 3, 16, 1,1), symb.cex = c(1.5, 1, 0.5, 1, 1.5),
                           adaptive = TRUE)
{
  if (!is.matrix(x) && !is.data.frame(x)) 
    stop("x must be matrix or data.frame")
  if (ncol(x) < 3) 
    stop("x must have at least 3 compositional parts")
  Z <- pivotCoord(x)
  V <- (orthbasis(ncol(x))$V)
  Zj <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  for (j in 1:ncol(x)) {
    Zj[, j] <- pivotCoord(cbind(x[, j], x[, -j]))[, 1]
  }
  dimnames(Zj)[[2]] <- names(x)
  matrizCovY1Final <- read.table("MatrizCovarianciaY1.txt", sep = "", dec = ".", header = TRUE) #Acrescentei
  mediaY1Final <- c(-1.5537044,-0.5070277,0.4915172,0.8187603) #Acrescentei
  mcd <- -13.01813 #Coloquei o valor fixo
  if (adaptive) {
    Zarw <- arw(Z, mediaY1Final, matrizCovY1Final, alpha = alpha)
    if (Zarw$cn != Inf) {
      alpha1 <- sqrt(c(Zarw$cn, qchisq(c(0.75, 0.5, 0.25), ncol(Z))))
    }
    else {
      alpha1 <- sqrt(qchisq(c(0.975, 0.75, 0.5, 0.25), ncol(Z)))
    }
    rd2 <- mahalanobis(Z, center = Zarw$m, cov = Zarw$c)
  }
  else {
    cutoff <- qchisq(1 - alpha, ncol(Z))
    rd2 <- mahalanobis(Z, center = mediaY2Final, cov = matrizCovY2Final)
    Zarw <- list(m = mediaY2Final, c = matrizCovY2Final, cn = cutoff, w = as.numeric(rd2 < cutoff))
    alpha1 <- sqrt(qchisq(c(0.975, 0.75, 0.5, 0.25), ncol(Z)))
  }
  rd <- sqrt(rd2)
  covobj <- list(center = Zarw$m, cov = Zarw$c, n.obs = length(rd), mah = rd)
  Z.pca <- suppressWarnings(princomp(Z, covmat = covobj, cor = FALSE))
  pcaclr <- Z.pca
  eval <- eigen(Zarw$c)$values
  pcaclr$sdev <- sqrt(eval)
  pcaclr$scores <- Z.pca$scores
  pcaclr$loadings <- V %*% Z.pca$loadings
  dimnames(pcaclr$loadings)[[1]] <- names(x)
  pcaobj <- list(method = "robust", eigenvalues = eval, princompOutputClr = pcaclr)
  class(pcaobj) <- "pcaCoDa"
  Zcent <- scale(Zj, center = apply(Zj, 2, median), scale = FALSE)
  eucl <- apply(abs(Zcent), 1, median)
  out <- (!Zarw$w)
  lq <- length(col.quantile)
  colcol <- rev(rainbow(lq - 1, start = 0,
                        end = 0.7))[as.integer(cut(eucl, quantile(eucl, col.quantile, labels = 1:(lq - 1))))]
  colbw <- rev(gray(seq(from = 0.1, to = 0.9,
                        length = lq - 1)))[as.integer(cut(eucl, quantile(eucl, col.quantile, labels = 1:(lq - 1))))]
  pchvec <- rep(symb.pch[1], nrow(Zj))
  cexvec <- rep(symb.cex[1], nrow(Zj))
  if (length(symb.pch) == 5 & length(symb.cex) == 5) {
    lalpha <- length(alpha1)
    for (j in 1:(lalpha)) {
      pchvec[rd < alpha1[j]] <- symb.pch[j + 1]
      cexvec[rd < alpha1[j]] <- symb.cex[j + 1]
    }
  }
  mvoutlierCoDa <- list(ilrvariables = Zj, outliers = out, pcaobj = pcaobj, colcol = colcol,
                        colbw = colbw, pchvec = pchvec, cexvec = cexvec,
                        cov = matrizCovY1Final, MCD4 = mcd, media = mediaY1Final)
  #Acrescentei os três últimos parâmetros
  class(mvoutlierCoDa) <- "mvoutlierCoDa"
  return(mvoutlierCoDa)
}

resY1Alt = mvoutlier.CoDa(dadosY1, quan = 0.75, alpha = 0.975, adaptive = TRUE)
str(resY1Alt)

#Confirmar os resultados
resY1Alt$cov
resY1Alt$MCD4
resY1Alt$media

#Outliers (confirmar)
outResY1Alt = which(resY1Alt$outliers == TRUE)
length(outResY1Alt)
Y1names[outResY1Alt,1]

#................................................................................#

#Biplot: para visualização dos outliers
plot(resY1Alt, which = "biplot", onlyout = TRUE, symb = TRUE, symbtxt = TRUE)
#Apenas os outliers são exibidos e identificados por números
#NOTA: os outliers estão ordenados

#Biplot: para visualização dos ângulos entre os links
biplotAngY1 = pcaCoDa(dadosY1, method = "robust")
summary(biplotAngY1)
biplot(biplotAngY1, xlabs = rep(".", nrow(dadosY1)))

#................................................................................#

#Gráfico de dispersão univariado
plot(resY1Alt, which = "uni", onlyout = TRUE, symb = TRUE, symbtxt = TRUE)
#Apenas os outliers são exibidos e identificados por números

#********************************************************************************#

#Matriz de variação (Tabela 4.3)
varY1=matrix(0,5,5)
varY1[1,2]=var(log(dadosY1[,1]/dadosY1[,2]))
varY1[1,3]=var(log(dadosY1[,1]/dadosY1[,3]))
varY1[1,4]=var(log(dadosY1[,1]/dadosY1[,4]))
varY1[1,5]=var(log(dadosY1[,1]/dadosY1[,5]))

varY1[2,3]=var(log(dadosY1[,2]/dadosY1[,3]))
varY1[2,4]=var(log(dadosY1[,2]/dadosY1[,4]))
varY1[2,5]=var(log(dadosY1[,2]/dadosY1[,5]))

varY1[3,4]=var(log(dadosY1[,3]/dadosY1[,4]))
varY1[3,5]=var(log(dadosY1[,3]/dadosY1[,5]))

varY1[4,5]=var(log(dadosY1[,4]/dadosY1[,5]))

varY1 #Imprimir a matriz

#********************************************************************************#

#Matriz de correlação (Tabela 4.4)
corrY1=matrix(0,4,4)
corrY1[1,2]=cor(log(dadosY1[,1]/dadosY1[,3]),log(dadosY1[,2]/dadosY1[,3]))
corrY1[1,3]=cor(log(dadosY1[,1]/dadosY1[,3]),log(dadosY1[,3]/dadosY1[,4]))
corrY1[1,4]=cor(log(dadosY1[,1]/dadosY1[,3]),log(dadosY1[,4]/dadosY1[,5]))

corrY1[2,3]=cor(log(dadosY1[,2]/dadosY1[,3]),log(dadosY1[,3]/dadosY1[,4]))
corrY1[2,4]=cor(log(dadosY1[,2]/dadosY1[,3]),log(dadosY1[,4]/dadosY1[,5]))

corrY1[3,4]=cor(log(dadosY1[,3]/dadosY1[,4]),log(dadosY1[,4]/dadosY1[,5]))

corrY1 #Imprimir a matriz

#------------------------------------------------------------------------------#
# Cartogramas ####

#Grupo de outliers final
outY1Final = c(3,11,34,63,67,70,85,86,100,103,104,121,132,138,141,144,145,148,
               156,159,188,200,212,233,267,270,284,294,296,307)

#Conversão dos nomes dos municípios da matriz (comum a todos os cartogramas)
Y1names[,1] = stri_trans_general(Y1names[,1], "Latin-ASCII")
Y1names[,1] = toupper(Y1names[,1]) #Colocar em letras maiúsculas

#********************************************************************************#
# Portugal ####

#Leitura do ficheiro shp de dados geoespaciais de Portugal
portugalY1 = readOGR(file.choose())
plot(portugalY1)
View(portugalY1)
str(portugalY1@data)

#Conversão dos nomes dos municípios de Portugal (sem caracteres estranhos)
portugalY1@data$Concelho = iconv(portugalY1@data$Concelho, from = "UTF-8", to = "ASCII//TRANSLIT")

#Construir o cartograma
cartPortY1 = data.frame(portugalY1@data$Concelho)
cartPortY1$classe = 0 #Acrescentar a coluna das classes
colnames(cartPortY1) = c("Municípios","Classe")
str(cartPortY1)

#Algoritmo que atribui à coluna 'classe' do cartograma os outliers com 1's
for(j in outY1Final)
{
  for(i in 1:3223)
  {
    if (as.logical(as.character(cartPortY1[i,1])==as.character(Y1names[j,1])))
    {
      cartPortY1[i,2]=1
    }
  }
}

#Construir o mapa de Portugal
base = data.frame(Municipios = portugalY1@data$Concelho, Classe = cartPortY1$Classe) #Data frame do cartograma
head(base) #Verifica-se as classes para cada município
portugalY1@data$X = as.factor(base$Classe)
spplot(portugalY1, "X", col.regions = c("black","#66CD00"), main = "Portugal (Continental)",
       colorkey = list(space = "left", height = 0.4), col = NA)

#********************************************************************************#
# Madeira ####

#Leitura do ficheiro shp de dados geoespaciais da Madeira
madeiraY1 = readOGR(file.choose())
plot(madeiraY1)
View(madeiraY1)
str(madeiraY1@data)

#Conversão dos nomes dos municípios da Madeira (sem caracteres estranhos)
madeiraY1@data$Concelho = iconv(madeiraY1@data$CONCELHO, from = "UTF-8", to = "ASCII//TRANSLIT")

#Construir o cartograma
cartMadY1= data.frame(madeiraY1@data$CONCELHO)
cartMadY1$classe = 0
colnames(cartMadY1) = c("Municípios","Classe")
str(cartMadY1)

#Algoritmo que atribui à coluna 'classe' do cartograma os outliers com 1's
for(j in outY1Final)
{
  for(i in 1:538)
  {
    if (as.logical(as.character(cartMadY1[i,1])==as.character(Y1names[j,1])))
    {
      cartMadY1[i,2]=1
    }
  }
}

#Construir o mapa da Madeira
base = data.frame(Municipios = madeiraY1@data$CONCELHO, Classe = cartMadY1$Classe)
head(base)
madeiraY1@data$X = as.factor(base$Classe)
spplot(madeiraY1, "X", col.regions = c("black","#66CD00"), main = "Madeira",
       colorkey = list(space="left", height = 0.4), col = NA, xlim=c(287859.9,389570.9),
       ylim=c(3580000,3665000))

#********************************************************************************#
# Açores/OCIDENTAL ####

#Leitura do ficheiro shp de dados geoespaciais dos Açores/Ocidental
acoresOciY1 = readOGR(file.choose())
plot(acoresOciY1)
View(acoresOciY1)
str(acoresOciY1@data)

#Conversão dos nomes dos municípios dos Açores/Ocidental (sem caracteres estranhos)
acoresOciY1@data$Concelho = iconv(acoresOciY1@data$Concelho, from = "UTF-8", to = "ASCII//TRANSLIT")

#Construir o cartograma
cartAcorOciY1 = data.frame(acoresOciY1@data$Concelho)
cartAcorOciY1$classe = 0
colnames(cartAcorOciY1) = c("Municípios","Classe")
str(cartAcorOciY1)

#Algoritmo que atribui à coluna 'classe' do cartograma os outliers com 1's
for(j in outY1Final)
{
  for(i in 1:12)
  {
    if (as.logical(as.character(cartAcorOciY1[i,1])==as.character(Y1names[j,1])))
    {
      cartAcorOciY1[i,2]=1
    }
  }
}

#Construir o mapa dos Açores/Ocidental
base = data.frame(Municipios = acoresOciY1@data$Concelho, Classe = cartAcorOciY1$Classe)
head(base)
acoresOciY1@data$X = as.factor(base$Classe)
spplot(acoresOciY1, "X", col.regions = c("black","#66CD00"), main = "Açores (Ocidental)",
       colorkey = list(space="left", height = 0.4), col = NA)

#********************************************************************************#
# Açores/CENTRAL ####

#Leitura do ficheiro shp de dados geoespaciais dos Açores/Central
acoresCentY1 = readOGR(file.choose())
plot(acoresCentY1)
View(acoresCentY1)
str(acoresCentY1@data)

#Conversão dos nomes dos municípios dos Açores/Central (sem caracteres estranhos)
acoresCentY1@data$Concelho = iconv(acoresCentY1@data$Concelho, from = "UTF-8", to = "ASCII//TRANSLIT")

#Construir o cartograma
cartAcorCentY1 = data.frame(acoresCentY1@data$Concelho)
cartAcorCentY1$classe = 0
colnames(cartAcorCentY1) = c("Municípios","Classe")
str(cartAcorCentY1)

#Algoritmo que atribui à coluna 'classe' do cartograma os outliers com 1's
for(j in outY1Final)
{
  for(i in 1:75)
  {
    if (as.logical(as.character(cartAcorCentY1[i,1])==as.character(Y1names[j,1])))
    {
      cartAcorCentY1[i,2]=1
    }
  }
}

#Construir o mapa dos Açores/Central
base = data.frame(Municipios = acoresCentY1@data$Concelho, Classe = cartAcorCentY1$Classe)
head(base)
acoresCentY1@data$X = as.factor(base$Classe)

#NOTA: este 'if' é importante para bases de dados onde não haja outliers
if(length(unique(acoresCentY1@data$X))==1)
{
  spplot(acoresCentY1, "X", col.regions = c("black"), main = "Açores (Central)",
         colorkey = list(space="left", height = 0.4), col = NA)
} else
{
  spplot(acoresCentY1, "X", col.regions = c("black","#66CD00"), main = "Açores (Central)",
         colorkey = list(space="left", height = 0.4), col = NA)
}

#********************************************************************************#
# Açores/ORIENTAL ####

#Leitura do ficheiro shp de dados geoespaciais dos Açores/Oriental
acoresOriY1 = readOGR(file.choose())
plot(acoresOriY1)
View(acoresOriY1)
str(acoresOriY1@data)

#Conversão dos nomes dos municípios dos Açores/Oriental (sem caracteres estranhos)
acoresOriY1@data$Concelho = iconv(acoresOriY1@data$Concelho, from = "UTF-8", to = "ASCII//TRANSLIT")

#Construir o cartograma
cartAcorOriY1 = data.frame(acoresOriY1@data$Concelho)
cartAcorOriY1$classe = 0
colnames(cartAcorOriY1) = c("Municípios","Classe")
str(cartAcorOriY1)

#Algoritmo que atribui à coluna 'classe' do cartograma os outliers com 1's
for(j in outY1Final)
{
  for(i in 1:70)
  {
    if (as.logical(as.character(cartAcorOriY1[i,1])==as.character(Y1names[j,1])))
    {
      cartAcorOriY1[i,2]=1
    }
  }
}

#Construir o mapa dos Açores/Oriental
base = data.frame(Municipios = acoresOriY1@data$Concelho, Classe = cartAcorOriY1$Classe)
head(base)
acoresOriY1@data$X = as.factor(base$Classe)

#NOTA: este 'if' é importante para bases de dados onde não haja outliers
if(length(unique(acoresOriY1@data$X))==1)
{
  spplot(acoresOriY1, "X", col.regions = c("black"), main = "Açores (Oriental)",
         colorkey = list(space="left", height = 0.4), col = NA)
} else
{
  spplot(acoresOriY1, "X", col.regions = c("black","#66CD00"), main = "Açores (Oriental)",
         colorkey = list(space="left", height = 0.4), col = NA)
}

###################***************************************###################
###################       Matriz HabilitaçõesI 2005       ###################
###################***************************************###################

# Leitura da matriz ####
Y2inicial = readWorksheet(ficheiro2005, sheet="MatrizHabilitaçõesI2005")
Y2 = Y2inicial
GrandTotal = Y2[-309, 12]
Y2 = Y2[-309,-c(1,12)] #Retirar 'GrandTotal'

dadosY2 = Y2/GrandTotal  #Distribuição das habilitações (académicas) dos habitantes internos que entraram em cada município
colnames(dadosY2) = c("Bacharelato","Doutoramento","EB 1ºCiclo","EB 2ºCiclo","EB 3ºCiclo",
                      "Ens Pos Secund", "Ens Secun", "Licenciatura","Mestrado","Nenhum")

#Variável que guarda só os nomes dos municípios
Y2names = Y2inicial[-309,-12] #Adicionar a coluna dos nomes e retirar a coluna do 'GrandTotal'
colnames(Y2names) = c("Municípios","Bacharelato","Doutoramento","EB 1ºCiclo","EB 2ºCiclo",
                      "EB 3ºCiclo","Ens Pos Secund", "Ens Secun", "Licenciatura","Mestrado","Nenhum")

#Método k-NN para atribuir valores aos espaços vazios
k = impKNNa(dadosY2, method = "knn", k = 3, metric = "Aitchison",
            agg = "median", primitive = FALSE, normknn = TRUE, das = TRUE, adj = "median")
str(impKNNa(dadosY2))
dadosY2C = impKNNa(dadosY2)$xImp #Matriz de dados completa (espaços preenchidos)

#------------------------------------------------------------------------------#
# Metodologias NUMÉRICAS ####

### PBS segundo Filzmoser ####
dadosY2ilr = pivotCoord(dadosY2C) #Coordenadas ilr-transformadas 
names(dadosY2ilr) = paste0("z", 1:ncol(dadosY2ilr))
Base = orthbasis(ncol(dadosY2)) #Base ortonormal
matrizContr = Base$V #Matriz de contrastes (colunas ortonormais)
PBS = Base$basisv #Tabela da PBS

#********************************************************************************#

### Estimador MCD ####
robustezY2 = covMcd(dadosY2ilr, alpha = 0.75, nsamp = "deterministic") #Algoritmo Fast-MCD
#NOTA: obtém-se as estimativas robustas da média e da matriz de covariância
summary(robustezY2)

dY2 = determinant(robustezY2$cov,logarithm = TRUE)$modulus[1] #-24.58819
ddY2 = determinant(robustezY2$cov,logarithm = FALSE)$modulus[1] #2.096441e-11
#NOTA: ln(2.096441e-11)=-24.58819

#Outliers de 'covMcd'
outCovY2 = which(robustezY2$mcd.wt==0) #69 outliers
Y2names[outCovY2,1] #Nomes dos municípios outliers

#Gráfico de 'covMcd'
plot(robustezY2, which = "distance", classic = TRUE)

#********************************************************************************#

### Distância de Mahalanobis ####
#Adjusted Quantile Plot
#Função aq.plot original alterada
aq.plot = function (x, delta = qchisq(0.975, df = ncol(x)), quan = 1/2,alpha = 0.05) 
{
  if (is.vector(x) == TRUE || ncol(x) == 1) {
    stop("x must be at least two-dimensional")
  }
  covr <- covMcd(x, alpha = quan)
  mcd <- determinant(covr$cov, logarithm = TRUE)$modulus[1] #Acrescentei esta linha
  dist <- mahalanobis(x, center = covr$center, cov = covr$cov)
  s <- sort(dist, index = TRUE)
  z <- x
  if (ncol(x) > 2) {
    p <- princomp(x, covmat = covr)
    z <- p$scores[, 1:2]
    sdprop <- (p$sd[1] + p$sd[2])/sum(p$sd)
    cat("Projection to the first and second robust principal components.\\n")
    cat("Proportion of total variation (explained variance): ")
    cat(sdprop)
    cat("\\n")
  }
  par(mfrow = c(2, 2), mai = c(0.8, 0.6, 0.2, 0.2), mgp = c(2.4,1, 0))
  plot(z, col = 3, type = "n", xlab = "", ylab = "")
  text(z, dimnames(as.data.frame(z))[[1]], col = 3, cex = 0.8)
  plot(s$x, (1:length(dist))/length(dist), col = 3, xlab = "Ordered squared robust distance", 
       ylab = "Cumulative probability", type = "n")
  text(s$x, (1:length(dist))/length(dist), as.character(s$ix),col = 3, cex = 0.8)
  t <- seq(0, max(dist), by = 0.01)
  lines(t, pchisq(t, df = ncol(x)), col = 6)
  abline(v = delta, col = 5)
  text(x = delta, y = 0.4, paste(100 * (pchisq(delta, df = ncol(x))), 
                                 "% Quantile", sep = ""), col = 5, pos = 2, srt = 90,cex = 0.8)
  xarw <- arw(x, covr$center, covr$cov, alpha = alpha)
  if (xarw$cn < Inf) {
    abline(v = xarw$cn, col = 4)
    text(x = xarw$cn, y = 0.4, "Adjusted Quantile", col = 4,pos = 4, srt = 90, cex = 0.8)
  }
  plot(z, col = 3, type = "n",
       main = paste("Outliers based on ",100 * (pchisq(delta, df = ncol(x))), "% quantile", sep = ""), 
       xlab = "", ylab = "")
  if (any(dist > delta)) {
    text(z[dist > delta, 1], z[dist > delta, 2],
         dimnames(as.data.frame(x)[dist > delta, ])[[1]], col = 2, cex = 0.8)
  }
  if (any(dist <= delta)) {
    text(z[dist <= delta, 1], z[dist <= delta, 2],
         dimnames(as.data.frame(x)[dist <= delta, ])[[1]], col = 3, cex = 0.8)
  }
  plot(z, col = 3, type = "n", main = "Outliers based on adjusted quantile", 
       xlab = "", ylab = "")
  if (xarw$cn < Inf) {
    text(z[dist > xarw$cn, 1], z[dist > xarw$cn, 2],
         dimnames(as.data.frame(x)[dist > xarw$cn, ])[[1]], col = 2, cex = 0.8)
  }
  text(z[dist <= xarw$cn, 1], z[dist <= xarw$cn, 2],
       dimnames(as.data.frame(x)[dist <= xarw$cn, ])[[1]], col = 3, cex = 0.8)
  
  o <- (sqrt(dist) > max(sqrt(xarw$cn), sqrt(qchisq(0.975, dim(x)[2]))))
  aqplot <- list(outliers = o, MCD1 = mcd) #Acrescentei o último parâmetro
  return(aqplot)
}

#Usar a função aqplot alterada
adjquanY2 = aq.plot(dadosY2ilr, delta = qchisq(0.975, df = ncol(dadosY2ilr)),
                    quan = 0.75, alpha = 0.975)

#Valor do estimador MCD
adjquanY2$MCD1

#Outliers
outAQplotY2 = which(adjquanY2$outliers == TRUE) #Obtém-se as observações que são outliers
length(outAQplotY2)
Y2names[outAQplotY2,1] #Nomes dos municípios outliers

#SOLUÇÃO: guardar os outliers associados ao menor valor de MCD
out1Y2Final = c(3,4,11,13,18,20,24,36,50,57,61,63,67,70,71,72,74,78,83,85,92,97,103,104,107,121,
                122,131,138,141,147,156,157,158,164,166,172,174,177,188,194,197,200,203,206,208,
                212,216,222,225,228,244,248,253,262,282,284,285,290,296,297,305)
# 62 outliers com MCD=-24.61327

#Nomes dos municípios outliers para a solução
Y2names[c(3,4,11,13,18,20,24,36,50,57,61,63,67,70,71,72,74,78,83,85,92,97,103,104,107,121,
          122,131,138,141,147,156,157,158,164,166,172,174,177,188,194,197,200,203,206,208,
          212,216,222,225,228,244,248,253,262,282,284,285,290,296,297,305),1]

#................................................................................#

#Distance-Distance Plot
#Função ddplot original alterada
dd.plot = function (x, quan = 1/2, alpha = 0.025, ...) 
{
  if (!is.matrix(x) && !is.data.frame(x)) 
    stop("x must be matrix or data.frame")
  rob <- covMcd(x, alpha = quan)
  mcd <- determinant(rob$cov, logarithm = TRUE)$modulus[1] #Acrescentei esta linha
  xarw <- arw(x, rob$center, rob$cov, alpha = alpha)
  distcla <- sqrt(mahalanobis(x, center = apply(x, 2, mean), 
                              cov = cov(x)))
  distrob <- sqrt(mahalanobis(x, center = rob$center, cov = rob$cov))
  plot(distcla, distrob, main = "Distance-Distance Plot", xlab = "Mahalanobis Distance", 
       ylab = "Robust Distance", type = "n", ...)
  if (xarw$cn != Inf) {
    alpha <- sqrt(c(xarw$cn, qchisq(c(0.75, 0.5, 0.25), ncol(x))))
  }
  else {
    alpha <- sqrt(qchisq(c(0.975, 0.75, 0.5, 0.25), ncol(x)))
  }
  abline(h = alpha[1])
  abline(v = alpha[1])
  abline(a = 0, b = 1)
  lpch <- c(3, 3, 16, 1, 1)
  lcex <- c(1.5, 1, 0.5, 1, 1.5)
  lalpha <- length(alpha)
  xs <- scale(x) - min(scale(x))
  eucl <- sqrt(apply(xs^2, 1, sum))
  rbcol <- rev(rainbow(nrow(x), start = 0,
                       end = 0.7))[as.integer(cut(eucl, nrow(x), labels = 1:nrow(x)))]
  rd <- distrob
  for (j in 1:lalpha) {
    if (j == 1) {
      points(distcla[rd >= alpha[j]], distrob[rd >= alpha[j]], 
             pch = lpch[j], cex = lcex[j], col = rbcol[rd >= alpha[j]])
    }
    if (j > 1 & j < lalpha) 
      points(distcla[rd < alpha[j - 1] & rd >= alpha[j]], 
             distrob[rd < alpha[j - 1] & rd >= alpha[j]], 
             cex = lcex[j], pch = lpch[j], col = rbcol[rd < alpha[j - 1] & rd >= alpha[j]])
    if (j == lalpha) {
      points(distcla[rd < alpha[j - 1] & rd >= alpha[j]], 
             distrob[rd < alpha[j - 1] & rd >= alpha[j]], 
             cex = lcex[j], pch = lpch[j], col = rbcol[rd < alpha[j - 1] & rd >= alpha[j]])
      points(distcla[rd < alpha[j]], distrob[rd < alpha[j]], 
             pch = lpch[j + 1], cex = lcex[j + 1], col = rbcol[rd < alpha[j]])
    }
  }
  o <- (rd > min(sqrt(xarw$cn), sqrt(qchisq(0.975, dim(x)[2]))))
  ddplot <- list(outliers = o, md.cla = distcla, md.rob = distrob, MCD2 = mcd) #Acrescentei o último parâmetro
  return(ddplot)
}

#Usar a função ddplot alterada
DDplotY2 = dd.plot(dadosY2ilr, quan = 0.75, alpha = 0.975)

#Valor do estimador MCD
DDplotY2$MCD2

#Outliers
outDDplotY2 = which(DDplotY2$outliers == TRUE)
length(outDDplotY2)
Y2names[outDDplotY2,1]

#SOLUÇÃO: guardar os outliers associados ao menor valor de MCD
out2Y2Final = c(3,4,11,13,18,20,24,36,39,50,57,61,63,67,70,71,72,74,78,83,85,92,97,103,104,107,
                121,122,131,138,141,147,156,157,158,164,166,171,172,174,177,188,194,197,200,203,
                206,208,212,216,217,222,225,228,230,243,244,248,253,262,264,282,283,284,285,290,
                296,297,305)
# 69 outliers com MCD=-24.58819

#Nomes dos municípios outliers para a solução
Y2names[c(3,4,11,13,18,20,24,36,39,50,57,61,63,67,70,71,72,74,78,83,85,92,97,103,104,107,
          121,122,131,138,141,147,156,157,158,164,166,171,172,174,177,188,194,197,200,203,
          206,208,212,216,217,222,225,228,230,243,244,248,253,262,264,282,283,284,285,290,
          296,297,305),1]

#................................................................................#

#Outlier CoDa
#Função outCoDa original alterada
outCoDa = function (x, quantile = 0.975, method = "robust", h = 1/2, coda = TRUE) 
{
  if (dim(x)[2] < 2) 
    stop("need data with at least 2 variables")
  covEst <- function(x, type) {
    standard <- function(x) {
      list(mean = colMeans(x, na.rm = TRUE), varmat = cov(x))
    }
    robust <- function(x) {
      v <- robustbase::covMcd(x)
      mcd <- base::determinant(v$cov, logarithm = TRUE)$modulus[1] #Acrescentei esta linha
      list(mean = v$center, varmat = v$cov, MCD = mcd) #Acrescentei o último parâmetro
    }
    switch(type, standard = standard(x), robust = robust(x))
  }
  if (!is.logical(coda) & !is.function(coda)) {
    stop("coda must be logical or function")
  }
  if (!is.logical(coda)) {
    x <- coda(x)
  }
  else if (coda) {
    x <- pivotCoord(x)
  }
  cv <- covEst(x, "robust")
  cvc <- covEst(x, "standard")
  dM <- sqrt(mahalanobis(x, center = cv$mean, cov = cv$varmat))
  dMc <- sqrt(mahalanobis(x, center = cvc$mean, cov = cvc$varmat))
  limit <- sqrt(qchisq(p = quantile, df = ncol(x) - 1))
  res <- list(mahalDist = dM, limit = limit, outlierIndex = dM > limit, method = method,
              om2 = dMc > limit, m2 = dMc, coda = coda, MCD3 = cv$MCD) #Acrescentei o último parâmetro
  class(res) <- "outCoDa"
  invisible(res)
}

#Usar a função outCoDa alterada
outRobY2 = outCoDa(dadosY2C, quantile = 0.975, method = "robust", h = 231)
str(outRobY2)

#Valor do estimador MCD
outRobY2$MCD3

#Outliers
outROBY2 = which(outRobY2$outlierIndex == TRUE)
Y2names[outROBY2,1]

#SOLUÇÃO: guardar os outliers associados ao menor valor de MCD
out3Y2Final = c(3,4,11,13,18,20,23,24,36,39,41,47,50,51,57,61,63,67,70,71,72,73,74,75,78,83,85,
                92,97,100,102,103,104,107,108,116,121,122,131,138,141,146,147,156,157,158,159,
                164,166,167,171,172,174,177,184,188,191,194,197,200,203,206,208,212,216,217,
                222,225,226,228,230,243,244,248,253,262,264,276,279,282,283,284,285,290,296,
                297,301,305)
# 88 outliers com MCD=-25.4167

#Nomes dos municípios outliers para a solução
Y2names[c(3,4,11,13,18,20,23,24,36,39,41,47,50,51,57,61,63,67,70,71,72,73,74,75,78,83,85,
          92,97,100,102,103,104,107,108,116,121,122,131,138,141,146,147,156,157,158,159,
          164,166,167,171,172,174,177,184,188,191,194,197,200,203,206,208,212,216,217,
          222,225,226,228,230,243,244,248,253,262,264,276,279,282,283,284,285,290,296,
          297,301,305),1]

#------------------------------------------------------------------------------#
# Metodologias GRÁFICAS ####

dadosY2 = data.frame(dadosY2C)
colnames(dadosY2) = c("Bacharelato","Doutoramento","EB 1ºCiclo","EB 2ºCiclo","EB 3ºCiclo",
                      "Ens Pos Secund", "Ens Secun", "Licenciatura","Mestrado","Nenhum")

#Função mvoutlier.CoDa original alterada
mvoutlier.CoDa = function (x, quan = 0.75, alpha = 0.025,
                           col.quantile = c(0,0.05, 0.1, 0.5, 0.9, 0.95, 1),
                           symb.pch = c(3, 3, 16, 1,1), symb.cex = c(1.5, 1, 0.5, 1, 1.5),
                           adaptive = TRUE)
{
  if (!is.matrix(x) && !is.data.frame(x)) 
    stop("x must be matrix or data.frame")
  if (ncol(x) < 3) 
    stop("x must have at least 3 compositional parts")
  Z <- pivotCoord(x)
  V <- (orthbasis(ncol(x))$V)
  Zj <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  for (j in 1:ncol(x)) {
    Zj[, j] <- pivotCoord(cbind(x[, j], x[, -j]))[, 1]
  }
  dimnames(Zj)[[2]] <- names(x)
  rob <- covMcd(Z, alpha = quan)
  mcd <- determinant(rob$cov, logarithm = TRUE)$modulus[1] #Acrescentei esta linha
  if (adaptive) {
    Zarw <- arw(Z, rob$center, rob$cov, alpha = alpha)
    if (Zarw$cn != Inf) {
      alpha1 <- sqrt(c(Zarw$cn, qchisq(c(0.75, 0.5, 0.25), ncol(Z))))
    }
    else {
      alpha1 <- sqrt(qchisq(c(0.975, 0.75, 0.5, 0.25), ncol(Z)))
    }
    rd2 <- mahalanobis(Z, center = Zarw$m, cov = Zarw$c)
  }
  else {
    cutoff <- qchisq(1 - alpha, ncol(Z))
    rd2 <- mahalanobis(Z, center = rob$center, cov = rob$cov)
    Zarw <- list(m = rob$center, c = rob$cov, cn = cutoff, w = as.numeric(rd2 < cutoff))
    alpha1 <- sqrt(qchisq(c(0.975, 0.75, 0.5, 0.25), ncol(Z)))
  }
  rd <- sqrt(rd2)
  covobj <- list(center = Zarw$m, cov = Zarw$c, n.obs = length(rd), mah = rd)
  Z.pca <- suppressWarnings(princomp(Z, covmat = covobj, cor = FALSE))
  pcaclr <- Z.pca
  eval <- eigen(Zarw$c)$values
  pcaclr$sdev <- sqrt(eval)
  pcaclr$scores <- Z.pca$scores
  pcaclr$loadings <- V %*% Z.pca$loadings
  dimnames(pcaclr$loadings)[[1]] <- names(x)
  pcaobj <- list(method = "robust", eigenvalues = eval, princompOutputClr = pcaclr)
  class(pcaobj) <- "pcaCoDa"
  Zcent <- scale(Zj, center = apply(Zj, 2, median), scale = FALSE)
  eucl <- apply(abs(Zcent), 1, median)
  out <- (!Zarw$w)
  lq <- length(col.quantile)
  colcol <- rev(rainbow(lq - 1, start = 0,
                        end = 0.7))[as.integer(cut(eucl, quantile(eucl, col.quantile, labels = 1:(lq - 1))))]
  colbw <- rev(gray(seq(from = 0.1, to = 0.9,
                        length = lq - 1)))[as.integer(cut(eucl, quantile(eucl, col.quantile, labels = 1:(lq - 1))))]
  pchvec <- rep(symb.pch[1], nrow(Zj))
  cexvec <- rep(symb.cex[1], nrow(Zj))
  if (length(symb.pch) == 5 & length(symb.cex) == 5) {
    lalpha <- length(alpha1)
    for (j in 1:(lalpha)) {
      pchvec[rd < alpha1[j]] <- symb.pch[j + 1]
      cexvec[rd < alpha1[j]] <- symb.cex[j + 1]
    }
  }
  mvoutlierCoDa <- list(ilrvariables = Zj, outliers = out, pcaobj = pcaobj, colcol = colcol,
                        colbw = colbw, pchvec = pchvec, cexvec = cexvec,
                        media = rob$center, cov = rob$cov, MCD4 = mcd)
  #Acrescentei os três últimos parâmetros
  class(mvoutlierCoDa) <- "mvoutlierCoDa"
  return(mvoutlierCoDa)
}

#Usar a função mvoutlierCoDa alterada
resY2 = mvoutlier.CoDa(dadosY2, quan = 0.75, alpha = 0.975, adaptive = TRUE)
str(resY2)

#Valor do estimador MCD
resY2$MCD4

#Outliers
outResY2 = which(resY2$outliers == TRUE)
length(outResY2)
Y2names[outResY2,1]

#Matriz de covariância
matrizCovY2 = resY2$cov
colnames(matrizCovY2) = c("z1","z2","z3","z4","z5","z6","z7","z8","z9")
row.names(matrizCovY2) = c("z1","z2","z3","z4","z5","z6","z7","z8","z9")
matrizCovY2

#Vetor da média aritmética
resY2$media

#................................................................................#

#SOLUÇÃO: guardar os outliers associados ao menor valor de MCD
out4Y2Final = c(3,4,11,13,18,20,24,36,50,57,61,63,67,70,71,72,74,78,83,85,92,97,103,104,
                107,121,122,131,138,141,147,156,157,158,164,166,171,172,174,177,188,194,
                197,200,203,206,208,212,216,222,225,228,244,248,253,262,282,284,285,290,
                296,297,305)
# 63 outliers com MCD=-24.58819

#Nomes dos municípios outliers para a solução
Y2names[c(3,4,11,13,18,20,24,36,50,57,61,63,67,70,71,72,74,78,83,85,92,97,103,104,
          107,121,122,131,138,141,147,156,157,158,164,166,171,172,174,177,188,194,
          197,200,203,206,208,212,216,222,225,228,244,248,253,262,282,284,285,290,
          296,297,305),1]

#Matriz de covariância associada ao menor valor de MCD
matrizCovY2Final = read.table("MatrizCovarianciaY2.txt", sep = "", dec = ".", header = TRUE)
# para MCD=-24.58819

#Média aritmética associada ao menor valor de MCD
mediaY2Final = c(-1.0643384,-3.2734310,0.8506010,0.6000983,1.0428625,-1.5563474,
                 0.7914861,0.8737796,-1.5389530)
# para MCD=-24.58819

#Novamente a função anterior, mas alterou-se o valor de MCD e as estimativas robustas dos parâmetros \mu e \Sigma
mvoutlier.CoDa = function (x, quan = 0.75, alpha = 0.025,
                           col.quantile = c(0,0.05, 0.1, 0.5, 0.9, 0.95, 1),
                           symb.pch = c(3, 3, 16, 1,1), symb.cex = c(1.5, 1, 0.5, 1, 1.5),
                           adaptive = TRUE)
{
  if (!is.matrix(x) && !is.data.frame(x)) 
    stop("x must be matrix or data.frame")
  if (ncol(x) < 3) 
    stop("x must have at least 3 compositional parts")
  Z <- pivotCoord(x)
  V <- (orthbasis(ncol(x))$V)
  Zj <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  for (j in 1:ncol(x)) {
    Zj[, j] <- pivotCoord(cbind(x[, j], x[, -j]))[, 1]
  }
  dimnames(Zj)[[2]] <- names(x)
  matrizCovY2Final <- read.table("MatrizCovarianciaY2.txt", sep = "", dec = ".", header = TRUE) #Acrescentei
  mediaY2Final <- c(-1.0643384,-3.2734310,0.8506010,0.6000983,1.0428625,-1.5563474,
                    0.7914861,0.8737796,-1.5389530) #Acrescentei
  mcd <- -24.58819 #Coloquei o valor fixo
  if (adaptive) {
    Zarw <- arw(Z, mediaY2Final, matrizCovY2Final, alpha = alpha)
    if (Zarw$cn != Inf) {
      alpha1 <- sqrt(c(Zarw$cn, qchisq(c(0.75, 0.5, 0.25), ncol(Z))))
    }
    else {
      alpha1 <- sqrt(qchisq(c(0.975, 0.75, 0.5, 0.25), ncol(Z)))
    }
    rd2 <- mahalanobis(Z, center = Zarw$m, cov = Zarw$c)
  }
  else {
    cutoff <- qchisq(1 - alpha, ncol(Z))
    rd2 <- mahalanobis(Z, center = mediaY2Final, cov = matrizCovY2Final)
    Zarw <- list(m = mediaY2Final, c = matrizCovY2Final, cn = cutoff, w = as.numeric(rd2 < cutoff))
    alpha1 <- sqrt(qchisq(c(0.975, 0.75, 0.5, 0.25), ncol(Z)))
  }
  rd <- sqrt(rd2)
  covobj <- list(center = Zarw$m, cov = Zarw$c, n.obs = length(rd), mah = rd)
  Z.pca <- suppressWarnings(princomp(Z, covmat = covobj, cor = FALSE))
  pcaclr <- Z.pca
  eval <- eigen(Zarw$c)$values
  pcaclr$sdev <- sqrt(eval)
  pcaclr$scores <- Z.pca$scores
  pcaclr$loadings <- V %*% Z.pca$loadings
  dimnames(pcaclr$loadings)[[1]] <- names(x)
  pcaobj <- list(method = "robust", eigenvalues = eval, princompOutputClr = pcaclr)
  class(pcaobj) <- "pcaCoDa"
  Zcent <- scale(Zj, center = apply(Zj, 2, median), scale = FALSE)
  eucl <- apply(abs(Zcent), 1, median)
  out <- (!Zarw$w)
  lq <- length(col.quantile)
  colcol <- rev(rainbow(lq - 1, start = 0,
                        end = 0.7))[as.integer(cut(eucl, quantile(eucl, col.quantile, labels = 1:(lq - 1))))]
  colbw <- rev(gray(seq(from = 0.1, to = 0.9,
                        length = lq - 1)))[as.integer(cut(eucl, quantile(eucl, col.quantile, labels = 1:(lq - 1))))]
  pchvec <- rep(symb.pch[1], nrow(Zj))
  cexvec <- rep(symb.cex[1], nrow(Zj))
  if (length(symb.pch) == 5 & length(symb.cex) == 5) {
    lalpha <- length(alpha1)
    for (j in 1:(lalpha)) {
      pchvec[rd < alpha1[j]] <- symb.pch[j + 1]
      cexvec[rd < alpha1[j]] <- symb.cex[j + 1]
    }
  }
  mvoutlierCoDa <- list(ilrvariables = Zj, outliers = out, pcaobj = pcaobj, colcol = colcol,
                        colbw = colbw, pchvec = pchvec, cexvec = cexvec,
                        cov = matrizCovY2Final, MCD4 = mcd, media = mediaY2Final)
  #Acrescentei os três últimos parâmetros
  class(mvoutlierCoDa) <- "mvoutlierCoDa"
  return(mvoutlierCoDa)
}

resY2Alt = mvoutlier.CoDa(dadosY2, quan = 0.75, alpha = 0.975, adaptive = TRUE)
str(resY2Alt)

#Confirmar os resultados
resY2Alt$cov
resY2Alt$MCD4
resY2Alt$media

#Outliers (confirmar)
outResY2Alt = which(resY2Alt$outliers == TRUE)
length(outResY2Alt)
Y2names[outResY2Alt,1]

#................................................................................#

#Biplot: para visualização dos outliers
plot(resY2Alt, which = "biplot", onlyout = TRUE, symb = TRUE, symbtxt = TRUE)
#Apenas os outliers são exibidos e identificados por números
#NOTA: os outliers estão ordenados

#Biplot: para visualização dos ângulos entre os links
biplotAngY2 = pcaCoDa(dadosY2, method = "robust")
summary(biplotAngY2)
biplot(biplotAngY2, xlabs = rep(".", nrow(dadosY2)))

#................................................................................#

#Gráfico de dispersão univariado
plot(resY2Alt, which = "uni", onlyout = TRUE, symb = TRUE, symbtxt = TRUE)
#Apenas os outliers são exibidos e identificados por números

#********************************************************************************#

#Matriz de variação (Tabela 4.8)
dadosY2CTabVar = dadosY2C[,c(10,3,4,5,7,6,1,8,9,2)] #Reordenar as colunas da matriz
varY2 = matrix(0,10,10)
varY2[1,2] = var(log(dadosY2CTabVar[,1]/dadosY2CTabVar[,2])); varY2[1,3] = var(log(dadosY2CTabVar[,1]/dadosY2CTabVar[,3]))
varY2[1,4] = var(log(dadosY2CTabVar[,1]/dadosY2CTabVar[,4])); varY2[1,5] = var(log(dadosY2CTabVar[,1]/dadosY2CTabVar[,5]))
varY2[1,6] = var(log(dadosY2CTabVar[,1]/dadosY2CTabVar[,6])); varY2[1,7] = var(log(dadosY2CTabVar[,1]/dadosY2CTabVar[,7]))
varY2[1,8] = var(log(dadosY2CTabVar[,1]/dadosY2CTabVar[,8])); varY2[1,9] = var(log(dadosY2CTabVar[,1]/dadosY2CTabVar[,9]))
varY2[1,10] = var(log(dadosY2CTabVar[,1]/dadosY2CTabVar[,10]))

varY2[2,3] = var(log(dadosY2CTabVar[,2]/dadosY2CTabVar[,3])); varY2[2,4] = var(log(dadosY2CTabVar[,2]/dadosY2CTabVar[,4]))
varY2[2,5] = var(log(dadosY2CTabVar[,2]/dadosY2CTabVar[,5])); varY2[2,6] = var(log(dadosY2CTabVar[,2]/dadosY2CTabVar[,6]))
varY2[2,7] = var(log(dadosY2CTabVar[,2]/dadosY2CTabVar[,7])); varY2[2,8] = var(log(dadosY2CTabVar[,2]/dadosY2CTabVar[,8]))
varY2[2,9] = var(log(dadosY2CTabVar[,2]/dadosY2CTabVar[,9])); varY2[2,10] = var(log(dadosY2CTabVar[,2]/dadosY2CTabVar[,10]))

varY2[3,4] = var(log(dadosY2CTabVar[,3]/dadosY2CTabVar[,4])); varY2[3,5] = var(log(dadosY2CTabVar[,3]/dadosY2CTabVar[,5]))
varY2[3,6] = var(log(dadosY2CTabVar[,3]/dadosY2CTabVar[,6])); varY2[3,7] = var(log(dadosY2CTabVar[,3]/dadosY2CTabVar[,7]))
varY2[3,8] = var(log(dadosY2CTabVar[,3]/dadosY2CTabVar[,8])); varY2[3,9] = var(log(dadosY2CTabVar[,3]/dadosY2CTabVar[,9]))
varY2[3,10] = var(log(dadosY2CTabVar[,3]/dadosY2CTabVar[,10]))

varY2[4,5] = var(log(dadosY2CTabVar[,4]/dadosY2CTabVar[,5])); varY2[4,6] = var(log(dadosY2CTabVar[,4]/dadosY2CTabVar[,6]))
varY2[4,7] = var(log(dadosY2CTabVar[,4]/dadosY2CTabVar[,7])); varY2[4,8] = var(log(dadosY2CTabVar[,4]/dadosY2CTabVar[,8]))
varY2[4,9] = var(log(dadosY2CTabVar[,4]/dadosY2CTabVar[,9])); varY2[4,10] = var(log(dadosY2CTabVar[,4]/dadosY2CTabVar[,10]))

varY2[5,6] = var(log(dadosY2CTabVar[,5]/dadosY2CTabVar[,6])); varY2[5,7] = var(log(dadosY2CTabVar[,5]/dadosY2CTabVar[,7]))
varY2[5,8] = var(log(dadosY2CTabVar[,5]/dadosY2CTabVar[,8])); varY2[5,9] = var(log(dadosY2CTabVar[,5]/dadosY2CTabVar[,9]))
varY2[5,10] = var(log(dadosY2CTabVar[,5]/dadosY2CTabVar[,10]))

varY2[6,7] = var(log(dadosY2CTabVar[,6]/dadosY2CTabVar[,7])); varY2[6,8] = var(log(dadosY2CTabVar[,6]/dadosY2CTabVar[,8]))
varY2[6,9] = var(log(dadosY2CTabVar[,6]/dadosY2CTabVar[,9])); varY2[6,10] = var(log(dadosY2CTabVar[,6]/dadosY2CTabVar[,10]))

varY2[7,8] = var(log(dadosY2CTabVar[,7]/dadosY2CTabVar[,8])); varY2[7,9] = var(log(dadosY2CTabVar[,7]/dadosY2CTabVar[,9]))
varY2[7,10] = var(log(dadosY2CTabVar[,7]/dadosY2CTabVar[,10]))

varY2[8,9] = var(log(dadosY2CTabVar[,8]/dadosY2CTabVar[,9])); varY2[8,10] = var(log(dadosY2CTabVar[,8]/dadosY2CTabVar[,10]))

varY2[9,10] = var(log(dadosY2CTabVar[,9]/dadosY2CTabVar[,10]))

varY2 #Imprimir a matriz

#********************************************************************************#

#Matriz de correlação (Tabela 4.9)
corrY2 = matrix(0,6,6)
corrY2[1,2] = cor(log(dadosY2C[,7]/dadosY2C[,2]), log(dadosY2C[,5]/dadosY2C[,7]))
corrY2[1,3] = cor(log(dadosY2C[,7]/dadosY2C[,2]), log(dadosY2C[,4]/dadosY2C[,5]))
corrY2[1,4] = cor(log(dadosY2C[,7]/dadosY2C[,2]), log(dadosY2C[,10]/dadosY2C[,4]))
corrY2[1,5] = cor(log(dadosY2C[,7]/dadosY2C[,2]), log(dadosY2C[,3]/dadosY2C[,10]))
corrY2[1,6] = cor(log(dadosY2C[,7]/dadosY2C[,2]), log(dadosY2C[,5]/dadosY2C[,6]))

corrY2[2,3] = cor(log(dadosY2C[,5]/dadosY2C[,7]), log(dadosY2C[,4]/dadosY2C[,5]))
corrY2[2,4] = cor(log(dadosY2C[,5]/dadosY2C[,7]), log(dadosY2C[,10]/dadosY2C[,4]))
corrY2[2,5] = cor(log(dadosY2C[,5]/dadosY2C[,7]), log(dadosY2C[,3]/dadosY2C[,10]))
corrY2[2,6] = cor(log(dadosY2C[,5]/dadosY2C[,7]), log(dadosY2C[,5]/dadosY2C[,6]))

corrY2[3,4] = cor(log(dadosY2C[,4]/dadosY2C[,5]), log(dadosY2C[,10]/dadosY2C[,4]))
corrY2[3,5] = cor(log(dadosY2C[,4]/dadosY2C[,5]), log(dadosY2C[,3]/dadosY2C[,10]))
corrY2[3,6] = cor(log(dadosY2C[,4]/dadosY2C[,5]), log(dadosY2C[,5]/dadosY2C[,6]))

corrY2[4,5] = cor(log(dadosY2C[,10]/dadosY2C[,4]), log(dadosY2C[,3]/dadosY2C[,10]))
corrY2[4,6] = cor(log(dadosY2C[,10]/dadosY2C[,4]), log(dadosY2C[,5]/dadosY2C[,6]))

corrY2[5,6] = cor(log(dadosY2C[,3]/dadosY2C[,10]), log(dadosY2C[,5]/dadosY2C[,6]))

corrY2 #Imprimir a matriz

#------------------------------------------------------------------------------#
# Cartogramas ####

#Grupo de outliers final
outY2Final = c(3,4,11,13,18,20,24,36,50,57,61,63,67,70,71,72,74,78,83,85,92,
               97,103,104,107,121,122,131,138,141,147,156,157,158,164,166,172,
               174,177,188,194,197,200,203,206,208,212,216,222,225,228,244,248,
               253,262,282,284,285,290,296,297,305)

#Conversão dos nomes dos municípios da matriz (comum a todos os cartogramas)
Y2names[,1] = stri_trans_general(Y2names[,1], "Latin-ASCII")
Y2names[,1] = toupper(Y2names[,1]) #Colocar em letras maiúsculas

#********************************************************************************#
# Portugal ####

#Leitura do ficheiro shp de dados geoespaciais de Portugal
portugalY2 = readOGR(file.choose())
plot(portugalY2)
View(portugalY2)
str(portugalY2@data)

#Conversão dos nomes dos municípios de Portugal (sem caracteres estranhos)
portugalY2@data$Concelho = iconv(portugalY2@data$Concelho, from = "UTF-8", to = "ASCII//TRANSLIT")

#Construir o cartograma
cartPortY2 = data.frame(portugalY2@data$Concelho)
cartPortY2$classe = 0
colnames(cartPortY2) = c("Municípios","Classe")
str(cartPortY2)

#Algoritmo que atribui à coluna 'classe' do cartograma os outliers com 1's
for(j in outY2Final)
{
  for(i in 1:3223)
  {
    if (as.logical(as.character(cartPortY2[i,1])==as.character(Y2names[j,1])))
    {
      cartPortY2[i,2]=1
    }
  }
}

#Construir o mapa de Portugal
base = data.frame(Municipios = portugalY2@data$Concelho, Classe = cartPortY2$Classe)
head(base)
portugalY2@data$X = as.factor(base$Classe)
spplot(portugalY2, "X", col.regions = c("black","#66CD00"), main = "Portugal (Continental)",
       colorkey = list(space = "left", height = 0.4), col = NA)

#********************************************************************************#
# Madeira ####

#Leitura do ficheiro shp de dados geoespaciais da Madeira
madeiraY2 = readOGR(file.choose())
plot(madeiraY2)
View(madeiraY2)
str(madeiraY2@data)

#Conversão dos nomes dos municípios da Madeira (sem caracteres estranhos)
madeiraY2@data$Concelho = iconv(madeiraY2@data$CONCELHO, from = "UTF-8", to = "ASCII//TRANSLIT")

#Construir o cartograma
cartMadY2= data.frame(madeiraY2@data$CONCELHO)
cartMadY2$classe = 0
colnames(cartMadY2) = c("Municípios","Classe")
str(cartMadY2)

#Algoritmo que atribui à coluna 'classe' do cartograma os outliers com 1's
for(j in outY2Final)
{
  for(i in 1:538)
  {
    if (as.logical(as.character(cartMadY2[i,1])==as.character(Y2names[j,1])))
    {
      cartMadY2[i,2]=1
    }
  }
}

#Construir o mapa da Madeira
base = data.frame(Municipios = madeiraY2@data$CONCELHO, Classe = cartMadY2$Classe)
head(base)
madeiraY2@data$X = as.factor(base$Classe)
spplot(madeiraY2, "X", col.regions = c("black","#66CD00"), main = "Madeira",
       colorkey = list(space="left", height = 0.4), col = NA, xlim=c(287859.9,389570.9),
       ylim=c(3580000,3665000))

#********************************************************************************#
# Açores/OCIDENTAL ####

#Leitura do ficheiro shp de dados geoespaciais dos Açores/Ocidental
acoresOciY2 = readOGR(file.choose())
plot(acoresOciY2)
View(acoresOciY2)
str(acoresOciY2@data)

#Conversão dos nomes dos municípios dos Açores/Ocidental (sem caracteres estranhos)
acoresOciY2@data$Concelho = iconv(acoresOciY2@data$Concelho, from = "UTF-8", to = "ASCII//TRANSLIT")

#Construir o cartograma
cartAcorOciY2 = data.frame(acoresOciY2@data$Concelho)
cartAcorOciY2$classe = 0
colnames(cartAcorOciY2) = c("Municípios","Classe")
str(cartAcorOciY2)

#Algoritmo que atribui à coluna 'classe' do cartograma os outliers com 1's
for(j in outY2Final)
{
  for(i in 1:12)
  {
    if (as.logical(as.character(cartAcorOciY2[i,1])==as.character(Y2names[j,1])))
    {
      cartAcorOciY2[i,2]=1
    }
  }
}

#Construir o mapa dos Açores/Ocidental
base = data.frame(Municipios = acoresOciY2@data$Concelho, Classe = cartAcorOciY2$Classe)
head(base)
acoresOciY2@data$X = as.factor(base$Classe)
spplot(acoresOciY2, "X", col.regions = c("black","#66CD00"), main = "Açores (Ocidental)",
       colorkey = list(space="left", height = 0.4), col = NA)

#********************************************************************************#
# Cartograma - Açores/CENTRAL ####

#Leitura do ficheiro shp de dados geoespaciais dos Açores/Central
acoresCentY2 = readOGR(file.choose())
plot(acoresCentY2)
View(acoresCentY2)
str(acoresCentY2@data)

#Conversão dos nomes dos municípios dos Açores/Central (sem caracteres estranhos)
acoresCentY2@data$Concelho = iconv(acoresCentY2@data$Concelho, from = "UTF-8", to = "ASCII//TRANSLIT")

#Construir o cartograma
cartAcorCentY2 = data.frame(acoresCentY2@data$Concelho)
cartAcorCentY2$classe = 0
colnames(cartAcorCentY2) = c("Municípios","Classe")
str(cartAcorCentY2)

#Algoritmo que atribui à coluna 'classe' do cartograma os outliers com 1's
for(j in outY2Final)
{
  for(i in 1:75)
  {
    if (as.logical(as.character(cartAcorCentY2[i,1])==as.character(Y2names[j,1])))
    {
      cartAcorCentY2[i,2]=1
    }
  }
}

#Construir o mapa dos Açores/Central
base = data.frame(Municipios = acoresCentY2@data$Concelho, Classe = cartAcorCentY2$Classe)
head(base)
acoresCentY2@data$X = as.factor(base$Classe)
spplot(acoresCentY2, "X", col.regions = c("black","#66CD00"), main = "Açores (Central)",
       colorkey = list(space="left", height = 0.4), col = NA)

#********************************************************************************#
# Cartograma - Açores/ORIENTAL ####

#Leitura do ficheiro shp de dados geoespaciais dos Açores/Oriental
acoresOriY2 = readOGR(file.choose())
plot(acoresOriY2)
View(acoresOriY2)
str(acoresOriY2@data)

#Conversão dos nomes dos municípios dos Açores/Oriental (sem caracteres estranhos)
acoresOriY2@data$Concelho = iconv(acoresOriY2@data$Concelho, from = "UTF-8", to = "ASCII//TRANSLIT")

#Construir o cartograma
cartAcorOriY2 = data.frame(acoresOriY2@data$Concelho)
cartAcorOriY2$classe = 0
colnames(cartAcorOriY2) = c("Municípios","Classe")
str(cartAcorOriY2)

#Algoritmo que atribui à coluna 'classe' do cartograma os outliers com 1's
for(j in outY2Final)
{
  for(i in 1:70)
  {
    if (as.logical(as.character(cartAcorOriY2[i,1])==as.character(Y2names[j,1])))
    {
      cartAcorOriY2[i,2]=1
    }
  }
}

#Construir o mapa dos Açores/Oriental
base = data.frame(Municipios = acoresOriY2@data$Concelho, Classe = cartAcorOriY2$Classe)
head(base)
acoresOriY2@data$X = as.factor(base$Classe)
spplot(acoresOriY2, "X", col.regions = c("black","#66CD00"), main = "Açores (Oriental)",
       colorkey = list(space="left", height = 0.4), col = NA)

###################******************************************###################
###################     Matriz Situação Profissinal 2005     ###################
###################******************************************###################

# Leitura da matriz ####
Y3inicial = readWorksheet(ficheiro2005, sheet="MatrizSituaçãoTrabI2005")
Y3 = Y3inicial
GrandTotal = Y3[-309, 5] #Número de habitantes internos que entraram, até à data de 2011, em cada município
Y3 = Y3[-309,-c(1,5)]
colnames(Y3) = c("Desempregado","Empregado","Inativo")

dadosY3 = Y3/GrandTotal #Distribuição da atividade económica dos habitantes internos que entraram em cada município
colnames(dadosY3) = colnames(Y3)

#Variável que guarda só os nomes dos municípios
Y3names = Y3inicial[-309,-5] #Adicionar a coluna dos nomes e retirar a coluna do 'GrandTotal'
colnames(Y3names) = c("Municípios","Desempregado","Empregado","Inativo")

#------------------------------------------------------------------------------#
# Metodologias NUMÉRICAS ####

### PBS segundo Filzmoser ####
dadosY3ilr = pivotCoord(dadosY3) #Coordenadas ilr-transformadas 
names(dadosY3ilr) = paste0("z", 1:ncol(dadosY3ilr))
Base = orthbasis(ncol(dadosY3)) #Base ortonormal
matrizContr = Base$V #Matriz de contrastes (colunas ortonormais)
PBS = Base$basisv #Tabela da PBS

#********************************************************************************#

### Estimador MCD ####
robustezY3 = covMcd(dadosY3ilr, alpha = 0.75, nsamp = "deterministic") #Algoritmo Fast-MCD
#NOTA: obtém-se as estimativas robustas da média e da matriz de covariância
summary(robustezY3)

dY3 = determinant(robustezY3$cov,logarithm = TRUE)$modulus[1] #-5.518727
ddY3 = determinant(robustezY3$cov,logarithm = FALSE)$modulus[1] #0.004010952
#NOTA: ln(0.004010952)=-5.518727

#Outliers de 'covMcd'
outCovY3 = which(robustezY3$mcd.wt==0) #11 outliers
Y3names[outCovY3,1] #Nomes dos municípios outliers

#Gráfico de 'covMcd'
plot(robustezY3, which = "distance", classic = TRUE)

#********************************************************************************#

### Distância de Mahalanobis ####
#Adjusted Quantile Plot
#Função aq.plot original alterada
aq.plot = function (x, delta = qchisq(0.975, df = ncol(x)), quan = 1/2,alpha = 0.05) 
{
  if (is.vector(x) == TRUE || ncol(x) == 1) {
    stop("x must be at least two-dimensional")
  }
  covr <- covMcd(x, alpha = quan)
  mcd <- determinant(covr$cov, logarithm = TRUE)$modulus[1] #Acrescentei esta linha
  dist <- mahalanobis(x, center = covr$center, cov = covr$cov)
  s <- sort(dist, index = TRUE)
  z <- x
  if (ncol(x) > 2) {
    p <- princomp(x, covmat = covr)
    z <- p$scores[, 1:2]
    sdprop <- (p$sd[1] + p$sd[2])/sum(p$sd)
    cat("Projection to the first and second robust principal components.\\n")
    cat("Proportion of total variation (explained variance): ")
    cat(sdprop)
    cat("\\n")
  }
  par(mfrow = c(2, 2), mai = c(0.8, 0.6, 0.2, 0.2), mgp = c(2.4,1, 0))
  plot(z, col = 3, type = "n", xlab = "", ylab = "")
  text(z, dimnames(as.data.frame(z))[[1]], col = 3, cex = 0.8)
  plot(s$x, (1:length(dist))/length(dist), col = 3, xlab = "Ordered squared robust distance", 
       ylab = "Cumulative probability", type = "n")
  text(s$x, (1:length(dist))/length(dist), as.character(s$ix),col = 3, cex = 0.8)
  t <- seq(0, max(dist), by = 0.01)
  lines(t, pchisq(t, df = ncol(x)), col = 6)
  abline(v = delta, col = 5)
  text(x = delta, y = 0.4, paste(100 * (pchisq(delta, df = ncol(x))), 
                                 "% Quantile", sep = ""), col = 5, pos = 2, srt = 90,cex = 0.8)
  xarw <- arw(x, covr$center, covr$cov, alpha = alpha)
  if (xarw$cn < Inf) {
    abline(v = xarw$cn, col = 4)
    text(x = xarw$cn, y = 0.4, "Adjusted Quantile", col = 4,pos = 4, srt = 90, cex = 0.8)
  }
  plot(z, col = 3, type = "n",
       main = paste("Outliers based on ",100 * (pchisq(delta, df = ncol(x))), "% quantile", sep = ""), 
       xlab = "", ylab = "")
  if (any(dist > delta)) {
    text(z[dist > delta, 1], z[dist > delta, 2],
         dimnames(as.data.frame(x)[dist > delta, ])[[1]], col = 2, cex = 0.8)
  }
  if (any(dist <= delta)) {
    text(z[dist <= delta, 1], z[dist <= delta, 2],
         dimnames(as.data.frame(x)[dist <= delta, ])[[1]], col = 3, cex = 0.8)
  }
  plot(z, col = 3, type = "n", main = "Outliers based on adjusted quantile", 
       xlab = "", ylab = "")
  if (xarw$cn < Inf) {
    text(z[dist > xarw$cn, 1], z[dist > xarw$cn, 2],
         dimnames(as.data.frame(x)[dist > xarw$cn, ])[[1]], col = 2, cex = 0.8)
  }
  text(z[dist <= xarw$cn, 1], z[dist <= xarw$cn, 2],
       dimnames(as.data.frame(x)[dist <= xarw$cn, ])[[1]], col = 3, cex = 0.8)
  
  o <- (sqrt(dist) > max(sqrt(xarw$cn), sqrt(qchisq(0.975, dim(x)[2]))))
  aqplot <- list(outliers = o, MCD1 = mcd) #Acrescentei o último parâmetro
  return(aqplot)
}

#Usar a função aqplot alterada
adjquanY3 = aq.plot(dadosY3ilr, delta = qchisq(0.975, df = ncol(dadosY3ilr)),
                    quan = 0.75, alpha = 0.975)

#Valor do estimador MCD
adjquanY3$MCD1

#Outliers
outAQplotY3 = which(adjquanY3$outliers == TRUE) #Obtém-se as observações que são outliers
length(outAQplotY3)
Y3names[outAQplotY3,1] #Nomes dos municípios outliers

#SOLUÇÃO: guardar os outliers associados ao menor valor de MCD
out1Y3Final = c(83,103,116,166,177,188,213,226,267,294,301)
# 11 outliers com MCD=-5.518727

#Nomes dos municípios outliers para a solução
Y3names[c(83,103,116,166,177,188,213,226,267,294,301),1]

#................................................................................#

#Distance-Distance Plot
#Função ddplot original alterada
dd.plot = function (x, quan = 1/2, alpha = 0.025, ...) 
{
  if (!is.matrix(x) && !is.data.frame(x)) 
    stop("x must be matrix or data.frame")
  rob <- covMcd(x, alpha = quan)
  mcd <- determinant(rob$cov, logarithm = TRUE)$modulus[1] #Acrescentei esta linha
  xarw <- arw(x, rob$center, rob$cov, alpha = alpha)
  distcla <- sqrt(mahalanobis(x, center = apply(x, 2, mean), 
                              cov = cov(x)))
  distrob <- sqrt(mahalanobis(x, center = rob$center, cov = rob$cov))
  plot(distcla, distrob, main = "Distance-Distance Plot", xlab = "Mahalanobis Distance", 
       ylab = "Robust Distance", type = "n", ...)
  if (xarw$cn != Inf) {
    alpha <- sqrt(c(xarw$cn, qchisq(c(0.75, 0.5, 0.25), ncol(x))))
  }
  else {
    alpha <- sqrt(qchisq(c(0.975, 0.75, 0.5, 0.25), ncol(x)))
  }
  abline(h = alpha[1])
  abline(v = alpha[1])
  abline(a = 0, b = 1)
  lpch <- c(3, 3, 16, 1, 1)
  lcex <- c(1.5, 1, 0.5, 1, 1.5)
  lalpha <- length(alpha)
  xs <- scale(x) - min(scale(x))
  eucl <- sqrt(apply(xs^2, 1, sum))
  rbcol <- rev(rainbow(nrow(x), start = 0,
                       end = 0.7))[as.integer(cut(eucl, nrow(x), labels = 1:nrow(x)))]
  rd <- distrob
  for (j in 1:lalpha) {
    if (j == 1) {
      points(distcla[rd >= alpha[j]], distrob[rd >= alpha[j]], 
             pch = lpch[j], cex = lcex[j], col = rbcol[rd >= alpha[j]])
    }
    if (j > 1 & j < lalpha) 
      points(distcla[rd < alpha[j - 1] & rd >= alpha[j]], 
             distrob[rd < alpha[j - 1] & rd >= alpha[j]], 
             cex = lcex[j], pch = lpch[j], col = rbcol[rd < alpha[j - 1] & rd >= alpha[j]])
    if (j == lalpha) {
      points(distcla[rd < alpha[j - 1] & rd >= alpha[j]], 
             distrob[rd < alpha[j - 1] & rd >= alpha[j]], 
             cex = lcex[j], pch = lpch[j], col = rbcol[rd < alpha[j - 1] & rd >= alpha[j]])
      points(distcla[rd < alpha[j]], distrob[rd < alpha[j]], 
             pch = lpch[j + 1], cex = lcex[j + 1], col = rbcol[rd < alpha[j]])
    }
  }
  o <- (rd > min(sqrt(xarw$cn), sqrt(qchisq(0.975, dim(x)[2]))))
  ddplot <- list(outliers = o, md.cla = distcla, md.rob = distrob, MCD2 = mcd) #Acrescentei o último parâmetro
  return(ddplot)
}

#Usar a função ddplot alterada
DDplotY3 = dd.plot(dadosY3ilr, quan = 0.75, alpha = 0.975)

#Valor do estimador MCD
DDplotY3$MCD2

#Outliers
outDDplotY3 = which(DDplotY3$outliers == TRUE)
length(outDDplotY3)
Y3names[outDDplotY3,1]

#SOLUÇÃO: guardar os outliers associados ao menor valor de MCD
out2Y3Final = c(83,103,116,166,177,188,195,213,222,226,267,294,297,301)
# 14 outliers com MCD=-5.518727

#Nomes dos municípios outliers para a solução
Y3names[c(83,103,116,166,177,188,195,213,222,226,267,294,297,301),1]

#................................................................................#

#Outlier CoDa
#Função outCoDa original alterada
outCoDa = function (x, quantile = 0.975, method = "robust", h = 1/2, coda = TRUE) 
{
  if (dim(x)[2] < 2) 
    stop("need data with at least 2 variables")
  covEst <- function(x, type) {
    standard <- function(x) {
      list(mean = colMeans(x, na.rm = TRUE), varmat = cov(x))
    }
    robust <- function(x) {
      v <- robustbase::covMcd(x)
      mcd <- base::determinant(v$cov, logarithm = TRUE)$modulus[1] #Acrescentei esta linha
      list(mean = v$center, varmat = v$cov, MCD = mcd) #Acrescentei o último parâmetro
    }
    switch(type, standard = standard(x), robust = robust(x))
  }
  if (!is.logical(coda) & !is.function(coda)) {
    stop("coda must be logical or function")
  }
  if (!is.logical(coda)) {
    x <- coda(x)
  }
  else if (coda) {
    x <- pivotCoord(x)
  }
  cv <- covEst(x, "robust")
  cvc <- covEst(x, "standard")
  dM <- sqrt(mahalanobis(x, center = cv$mean, cov = cv$varmat))
  dMc <- sqrt(mahalanobis(x, center = cvc$mean, cov = cvc$varmat))
  limit <- sqrt(qchisq(p = quantile, df = ncol(x) - 1))
  res <- list(mahalDist = dM, limit = limit, outlierIndex = dM > limit, method = method,
              om2 = dMc > limit, m2 = dMc, coda = coda, MCD3 = cv$MCD) #Acrescentei o último parâmetro
  class(res) <- "outCoDa"
  invisible(res)
}

#Usar a função outCoDa alterada
outRobY3 = outCoDa(dadosY3, quantile = 0.975, method = "robust", h = 231)
str(outRobY3)

#Valor do estimador MCD
outRobY3$MCD3

#Outliers
outROBY3 = which(outRobY3$outlierIndex == TRUE)
Y3names[outROBY3,1]

#SOLUÇÃO: guardar os outliers associados ao menor valor de MCD
out3Y3Final = c(34,83,103,116,145,166,177,188,195,196,213,222,226,242,267,284,294,297,301)
# 19 outliers com MCD=-5.539568

#Nomes dos municípios outliers para a solução
Y3names[c(34,83,103,116,145,166,177,188,195,196,213,222,226,242,267,284,294,297,301),1]

#------------------------------------------------------------------------------#
# Metodologias GRÁFICAS ####
#Função mvoutlier.CoDa original alterada
mvoutlier.CoDa = function (x, quan = 0.75, alpha = 0.025,
                           col.quantile = c(0,0.05, 0.1, 0.5, 0.9, 0.95, 1),
                           symb.pch = c(3, 3, 16, 1,1), symb.cex = c(1.5, 1, 0.5, 1, 1.5),
                           adaptive = TRUE)
{
  if (!is.matrix(x) && !is.data.frame(x)) 
    stop("x must be matrix or data.frame")
  if (ncol(x) < 3) 
    stop("x must have at least 3 compositional parts")
  Z <- pivotCoord(x)
  V <- (orthbasis(ncol(x))$V)
  Zj <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  for (j in 1:ncol(x)) {
    Zj[, j] <- pivotCoord(cbind(x[, j], x[, -j]))[, 1]
  }
  dimnames(Zj)[[2]] <- names(x)
  rob <- covMcd(Z, alpha = quan)
  mcd <- determinant(rob$cov, logarithm = TRUE)$modulus[1] #Acrescentei esta linha
  if (adaptive) {
    Zarw <- arw(Z, rob$center, rob$cov, alpha = alpha)
    if (Zarw$cn != Inf) {
      alpha1 <- sqrt(c(Zarw$cn, qchisq(c(0.75, 0.5, 0.25), ncol(Z))))
    }
    else {
      alpha1 <- sqrt(qchisq(c(0.975, 0.75, 0.5, 0.25), ncol(Z)))
    }
    rd2 <- mahalanobis(Z, center = Zarw$m, cov = Zarw$c)
  }
  else {
    cutoff <- qchisq(1 - alpha, ncol(Z))
    rd2 <- mahalanobis(Z, center = rob$center, cov = rob$cov)
    Zarw <- list(m = rob$center, c = rob$cov, cn = cutoff, w = as.numeric(rd2 < cutoff))
    alpha1 <- sqrt(qchisq(c(0.975, 0.75, 0.5, 0.25), ncol(Z)))
  }
  rd <- sqrt(rd2)
  covobj <- list(center = Zarw$m, cov = Zarw$c, n.obs = length(rd), mah = rd)
  Z.pca <- suppressWarnings(princomp(Z, covmat = covobj, cor = FALSE))
  pcaclr <- Z.pca
  eval <- eigen(Zarw$c)$values
  pcaclr$sdev <- sqrt(eval)
  pcaclr$scores <- Z.pca$scores
  pcaclr$loadings <- V %*% Z.pca$loadings
  dimnames(pcaclr$loadings)[[1]] <- names(x)
  pcaobj <- list(method = "robust", eigenvalues = eval, princompOutputClr = pcaclr)
  class(pcaobj) <- "pcaCoDa"
  Zcent <- scale(Zj, center = apply(Zj, 2, median), scale = FALSE)
  eucl <- apply(abs(Zcent), 1, median)
  out <- (!Zarw$w)
  lq <- length(col.quantile)
  colcol <- rev(rainbow(lq - 1, start = 0,
                        end = 0.7))[as.integer(cut(eucl, quantile(eucl, col.quantile, labels = 1:(lq - 1))))]
  colbw <- rev(gray(seq(from = 0.1, to = 0.9,
                        length = lq - 1)))[as.integer(cut(eucl, quantile(eucl, col.quantile, labels = 1:(lq - 1))))]
  pchvec <- rep(symb.pch[1], nrow(Zj))
  cexvec <- rep(symb.cex[1], nrow(Zj))
  if (length(symb.pch) == 5 & length(symb.cex) == 5) {
    lalpha <- length(alpha1)
    for (j in 1:(lalpha)) {
      pchvec[rd < alpha1[j]] <- symb.pch[j + 1]
      cexvec[rd < alpha1[j]] <- symb.cex[j + 1]
    }
  }
  mvoutlierCoDa <- list(ilrvariables = Zj, outliers = out, pcaobj = pcaobj, colcol = colcol,
                        colbw = colbw, pchvec = pchvec, cexvec = cexvec,
                        cov = rob$cov, MCD4 = mcd) #Acrescentei os dois últimos parâmetros
  class(mvoutlierCoDa) <- "mvoutlierCoDa"
  return(mvoutlierCoDa)
}

#Usar a função mvoutlierCoDa alterada
resY3 = mvoutlier.CoDa(dadosY3, quan = 0.75, alpha = 0.975, adaptive = TRUE)
str(resY3)

#Valor do estimador MCD
resY3$MCD4

#Outliers
outResY3 = which(resY3$outliers == TRUE)
length(outResY3)
Y2names[outResY3,1]

#SOLUÇÃO: guardar os outliers associados ao menor valor de MCD
out4Y3Final = c(83,103,116,166,177,188,195,196,213,222,226,267,294,297,301)
# 15 outliers com MCD=-5.518727

#Nomes dos municípios outliers para a solução
Y3names[c(83,103,116,166,177,188,195,196,213,222,226,267,294,297,301),1]

#................................................................................#

#Biplot: para visualização dos outliers
plot(resY3, which = "biplot", onlyout = TRUE, symb = TRUE, symbtxt = TRUE)
#Apenas os outliers são exibidos e identificados por números
#NOTA: os outliers estão ordenados

#Biplot: para visualização dos ângulos entre os links
biplotAngY3 = pcaCoDa(dadosY3, method = "robust")
summary(biplotAngY3)
biplot(biplotAngY3, xlabs = rep(".", nrow(dadosY3)))

#................................................................................#

#Gráfico de dispersão univariado
plot(resY3, which = "uni", onlyout = TRUE, symb = TRUE, symbtxt = TRUE)
#Apenas os outliers são exibidos e identificados por números

#********************************************************************************#

#Diagrama Ternário
ternaryDiag(dadosY3, colnames(dadosY3), text = NULL, grid = FALSE,
            gridCol = grey(0.5), mcex = 1.2, line = "ellipse", robust = TRUE,
            group = NULL, tol = 0.975, col = "#458B74")

#Exemplo 3.2.1
X = matrix(c(0.104,0.510,0.386, 0.068,0.546,0.386, 0.094,0.561,0.345,
             0.072,0.571,0.357, 0.094,0.502,0.404), nrow = 5, ncol = 3, byrow = TRUE)
colnames(X) = c("Desempregado","Empregado","Inativo")
ternaryDiag(X, colnames(X), text = NULL, grid = FALSE, gridCol = grey(0.5),
            mcex = 1.2, robust = TRUE, group = NULL, col = "black")

Y = X%*%diag(c(0.723,0.114,0.163)); Y
clo(Y) #Operação de fecho para a matriz Y
ternaryDiag(clo(Y), colnames(X), text = NULL, grid = FALSE, gridCol = grey(0.5),
            mcex = 1.2, robust = TRUE, group = NULL, col = "black")

#********************************************************************************#

#Matriz de variação (Tabela 4.14)
dadosY3TabVar = dadosY3[,c(2,1,3)] #Reordenar as colunas da matriz
varY3 = matrix(0,3,3)
varY3[1,2] = var(log(dadosY3TabVar[,1]/dadosY3TabVar[,2]))
varY3[1,3] = var(log(dadosY3TabVar[,1]/dadosY3TabVar[,3]))
varY3[2,3] = var(log(dadosY3TabVar[,2]/dadosY3TabVar[,3]))
varY3 #Imprimir a matriz

#********************************************************************************#

#Matriz de correlação (Tabela 4.15)
corrY3 = matrix(0,3,3)
corrY3[1,2] = cor(log(dadosY3[,2]/dadosY3[,1]), log(dadosY3[,3]/dadosY3[,1]))
corrY3[1,3] = cor(log(dadosY3[,2]/dadosY3[,1]), log(dadosY3[,2]/dadosY3[,3]))
corrY3[2,3] = cor(log(dadosY3[,3]/dadosY3[,1]), log(dadosY3[,2]/dadosY3[,3]))
corrY3 #Imprimir a matriz

#------------------------------------------------------------------------------#
# Cartogramas ####

#Grupo de outliers final
outY3Final = c(83,103,116,166,177,188,213,226,267,294,301)

#Conversão dos nomes dos municípios da matriz (comum a todos os cartogramas)
Y3names[,1] = stri_trans_general(Y3names[,1], "Latin-ASCII")
Y3names[,1] = as.factor(toupper(Y3names[,1])) #Colocar em letras maiúsculas

#********************************************************************************#
# Portugal ####

#Leitura do ficheiro shp de dados geoespaciais de Portugal
portugalY3 = readOGR(file.choose())
plot(portugalY3)
View(portugalY3)
str(portugalY3@data)

#Conversão dos nomes dos municípios de Portugal (sem caracteres estranhos)
portugalY3@data$Concelho = iconv(portugalY3@data$Concelho, from = "UTF-8", to = "ASCII//TRANSLIT")

#Construir o cartograma
cartPortY3 = data.frame(portugalY3@data$Concelho)
cartPortY3$classe = 0
colnames(cartPortY3) = c("Municípios","Classe")
str(cartPortY3)

#Algoritmo que atribui à coluna 'classe' do cartograma os outliers com 1's
for(j in outY3Final)
{
  for(i in 1:3223)
  {
    if (as.logical(as.character(cartPortY3[i,1])==as.character(Y3names[j,1])))
    {
      cartPortY3[i,2]=1
    }
  }
}

#Construir o mapa de Portugal
base = data.frame(Municipios = portugalY3@data$Concelho, Classe = cartPortY3$Classe)
head(base)
portugalY3@data$X = as.factor(base$Classe)
spplot(portugalY3, "X", col.regions = c("black","#66CD00"), main = "Portugal (Continental)",
       colorkey = list(space = "left", height = 0.4), col = NA)

#********************************************************************************#
# Madeira ####

#Leitura do ficheiro shp de dados geoespaciais da Madeira
madeiraY3 = readOGR(file.choose())
plot(madeiraY3)
View(madeiraY3)
str(madeiraY3@data)

#Conversão dos nomes dos municípios da Madeira (sem caracteres estranhos)
madeiraY3@data$Concelho = iconv(madeiraY3@data$CONCELHO, from = "UTF-8", to = "ASCII//TRANSLIT")

#Construir o cartograma
cartMadY3= data.frame(madeiraY3@data$CONCELHO)
cartMadY3$classe = 0
colnames(cartMadY3) = c("Municípios","Classe")
str(cartMadY3)

#Algoritmo que atribui à coluna 'classe' do cartograma os outliers com 1's
for(j in outY3Final)
{
  for(i in 1:538)
  {
    if (as.logical(as.character(cartMadY3[i,1])==as.character(Y3names[j,1])))
    {
      cartMadY3[i,2]=1
    }
  }
}

#Construir o mapa da Madeira
base = data.frame(Municipios = madeiraY3@data$CONCELHO, Classe = cartMadY3$Classe)
head(base)
madeiraY3@data$X = as.factor(base$Classe)
spplot(madeiraY3, "X", col.regions = c("black","#66CD00"), main = "Madeira",
       colorkey = list(space="left", height = 0.4), col = NA, xlim=c(287859.9,389570.9),
       ylim=c(3580000,3665000))

#********************************************************************************#
# Açores/OCIDENTAL ####

#Leitura do ficheiro shp de dados geoespaciais dos Açores/Ocidental
acoresOciY3 = readOGR(file.choose())
plot(acoresOciY3)
View(acoresOciY3)
str(acoresOciY3@data)

#Conversão dos nomes dos municípios dos Açores/Ocidental (sem caracteres estranhos)
acoresOciY3@data$Concelho = iconv(acoresOciY3@data$Concelho, from = "UTF-8", to = "ASCII//TRANSLIT")

#Construir o cartograma
cartAcorOciY3 = data.frame(acoresOciY3@data$Concelho)
cartAcorOciY3$classe = 0
colnames(cartAcorOciY3) = c("Municípios","Classe")
str(cartAcorOciY3)

#Algoritmo que atribui à coluna 'classe' do cartograma os outliers com 1's
for(j in outY3Final)
{
  for(i in 1:12)
  {
    if (as.logical(as.character(cartAcorOciY3[i,1])==as.character(Y3names[j,1])))
    {
      cartAcorOciY3[i,2]=1
    }
  }
}

#Construir o mapa dos Açores/Ocidental
base = data.frame(Municipios = acoresOciY3@data$Concelho, Classe = cartAcorOciY3$Classe)
head(base)
acoresOciY3@data$X = as.factor(base$Classe)
spplot(acoresOciY3, "X", col.regions = c("black","#66CD00"), main = "Açores (Ocidental)",
       colorkey = list(space="left", height = 0.4), col = NA)

#********************************************************************************#
# Açores/CENTRAL ####

#Leitura do ficheiro shp de dados geoespaciais dos Açores/Central
acoresCentY3 = readOGR(file.choose())
plot(acoresCentY3)
View(acoresCentY3)
str(acoresCentY3@data)

#Conversão dos nomes dos municípios dos Açores/Central (sem caracteres estranhos)
acoresCentY3@data$Concelho = iconv(acoresCentY3@data$Concelho, from = "UTF-8", to = "ASCII//TRANSLIT")

#Construir o cartograma
cartAcorCentY3 = data.frame(acoresCentY3@data$Concelho)
cartAcorCentY3$classe = 0
colnames(cartAcorCentY3) = c("Municípios","Classe")
str(cartAcorCentY3)

#Algoritmo que atribui à coluna 'classe' do cartograma os outliers com 1's
for(j in outY3Final)
{
  for(i in 1:75)
  {
    if (as.logical(as.character(cartAcorCentY3[i,1])==as.character(Y3names[j,1])))
    {
      cartAcorCentY3[i,2]=1
    }
  }
}

#Construir o mapa dos Açores/Central
base = data.frame(Municipios = acoresCentY3@data$Concelho, Classe = cartAcorCentY3$Classe)
head(base)
acoresCentY3@data$X = as.factor(base$Classe)

#NOTA: este 'if' é importante para bases de dados onde não haja outliers
if(length(unique(acoresCentY3@data$X))==1)
{
  spplot(acoresCentY3, "X", col.regions = c("black"), main = "Açores (Central)",
         colorkey = list(space="left", height = 0.4), col = NA)
} else
{
  spplot(acoresCentY3, "X", col.regions = c("black","#66CD00"), main = "Açores (Central)",
         colorkey = list(space="left", height = 0.4), col = NA)
}

#********************************************************************************#
# Açores/ORIENTAL ####

#Leitura do ficheiro shp de dados geoespaciais dos Açores/Oriental
acoresOriY3 = readOGR(file.choose())
plot(acoresOriY3)
View(acoresOriY3)
str(acoresOriY3@data)

#Conversão dos nomes dos municípios dos Açores/Oriental (sem caracteres estranhos)
acoresOriY3@data$Concelho = iconv(acoresOriY3@data$Concelho, from = "UTF-8", to = "ASCII//TRANSLIT")

#Construir o cartograma
cartAcorOriY3 = data.frame(acoresOriY3@data$Concelho)
cartAcorOriY3$classe = 0
colnames(cartAcorOriY3) = c("Municípios","Classe")
str(cartAcorOriY3)

#Algoritmo que atribui à coluna 'classe' do cartograma os outliers com 1's
for(j in outY3Final)
{
  for(i in 1:70)
  {
    if (as.logical(as.character(cartAcorOriY3[i,1])==as.character(Y3names[j,1])))
    {
      cartAcorOriY3[i,2]=1
    }
  }
}

#Construir o mapa dos Açores/Oriental
base = data.frame(Municipios = acoresOriY3@data$Concelho, Classe = cartAcorOriY3$Classe)
head(base)
acoresOriY3@data$X = as.factor(base$Classe)

#NOTA: este 'if' é importante para bases de dados onde não haja outliers
if(length(unique(acoresOriY3@data$X))==1)
{
  spplot(acoresOriY3, "X", col.regions = c("black"), main = "Açores (Oriental)",
         colorkey = list(space="left", height = 0.4), col = NA)
} else
{
  spplot(acoresOriY3, "X", col.regions = c("black","#66CD00"), main = "Açores (Oriental)",
         colorkey = list(space="left", height = 0.4), col = NA)
}