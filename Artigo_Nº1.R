#---
#title: " o meu trabalho "
#author: "Letícia"
#date: "28/maio/2019"
#output: pdf_document
#keep_tex: true
#---

#Artigo 1: Exploratory compositional data analysis using the R-package robCompositions

#help(package="robCompositions")

library(robCompositions)
dados=expenditures #Despesas
str(dados)

#Deteta quantas observações são outliers
outlierRob = outCoDa(dados, quantile=0.975, method="robust")
outlierRob

#outlierCla = outCoDa(dados, quantile=0.975, method="standard")
#outlierCla

#Resultados da distância de Mahalanobis para cada observação
outlierRob$mahalDist

#Indica quais as observações que são outliers (True ou False)
outlierRob$outlierIndex

#Gráfico onde se visualiza os outliers
plot(outlierRob)

#Análise de Componentes Principais Robusta
PrintComRob = pcaCoDa(dados, method = "robust")
plot(PrintComRob) #Não dá um biplot!

#dadosx = data.frame(dados)
#dadosx
#pcaCoDa(dadosx, method = "robust")
#biplot(pcaCoDa, choices = 1:2, scale = 1,pch=15, cex=0.8,cex.axis=0.8,arrow.len = 0.05,
#       xlab=paste(" CP1  (", (round(100*summary(pcaCoDa)$importance[2,1],digits=1)), " % )"),
#       ylab=paste(" CP2  (", (round(100*summary(pcaCoDa)$importance[2,2],digits=1)), " % )"))
