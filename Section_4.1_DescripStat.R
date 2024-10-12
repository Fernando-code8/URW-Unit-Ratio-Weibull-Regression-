###################################################################################
# PAPER: The COVID-19 mortality rate in Latin America: a cross-country analysis
# SUBSECTION: 4.1. Descriptive summary and correlation analysis
# GOAL: Doing a descriptive statistical analysis of the data set.
# AUTHOR: Fernando Arturo Pe\~na-Ram\'irez, Renata Rojas Guerra and 
#         Fernando Jos\'e Monteiro de Ara\'ujo
# LAST UPDATE: October 09, 2024
###################################################################################

library(readxl)
library(fBasics)
library(xtable)
library(GGally) 
library(ggplot2)
library(ggcorrplot)
library(knitr)
library(tidyverse)


dados<- read_excel("dados_final.xlsx") %>%
  mutate(GDP=GDP/1000000)
View(dados)
class(dados)
attach(dados)
str(dados)

### Descriptive Summary
stat<-basicStats(dados[,c(4,6:13)])
stat<-as.matrix(stat)
class(stat)
stat<-t(stat)
stat<-as.data.frame(stat)
stat["CV"]<-(stat$Stdev/stat$Mean)*100
xtable(stat,digits=4)

### Correlation
round(cor(dados[,c(4,6:13)],method="spearman"),4) 
dados1<-dados[,c(4,6:13)]

### histogram
hist(dados1$MR,xlab="MR",ylab="Density", main="")

### boxplot
boxplot(dados1$MR,ylab="MR")

### scatter plot
plot(dados1$PD,dados1$MR,xlab="PD",ylab="MR")
plot(dados1$GDP,dados1$MR,xlab="GDP",ylab="MR")
plot(dados1$UP,dados1$MR,xlab="UP",ylab="MR")
plot(dados1$HB,dados1$MR,xlab="HB",ylab="MR")
plot(dados1$HDI,dados1$MR,xlab="HDI",ylab="MR")
plot(dados1$CHE,dados1$MR,xlab="CHE",ylab="MR")
plot(dados1$DP,dados1$MR,xlab="DP",ylab="MR")
plot(dados1$DGGHE,dados1$MR,xlab="DGGHE",ylab="MR")

