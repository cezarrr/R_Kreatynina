#!/usr/bin/env Rscript

#dolaczamy biblioteki
library(dmr.regeval)
library(dmr.trans)
library(dmr.util)
library(dmr.linreg)
library(rpart)
library(dmr.claseval)
library(FSelector)


#wczytujemy baze - pierwszy wiersz naglowkowy
df_base <- read.csv("BAZA_kreatynina.csv",header=TRUE)

#ustalamy nazwy funkcji celu
kreat1_name <- ("kreat_1")
kreat_names <- c("kreat_2", "kreat_3", "kreat_7", "kreat_14", "kreat_30")

#tworzymy baze bez dodatkowych wartosci wyjsciowych
reduced_base<- df_base[, -which(names(df_base) %in% kreat_names)]

#ustalamy kilka "naiwnych" licznosci atrybutow
n_lengths <- c(5,10,20,40)


db_all_cv10<-dmr.claseval::crossval(rpart,kreat_1~.,df_base,k=10)
dmr.regeval::r2(db_all_cv10$pred, db_all_cv10$true)
mse(db_all_cv10$pred, db_all_cv10$true)

namesOfTests <- c("linear.correlation","rank.correlation","random.forest")

#selectionResults <- as.data.frame(matrix(0,0,ncol=0))
selectionResults<-as.data.frame(setNames(replicate(length(namesOfTests),numeric(0), simplify = F), namesOfTests))

for (algor in namesOfTests){
  if (algor==namesOfTests[1]){
    weights <-linear.correlation(formula=kreat_1~.,data=reduced_base)
  }
  else if(algor==namesOfTests[2]){
    weights <-rank.correlation(formula=kreat_1~.,data=reduced_base)
  }
  else if(algor==namesOfTests[3]){
    weights <- random.forest.importance(formula=kreat_1~.,data=reduced_base,importance.type=1)
  }
  for (i in n_lengths){
    
    subset <- cutoff.k(weights, i)
    sub_base <- reduced_base[, which(names(reduced_base) %in% c(subset,"kreat_1"))]
    testRes<-dmr.claseval::crossval(rpart,kreat_1~.,sub_base,k=10,n=3)
    selectionResults[toString(i),algor]<-dmr.regeval::r2(testRes$pred, testRes$true)
  }
}