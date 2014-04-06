#!/usr/bin/env Rscript

#dolaczamy biblioteki
library(dmr.regeval)
library(dmr.trans)
library(dmr.util)
library(dmr.linreg)
library(rpart)
library(dmr.claseval)
library(FSelector)


source("GetWhichHighestNumber.r")
source("DrawSubsetCorrgram.r")

#wczytujemy baze - pierwszy wiersz naglowkowy
df_base <- read.csv("BAZA_kreatynina.csv",header=TRUE)

#ustalamy nazwy funkcji celu
kreat1_name <- ("kreat_1")
kreat_names <- c("kreat_2", "kreat_3", "kreat_7", "kreat_14", "kreat_30")

#tworzymy baze bez dodatkowych wartosci wyjsciowych
reduced_base<- df_base[, -which(names(df_base) %in% kreat_names)]

#ustalamy kilka "naiwnych" licznosci atrybutow
n_lengths <- c(2:10,15,20,40,60)




db_all_cv10<-dmr.claseval::crossval(rpart,kreat_1~.,df_base,k=10,n=3)
dmr.regeval::r2(db_all_cv10$pred, db_all_cv10$true)
mse(db_all_cv10$pred, db_all_cv10$true)

namesOfTests <- c("linear.correlation","rank.correlation","random.forest")


#selectionResults <- as.data.frame(matrix(0,0,ncol=0))
selectionResults<-as.data.frame(setNames(replicate(length(namesOfTests),numeric(0), simplify = F), namesOfTests))
weights <- list(namesOfTests)
for (algor in namesOfTests){
  if (algor==namesOfTests[1]){
    weights[[algor]] <-linear.correlation(formula=kreat_1~.,data=reduced_base)
  }
  else if(algor==namesOfTests[2]){
    weights[[algor]] <-rank.correlation(formula=kreat_1~.,data=reduced_base)
  }
  else if(algor==namesOfTests[3]){
    weights[[algor]] <- random.forest.importance(formula=kreat_1~.,data=reduced_base,importance.type=1)
  }
  for (i in n_lengths){
    sub_set <- cutoff.k(weights[[algor]], i)
    sub_base <- reduced_base[, which(names(reduced_base) %in% c(sub_set,"kreat_1"))]
    testRes<-dmr.claseval::crossval(rpart,kreat_1~.,sub_base,k=10,n=1)
    selectionResults[toString(i),algor]<-dmr.regeval::r2(testRes$pred, testRes$true)
  }
}

for (algor in namesOfTests){
  png(paste("plots/","R2_",algor,".png",sep=""),width=600,height=600)
  plot(n_lengths,selectionResults[,algor],
       type="o",xlab="Moc podzbioru atrybutow",ylab="Wspolczynnik determinacji",
       main=paste("Zaleznosc wspolczynnika determinacji od mocy podzbioru atrybutow \ndla ",algor,sep=""))
  dev.off()
}


cfsAtts <- cfs(formula=kreat_1~.,data=reduced_base)

selectedAtts<-GetWhichHighestNumber(selectionResults)
selectedSubsets<-list()
for (t in namesOfTests){
  attVec<-cutoff.k(weights[[t]],selectedAtts[1,t])
  selectedSubsets[[t]]<-reduced_base[, which(names(reduced_base) %in% c(attVec,"kreat_1"))]
}
selectedSubsets[["cfs"]]<- reduced_base[, which(names(reduced_base) %in% c(cfsAtts,"kreat_1"))]

extendedNamesOfTests<-c(namesOfTests,"cfs")

for(corGraph in extendedNamesOfTests){
  png(paste("plots/","corrgram",corGraph,".png",sep=""),width=600,height=600)
  DrawSubsetCorrgram(selectedSubsets[[corGraph]],corGraph)
  dev.off()
}

#szkic histogramu wartosci kreat_1
png(paste("plots/","hist_kreat_1",".png",sep=""),width=600,height=600)
x <- df_base$kreat_1
h<-hist(x, breaks=15, col="red", xlab="Kreatynina w 1. dniu po operacji [mg/dl]",
        main="Histogram z krzywa rozkladu Gaussa dla kreatyniny_1",
        ylab="Czestosc wystepowania")
xfit<-seq(min(x),max(x),length=40)
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
yfit <- yfit*diff(h$mids[1:2])*length(x)
lines(xfit, yfit, col="blue", lwd=2)  
dev.off()

#szkic histogramu log wartosci kreat_1
png(paste("plots/","hist_log_kreat_1",".png",sep=""),width=600,height=600)
x <- log(df_base$kreat_1)
h<-hist(x, breaks=15, col="red", xlab="Logarytm z poziomu kreatyniny w 1. dniu po operacji [mg/dl]",
        main="Histogram z krzywa rozkladu Gaussa dla logarytmu z kreatyniny_1",
        ylab="Czestosc wystepowania")
xfit<-seq(min(x),max(x),length=40)
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
yfit <- yfit*diff(h$mids[1:2])*length(x)
lines(xfit, yfit, col="blue", lwd=2)  
dev.off()

#podejrzewamy zaleznosc kreat_1 od kreat_w - narysujmy wykres
png(paste("plots/","scatterPlotKreat_w_1",".png",sep=""),width=600,height=600)
plot(df_base$kreat_w,df_base$kreat_1,xlab="Kreatynina wyjsciowa",ylab="Kreatynina 1. dnia po operacji",
     main="Zaleznosc poziomu kreatyniny 1. dnia po operacji\nod poziomu wyjsciowego")
abline(lm(df_base$kreat_1~df_base$kreat_w), col="red")
dev.off()
