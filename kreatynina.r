#!/usr/bin/env Rscript

#####zmienne do edytowania#####

#parametr n w walidacji krzyzowej
#ma najwiekszy wplyw na czas wykonania skryptu
cv_n = 1

#czy liczyc dodatkowo regresje (z tuningiem) dla wszystkich atrybutow
#bardzo powolna operacja z relatywnie slabymi wynikami
useAllAttributess=FALSE

#####koniec zmiennych do (czestego) edytowania######

#dolaczamy biblioteki
library(dmr.regeval)
library(dmr.trans)
library(dmr.util)
library(dmr.linreg)
library(rpart)
library(dmr.claseval)
library(FSelector)
library(dmr.regtree)
library(dmr.kernel)
library(e1071)


source("GetWhichHighestNumber.r")
source("DrawSubsetCorrgram.r")
source("CreateDfFromAttributes.r")

#wczytujemy baze - pierwszy wiersz naglowkowy
df_base <- read.csv("BAZA_kreatynina.csv",header=TRUE)

#ustalamy nazwy funkcji celu
kreat1_name <- ("kreat_1")
kreat_names <- c("kreat_2", "kreat_3", "kreat_7", "kreat_14", "kreat_30")

#tworzymy baze bez dodatkowych wartosci wyjsciowych
reduced_base<- df_base[, -which(names(df_base) %in% kreat_names)]

#ustalamy kilka "naiwnych" licznosci atrybutow
n_lengths <- c(2:20,25,30,40,50,60,70,80,90)

#tworzy katalog na wykresy
dir.create(file.path("plots"),showWarnings=FALSE)

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

#sprawdzamy jakosc regresji walidacja krzyzowa dla calej bazy - jako punkt odniesienia dla selekcji atrybutow
db_all_cv10<-dmr.claseval::crossval(rpart,kreat_1~.,as.data.frame(impute(reduced_base,what="median")),k=10,n=cv_n)
#liczymy wspolczynnik determinacji - R2
db_all_cv10_r2<- dmr.regeval::r2(db_all_cv10$pred, db_all_cv10$true)
mse(db_all_cv10$pred, db_all_cv10$true)

#metody selekcji atrybutow
#usuwamy i tworzymy na nowo, zeby uniknac starych wartosci
if(exists("namesOfTests")){
  rm(namesOfTests)
}
namesOfTests <- c("linear.correlation","rank.correlation","random.forest")


#tworzymy dataframe z wynikami selekcji
selectionResults<-as.data.frame(setNames(replicate(length(namesOfTests),numeric(0), simplify = F), namesOfTests))
selectionResultsR2<-as.data.frame(setNames(replicate(length(namesOfTests),numeric(0), simplify = F), namesOfTests))
weights <- list(namesOfTests)

#iterujemy po metodach selekcji atrybutow
# i wykonujemy przypisanie wag
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
  # dla wybranych "naiwnych" dlugosci przeprowadzamy sprawdzenie
  # tzn. dokonujemy walidacji krzyzowej i mierzymy R2
  # zakładamy, że w ten sposób dobierzemy najlepszy i-elementowy podzbior atrybutow
  for (i in n_lengths){
    sub_set <- cutoff.k(weights[[algor]], i)
    sub_base <- reduced_base[, which(names(reduced_base) %in% c(sub_set,"kreat_1"))]
    testRes<-dmr.claseval::crossval(rpart,kreat_1~.,as.data.frame(impute(sub_base,what="median")),k=10,n=cv_n)
    selectionResults[toString(i),algor]<-dmr.regeval::r2(testRes$pred, testRes$true)
    #sprawdzmy takze dodatkowo korelacje przewidywanych wynikow z rzeczywistymi
    selectionResultsR2[toString(i),algor]<-cor(testRes$pred, testRes$true)
  }
}

#rysujemy zaleznosc R2 od mocy podzbioru atrybutow
for (algor in namesOfTests){
  png(paste("plots/","R2_",algor,".png",sep=""),width=600,height=600)
  plot(n_lengths,selectionResults[,algor],
       type="o",xlab="Moc podzbioru atrybutow",ylab="Wspolczynnik determinacji",
       main=paste("Zaleznosc wspolczynnika determinacji od mocy podzbioru atrybutow \ndla ",algor,sep=""))
  dev.off()
}

#rysujemy zaleznosc korelacji rzeczywistych wynikow z przewidywaniami od mocy podzbioru atrybutow
for (algor in namesOfTests){
  png(paste("plots/","RegCor_",algor,".png",sep=""),width=600,height=600)
  plot(n_lengths,selectionResultsR2[,algor],
       type="o",xlab="Moc podzbioru atrybutow",ylab="Korelacja prawdziwych oraz przewidywanych wartosci funkcji docelowej",
       main=paste("Zaleznosc korelacji prawdziwych oraz przewidywanych wartosci funkcji docelowej \nod mocy podzbioru atrybutow dla ",algor,sep=""))
  dev.off()
}

# wybieramy atrybuty za pomoca CFS
cfsAtts <- cfs(formula=kreat_1~.,data=reduced_base)

# wybieramy najlepsze testy (kryterium - najwieksze R2)
#z pelnej rozpatrywanej puli 
selectedAtts<-GetWhichHighestNumber(selectionResults)
selectedSubsets<-list()
for (t in namesOfTests){
  attVec<-cutoff.k(weights[[t]],selectedAtts[1,t])
  selectedSubsets[[t]]<-CreateDfFromAttributes(reduced_base,attVec,"kreat_1")
}
#dodajemy takze atrybuty z cfs
selectedSubsets[["cfs"]]<- CreateDfFromAttributes(reduced_base,cfsAtts,"kreat_1")
namesOfTests<-c(namesOfTests,"cfs")

#rysujemy korelogramy dla wybranych metod selekcji atrybutow 
for(corGraph in namesOfTests){
  png(paste("plots/","corrgram",corGraph,".png",sep=""),width=600,height=600)
  DrawSubsetCorrgram(selectedSubsets[[corGraph]],corGraph)
  dev.off()
}

#dodajemy pelna baze jako punkt odniesienia
if(useAllAttributes==TRUE){
  selectedSubsets[["all.attributes"]] <-reduced_base 
  namesOfTests<-c(namesOfTests,"all.attributes")
}
#tworzymy struktury, potrzebne do testu regresji
regressionResults<-list()
regressionR2<-as.data.frame(setNames(replicate(length(namesOfTests),numeric(0), simplify = F), namesOfTests))
bestR2<- -5

#przeprowadzamy tuning modeli (wraz z ich ocena)
for(selection in namesOfTests){
  #regression tree
  tmpResults<-tune(rpart,kreat_1~.,data=as.data.frame(impute(selectedSubsets[[selection]]),what="median"),ranges = list(minsplit = c(5,10,20,30,40,50,60,80),cp=c(0.01,0.05,0.08,0.1,0.15)),
                   tunecontrol=tune.control(nrepeat=cv_n, repeat.aggregate=mean,sampling="cross", cross=10))
  tmpR2 <- (1 -( length(selectedSubsets[[selection]]$kreat_1)*tmpResults$best.performance/((length(selectedSubsets[[selection]]$kreat_1)-1)*var(selectedSubsets[[selection]]$kreat_1))))
  print(paste(selection," regression tree: ",tmpR2))
  regressionResults[[selection]][["regression tree"]] <- tmpResults
  if(bestR2<tmpR2){
    bestR2<-tmpR2
    bestModel <- tmpResults
  }
  regressionR2["regression tree",selection] <- tmpR2

  #linear regression
  tmpResults<-tune(lm,kreat_1~.,data=as.data.frame(impute(selectedSubsets[[selection]]),what="median"),
                   tunecontrol=tune.control(nrepeat=cv_n, repeat.aggregate=mean,sampling="cross", cross=10))
  
  tmpR2 <- (1 -( length(selectedSubsets[[selection]]$kreat_1)*tmpResults$best.performance/((length(selectedSubsets[[selection]]$kreat_1)-1)*var(selectedSubsets[[selection]]$kreat_1))))
  print(paste(selection," linear regression: ",tmpR2))
  regressionResults[[selection]][["linear regression"]] <- tmpResults
  if(bestR2<tmpR2){
    bestR2<-tmpR2
    bestModel <- tmpResults
  }
  regressionR2["linear regression",selection] <- tmpR2
  
  #support vector regression - radial kernel
  tmpResults<-tune(svm,kreat_1~., data=as.data.frame(impute(selectedSubsets[[selection]]),what="median"), ranges = list(gamma = 2^(-5:1), cost = 2^(-2:1),kernel=c("polynomial","radial")),
                   tunecontrol=tune.control(nrepeat=cv_n, repeat.aggregate=mean,sampling="cross", cross=10))
  
  tmpR2 <- (1 -( length(selectedSubsets[[selection]]$kreat_1)*tmpResults$best.performance/((length(selectedSubsets[[selection]]$kreat_1)-1)*var(selectedSubsets[[selection]]$kreat_1))))
  print(paste(selection," svr: ",tmpR2))
  regressionResults[[selection]][["svr"]] <- tmpResults
  if(bestR2<tmpR2){
    bestR2<-tmpR2
    bestModel <- tmpResults
  }
  regressionR2["svr",selection] <- tmpR2
  
}

#zapis wynikow do pliku, do pozniejszej analizy
dir.create(file.path("results"),showWarnings=FALSE)
resultsFile <- file("results/regressionResults.txt", open = "wt")
sink(resultsFile)
print(regressionResults)

close(resultsFile)
#powrot do zwyklego ujscia
sink(type = "message")
sink()

#wyniki R2
write.csv(regressionR2,"results/regressionR2.csv")

bestRegressionResultIndex<- which(regressionR2== max(regressionR2,na.rm=TRUE), arr.ind = TRUE)

#wprowadzamy wektor nazw modeli
modelNames<- c("regression tree","linear regression","svr")

bestRegressionDataFrame <- as.data.frame(selectedSubsets[[namesOfTests[bestRegressionResultIndex[2]]]])
bestRegressionModel <- modelNames[bestRegressionResultIndex[1]]
bestParams<- regressionResults[namesOfTests[bestRegressionResultIndex[2]]]

testModels <-vector()
if(bestRegressionModel=="svr"){
  model1 <- svm(kreat_1~.,data=bestRegressionDataFrame,gamma=bestModel$best.parameters$gamma,cost=bestModel$best.parameters$cost,kernel=bestModel$best.parameters$kernel)
  testModel1<-tune(svm,kreat_1~.,data=bestRegressionDataFrame,gamma=bestModel$best.parameters$gamma,cost=bestModel$best.parameters$cost,kernel=bestModel$best.parameters$kernel,
                   tunecontrol=tune.control(nrepeat=cv_n, repeat.aggregate=mean,sampling="cross", cross=10))
  testModels<- c(testModel1$best.performance)
 
  for(kreat in kreat_names){
    
    kreatBase <- bestRegressionDataFrame[, -which(names(bestRegressionDataFrame) %in% kreat1_name)]
    kreatBase <- CreateDfFromAttributes(fullDF=df_base,atts=colnames(kreatBase),additionalAtt=kreat)
    testModel<- tune(svm,as.formula(paste(kreat,"~.")),data=kreatBase,gamma=bestModel$best.parameters$gamma,cost=bestModel$best.parameters$cost,kernel=bestModel$best.parameters$kernel,
                   tunecontrol=tune.control(nrepeat=cv_n, repeat.aggregate=mean,sampling="cross", cross=10))
    testModels<-c(testModels,testModel$best.performance)
  }
}
# else if(bestRegressionModel=="regression tree"){
#   model1<- rpart(kreat_1~.,data=bestRegressionDataFrame,minsplit=kreat_1~.,bestModel$best.parameters$minsplit,cp=bestModel$best.parameters$cp)
# }


