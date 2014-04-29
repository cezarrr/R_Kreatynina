#!/usr/bin/env Rscript

#uwaga - skrypt tworzy katalogi w wybranym workspace'ie

#####zmienne do edytowania#####

#parametr n w walidacji krzyzowej
#ma najwiekszy wplyw na czas wykonania skryptu
cv_n = 5

#czy liczyc dodatkowo regresje (z tuningiem) dla wszystkich atrybutow
#bardzo powolna operacja z relatywnie slabymi wynikami
useAllAttributes=FALSE

#####koniec zmiennych do (czestego) edytowania#####

#####dolaczanie bibliotek i plikow z funkcjami#####

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

#####poczatek wlasciwego skryptu#####

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
db_all_cv10<-dmr.claseval::crossval(rpart,kreat_1~.,as.data.frame(impute(reduced_base,what="mean")),k=10,n=cv_n)
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
  # zakladamy, ze w ten sposob dobierzemy najlepszy i-elementowy podzbior atrybutow
  for (i in n_lengths){
    sub_set <- cutoff.k(weights[[algor]], i)
    sub_base <- reduced_base[, which(names(reduced_base) %in% c(sub_set,"kreat_1"))]
    testRes<-dmr.claseval::crossval(rpart,kreat_1~.,as.data.frame(impute(sub_base,what="mean")),k=10,n=cv_n)
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

#####tuning#####

#przeprowadzamy tuning modeli (wraz z ich ocena)
#najdluzszy czas dzialania skryptu
#najlepszy (=o najwiekszym R2) model zostanie nastepnie poddany ocenie
#a takze poddany eksperymentowi zaaplikowania takich samych parametrow do
#przewidywania poziomow kreatyniny w kolejnych dniach

#wprowadzamy wektor nazw modeli
modelNames<- c("regression tree","linear regression","svr")

for(selection in namesOfTests){
  #regression tree
  tmpResults<-tune(rpart,kreat_1~.,data=as.data.frame(impute(selectedSubsets[[selection]]),what="mean"),ranges = list(minsplit = c(5,10,20,30,40,50,60,80),cp=c(0.01,0.05,0.08,0.1,0.15)),
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
  tmpResults<-tune(lm,kreat_1~.,data=as.data.frame(impute(selectedSubsets[[selection]]),what="mean"),
                   tunecontrol=tune.control(nrepeat=cv_n, repeat.aggregate=mean,sampling="cross", cross=10))
  
  tmpR2 <- (1 -( length(selectedSubsets[[selection]]$kreat_1)*tmpResults$best.performance/((length(selectedSubsets[[selection]]$kreat_1)-1)*var(selectedSubsets[[selection]]$kreat_1))))
  print(paste(selection," linear regression: ",tmpR2))
  regressionResults[[selection]][["linear regression"]] <- tmpResults
  if(bestR2<tmpR2){
    bestR2<-tmpR2
    bestModel <- tmpResults
  }
  regressionR2["linear regression",selection] <- tmpR2
  
  #support vector regression 
  tmpResults<-tune(svm,kreat_1~., data=as.data.frame(impute(selectedSubsets[[selection]]),what="mean"), ranges = list(gamma = 2^(-5:1), cost = 2^(-2:1),kernel=c("polynomial","radial")),
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
# sink()

#wyniki R2
write.csv(regressionR2,"results/regressionR2.csv")

bestRegressionResultIndex<- which(regressionR2== max(regressionR2,na.rm=TRUE), arr.ind = TRUE)

#wykresy, obrazujace tuning
dir.create(file.path("resultPlots"),showWarnings=FALSE)

for (selectionType in namesOfTests){
  counts <- c(regressionR2[,selectionType])
  png(paste("resultPlots/",gsub("[.]", "_", selectionType),"_R2",".png",sep=""),width=600,height=600)
  barplot(counts,names.arg=modelNames,col=c("darkblue","red","orange"),main=paste("Zależność R2 od wybranego modelu dla",selectionType ,sep=" "),ylab="Wartosc R2")
  dev.off()
}

countsR2<-numeric(0)
R2names<-numeric(0)
for (selectionType in namesOfTests){
  for (modelType in modelNames){
    countsR2<-c(countsR2,regressionR2[modelType,selectionType])
    R2names<- c(R2names, paste(modelType,selectionType,sep=" "))
  }
}


png(paste("resultPlots/","full","_R2",".png",sep=""),width=800,height=600)
par(las=2) # make label text perpendicular to axis
par(mar=c(5,16,4,2)) # increase y-axis margin.
barplot(countsR2,names.arg=R2names,horiz=TRUE,col=c("darkblue","red","orange"),main="Zależność R2 od wybranego modelu oraz metody selekcji atrybutow",xlab="Wartosc R2")
dev.off()

#####ocena najlepszego modelu oraz eksperyment#####


bestRegressionDataFrame <- as.data.frame(selectedSubsets[[namesOfTests[bestRegressionResultIndex[2]]]])
bestRegressionModel <- modelNames[bestRegressionResultIndex[1]]
bestParams<- regressionResults[namesOfTests[bestRegressionResultIndex[2]]]

testCvR2 <-vector()
testBootR2<-vector()
if(bestRegressionModel=="svr"){
  #model1 <- svm(kreat_1~.,data=bestRegressionDataFrame,gamma=bestModel$best.parameters$gamma,cost=bestModel$best.parameters$cost,kernel=bestModel$best.parameters$kernel)
  
  testData<-as.data.frame(impute(bestRegressionDataFrame,what="mean"))
  #walidacja krzyzowa
  testModel1<-tune(svm,kreat_1~.,data=testData,
                   gamma=bestModel$best.parameters$gamma,cost=bestModel$best.parameters$cost,kernel=bestModel$best.parameters$kernel,
                   tunecontrol=tune.control(nrepeat=cv_n, repeat.aggregate=median,sampling="cross", cross=10))
  testCvR2 <- (1 -( length(testData$kreat_1)*testModel1$best.performance/((length(testData$kreat_1)-1)*var(testData$kreat_1))))
  
  #bootstraping
  testModel1B<-tune(svm,kreat_1~.,data=testData,
                   gamma=bestModel$best.parameters$gamma,cost=bestModel$best.parameters$cost,kernel=bestModel$best.parameters$kernel,
                   tunecontrol=tune.control(sampling="boot", nboot = cv_n*10, boot.size = 1))
  testBootR2 <- (1 -( length(testData$kreat_1)*testModel1B$best.performance/((length(testData$kreat_1)-1)*var(testData$kreat_1))))
  
  
  
  for(kreat in kreat_names){
    
    kreatBase <- bestRegressionDataFrame[, -which(names(bestRegressionDataFrame) %in% kreat1_name)]
    testData <- as.data.frame(impute(CreateDfFromAttributes(fullDF=df_base,atts=colnames(kreatBase),additionalAtt=kreat),what="mean"))
    
    #walidacja krzyzowa
    testModel<- tune(svm,as.formula(paste(kreat,"~.")),data=testData,gamma=bestModel$best.parameters$gamma,cost=bestModel$best.parameters$cost,kernel=bestModel$best.parameters$kernel,
                     tunecontrol=tune.control(nrepeat=cv_n, repeat.aggregate=median,sampling="cross", cross=10))
    tmpR2<-(1 -( (dim(testData)[1])*testModel$best.performance/((dim(testData)[1]-1)*var(testData[,length(testData)]))))
    testCvR2<-c(testCvR2,tmpR2)
    
    #bootstraping
    testModelB<- tune(svm,as.formula(paste(kreat,"~.")),data=testData,gamma=bestModel$best.parameters$gamma,cost=bestModel$best.parameters$cost,kernel=bestModel$best.parameters$kernel,
                     tunecontrol=tune.control(sampling="boot", nboot = cv_n*10, boot.size = 1))
    tmpBootR2<-(1 -( (dim(testData)[1])*testModelB$best.performance/((dim(testData)[1]-1)*var(testData[,length(testData)]))))
    testBootR2<-c(testBootR2,tmpBootR2)
    
  }
}else if(bestRegressionModel=="regression tree"){
  testData<-as.data.frame(impute(bestRegressionDataFrame,what="mean"))
  #walidacja krzyzowa
  testModel1<-tune(svm,kreat_1~.,data=testData,
                   minsplit=bestModel$best.parameters$minsplit, cp=bestModel$best.parameters$cp,
                   tunecontrol=tune.control(nrepeat=cv_n, repeat.aggregate=median,sampling="cross", cross=10))
  testCvR2 <- (1 -( length(testData$kreat_1)*testModel1$best.performance/((length(testData$kreat_1)-1)*var(testData$kreat_1))))
  
  #bootstraping
  testModel1<-tune(svm,kreat_1~.,data=testData,
                   minsplit=bestModel$best.parameters$minsplit, cp=bestModel$best.parameters$cp,
                   tunecontrol=tune.control(sampling="boot", nboot = cv_n*10, boot.size = 1))
  testBootR2 <- (1 -( length(testData$kreat_1)*testModel1$best.performance/((length(testData$kreat_1)-1)*var(testData$kreat_1))))
  
  
  
  for(kreat in kreat_names){
    
    kreatBase <- bestRegressionDataFrame[, -which(names(bestRegressionDataFrame) %in% kreat1_name)]
    testData <- as.data.frame(impute(CreateDfFromAttributes(fullDF=df_base,atts=colnames(kreatBase),additionalAtt=kreat),what="mean"))
    
    #walidacja krzyzowa
    testModel<- tune(svm,as.formula(paste(kreat,"~.")),data=testData,minsplit=bestModel$best.parameters$minsplit, cp=bestModel$best.parameters$cp,
                     tunecontrol=tune.control(nrepeat=cv_n, repeat.aggregate=median,sampling="cross", cross=10))
    tmpR2<-(1 -( (dim(testData)[1])*testModel$best.performance/((dim(testData)[1]-1)*var(testData[,length(testData)]))))
    testCvR2<-c(testCvR2,tmpR2)
    
    #bootstraping
    testModel<- tune(svm,as.formula(paste(kreat,"~.")),data=testData,minsplit=bestModel$best.parameters$minsplit, cp=bestModel$best.parameters$cp,
                     tunecontrol=tune.control(sampling="boot", nboot = cv_n*10, boot.size = 1))
    tmpBootR2<-(1 -( (dim(testData)[1])*testModel$best.performance/((dim(testData)[1]-1)*var(testData[,length(testData)]))))
    testBootR2<-c(testBootR2,tmpBootR2)
    
  }
}

#wykresy dla eksperymentu

png(paste("resultPlots/","cvExperiment","_R2",".png",sep=""),width=600,height=600)
barplot(testCvR2,names.arg=c(kreat1_name,kreat_names),main="Zależność R2 dla wartosci kreatyniny w kolejnych dniach (dla CV-10)",ylab="Wartosc R2")
dev.off()
write.csv(testCvR2,"results/testCvR2.csv")

png(paste("resultPlots/","bootstrapExperiment","_R2",".png",sep=""),width=600,height=600)
barplot(testBootR2,names.arg=c(kreat1_name,kreat_names),main="Zależność R2 dla wartosci kreatyniny w kolejnych dniach (dla bootstrapu)",ylab="Wartosc R2")
dev.off()
write.csv(testBootR2,"results/testBootR2.csv")
