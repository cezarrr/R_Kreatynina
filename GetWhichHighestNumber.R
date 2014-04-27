#GetWhichHighestNumber
#funkcja zwracajaca indeks wystapienia wartosci maksymalnej dla kazdej kolumny data frame
#parametry: df - dataframe
#zwraca dataframe z wierszem, w ktorym dla kazdej kolumny wypisany jest indeks wystapienia najwiekszej wartosci
GetWhichHighestNumber<-function (df){
  toReturn<-as.data.frame(setNames(replicate(length(colnames(df)),numeric(0), simplify = F), colnames(df)))
  for (i in 1:length(colnames(df))){
    if(any(is.na(df[,i]))){
      toReturn[1,i]<-"NA value"
    }
    toReturn[1,i]<-row.names(df)[which.max(df[,i])]
  }
  return (toReturn)
}