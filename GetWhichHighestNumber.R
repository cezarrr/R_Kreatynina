GetWhichHighestNumber<-function (df){
  toReturn<-as.data.frame(setNames(replicate(length(colnames(df)),numeric(0), simplify = F), colnames(df)))
  for (i in 1:length(colnames(df))){
    toReturn[1,i]<-row.names(df)[which.max(df[,i])]
  }
  return (toReturn)
}