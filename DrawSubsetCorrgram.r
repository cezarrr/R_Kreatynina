library(corrgram)

DrawSubsetCorrgram <- function(sub_set,name){
  corrgram(sub_set, order=NULL, lower.panel=panel.ellipse,
           upper.panel=panel.pie, text.panel=panel.txt,
           diag.panel=panel.minmax,
           main=paste("Korelogram podzbiorow bazy Kreatynina (",name,")",sep=""))
}