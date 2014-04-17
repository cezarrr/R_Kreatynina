#CreateDfFromAttributes
#tworzy nowy dataframe zawierajacy wybrane atrybuty - na podstawie przekazanego
#parametry: fullDF - pelny dataframe, na podstawie ktorego tworzymy nowy
#           atts - nazwy wybranych atrybutow
#           additionalAtt - opcjonalny dodatkowy atrybut
#zwraca stworzony dataframe
CreateDfFromAttributes<-function(fullDF,atts,additionalAtt=""){
  if(additionalAtt==""){
    toReturn <- fullDF[, which(names(fullDF) %in% atts)]
  }
  else{
    toReturn <-fullDF[, which(names(fullDF) %in% c(atts,additionalAtt))]
  }
  return (toReturn)
}