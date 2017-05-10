
importEnviroCar2 = function(file) {
  
  require(rjson) # fromJSON
  require(maptools) # spCbind
  require(rgdal) #readOGR
  require(RCurl) #getURL
  require(stringr) #str_replace_all
  
  #file <- url
  # read data as spatial object:
  layer = readOGR(getURL(file,ssl.verifypeer = FALSE), layer = "OGRGeoJSON")
  
  # convert time from text to POSIXct:
  #layer$time = as.POSIXct(layer$time, format="%Y/%m/%d %H:%M:%S", usetz = TRUE)
  layer$time = as.POSIXct(layer$time, format="%Y-%m-%dT%H:%M:%SZ")
  
  # the third column is JSON, we want it in a table (data.frame) form:
  # 1. form a list of lists
  l1 = lapply(as.character(layer[[3]]), fromJSON)
  # 2. parse the $value elements in the sublist:
  l2 = lapply(l1,function(x) as.data.frame(lapply(x, function(X) X$value)))
  # dynamic parsing of phenomenon names and units
  phenomenonsUrl = "https://envirocar.org/api/stable/phenomenons"
  phenomenons = fromJSON(getURL(phenomenonsUrl,ssl.verifypeer = FALSE))
  
  colNames <- c("GPS.Bearing", "GPS.HDOP", "GPS.Speed")
  if (!all(colNames %in% names(l2[[2]]))){ 
    stop("Trajectory does not contain all the necessary data (GPS.Bearing, GPS.HDOP, GPS.Speed)")
  }else{
    colNames <- names(l2[[2]])
  }
  
  
  resultMatrix = matrix(nrow = length(l2),ncol = length(colNames))
  dimnames(resultMatrix)[[2]] = colNames
  for (i in seq(along = l2))
    resultMatrix[i,colNames] = as.numeric(l2[[i]])[match(colNames, names(l2[[i]]))]
  result = as.data.frame(resultMatrix)
  
  # set the units:
  units <- sapply(phenomenons[[1]], "[[", "unit")
  names(units)=colNames
  
  # add a units attribute to layer
  layer[[3]] = NULL
  # add the table as attributes to the spatial object
  if (length(layer) == nrow(result)) {
    layer = spCbind(layer, result)
    attr(layer, "units") = units
    return(layer)
  } else
    return(NULL)
}
