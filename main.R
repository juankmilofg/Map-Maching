#Universidade Federal de Minas Gerais - Brasil
#LITC - Computational Intelligence and Technology Laboratory. http://litc.cpdee.ufmg.br.
#Code developed by:
#Jose Maia Neto - jmnt@ufmg.br, Juan Camilo Fonseca Galindo - juankmilofg@gmail.com

# clear all
rm(list = ls())

# Load required packages
library(XML)
library(OpenStreetMap)
library(lubridate)
library(raster)
library(osmar)
library(dtw)

# File address
setwd('~/Dropbox/Trajectory Data Mining/repositorios')

# loading functions
source('map_maching.R')
source('Quddus.R')
source('Jagadeesh.R')
source('RenKarimi.R')


# Choose the route
N_trajectoria <- 1

# loadind ID route reference
Reftraj <- "RefTraj"
Reftraj <- paste(Reftraj,as.character(N_trajectoria),sep = "")
Reftraj <- paste(Reftraj,"gpx",sep = ".")

pfileRef <- htmlTreeParse(Reftraj,
                          error = function (...) {},
                          useInternalNodes = T)

Ref.times <- xpathSApply(pfileRef, path = '//trkpt/time', xmlValue)
Ref.coords <- xpathSApply(pfileRef, path = '//trkpt', xmlAttrs)

# longitude and latitude coordinates route
lon_ref <- Ref.coords[2,]
lat_ref <- Ref.coords[1,]

# Remove variables
rm(list=c('Reftraj', 'pfileRef','Ref.times','Ref.coords'))

load("NovasTrajetorias.RData")

idx.Traj <- TrajsIdx[[N_trajectoria]]
# URL of the trajectory
url = paste("https://envirocar.org/api/stable/tracks/", idx.Traj, sep = "")

# Import the trajectory
source('importEnviroCar2.R')
trajetoria <- tryCatch(importEnviroCar2(url),error=function(c)"Error")

if (is.object(trajetoria))
{
  # Limit the number of points on the route  (<= 150 point GPS)
  trajetoria <- trajetoria[!is.na(trajetoria$GPS.Bearing),]
  if (length(trajetoria) > 150){    tam <- 150
  }else{tam <- length(trajetoria)}
  trajetoria <- trajetoria[!is.na(trajetoria$GPS.Bearing),][1:tam,] 
  
  # loading data extra
  load("Trajetorias_DadosExtras.RData")
  trajetoria_full <- as.data.frame(trajetoria)
  
  # Extract times,latitude, longitude and heading angle
  times   <- as.numeric(trajetoria_full$time)
  lats    <- as.numeric(as.numeric(trajetoria_full$coords.x2))
  lons    <- as.numeric(as.numeric(trajetoria_full$coords.x1))
  if (trajetoria_full$GPS.Bearing<180){heading <- as.numeric(trajetoria_full$GPS.Bearing)
  }else{heading <- as.numeric(trajetoria_full$GPS.Bearing) - 360}
  
  # Put everything in a dataframe and get rid of old variables
  geodf <- data.frame(lat = lats, lon = lons, time = times, heading = heading)
  geodf <- subset(geodf, is.finite(heading))
  
  # Remove variables
  rm(list=c('lats', 'lons', 'times', 'heading'))
  
  # limit road network 
  cent_lat<-(max(geodf$lat)+min(geodf$lat))/2
  cent_lon<-(max(geodf$lon)+min(geodf$lon))/2
  dist_lat<-geodist(as.numeric(min(geodf$lat)),as.numeric(cent_lon),
                    as.numeric(max(geodf$lat)),as.numeric(cent_lon)) 
  dist_lon<-geodist(as.numeric(cent_lat),as.numeric(min(geodf$lon)),
                    as.numeric(cent_lat),as.numeric(max(geodf$lon))) 
  
  muc_bbox <- center_bbox(cent_lon, cent_lat, dist_lon*1.05, dist_lat*1.05)
  
  # Remove variables
  rm(list=c('cent_lat', 'cent_lon', 'dist_lat', 'dist_lon'))
  
  muc <- tryCatch(get_osm(muc_bbox),error=function(c)"Error")
  if (muc == "Error")
  {
    print("mais de 5000 points")
  }else{
    hways_muc <- subset(muc, way_ids = find(muc, way(tags(k == "highway"))))
    hways <- find(hways_muc, way(tags(k == "name")))
    hways <- find_down(muc, way(hways))
    hways_muc <- subset(muc, ids = hways)
    
    # Remove variables
    rm(list=c('muc'))
    
    # Road network way.map(from_node_id, to_node_id, way_id, weights)
    way.map <- tryCatch(as.data.frame(as_igraph2(hways_muc)),error=function(c)"Error")
    if (way.map == "Error")
    {
      print("way.map null")
    }else{
      # Nodes listof road network  nodes.map(ids, lat, lon)
      nodes.map <- data.frame(id = hways_muc$nodes$attrs$id, lat = hways_muc$nodes$attrs$lat,
                              lon = hways_muc$nodes$attrs$lon)
      
      # Clean reference route
      ID_true <- way_ref(lat_ref,lon_ref,nodes.map,way.map,geodf)
      
      
      ##### Original Data ####
      print('Original Data')
      rota.dados<- data.frame(lat_gps=geodf$lat, lon_gps=geodf$lon ,lat_cal=geodf$lat, lon_cal=geodf$lon)
      
      ## Teste Average distance
      Test_dist_media_dadosBrutos <- test_dist(lat_ref, lon_ref, rota.dados,0)
      print(paste("Average distance:",Test_dist_media_dadosBrutos))
      ## Teste Medium distance
      Test_dist_mediana_dadosBrutos <- test_dist(lat_ref, lon_ref, rota.dados,1)
      print(paste("Medium distance:",Test_dist_mediana_dadosBrutos))
      ## Teste Summation distance
      Test_dist_sum_dadosBrutos <- test_dist(lat_ref, lon_ref, rota.dados,3)
      print(paste("Summation distance:",Test_dist_sum_dadosBrutos))
      ## Teste Exponential distance
      Test_dist_exp_dadosBrutos <- test_dist(lat_ref, lon_ref, rota.dados,4)
      print(paste("Exponential distance:",Test_dist_exp_dadosBrutos))
      ## Teste DTW
      Test_dist_DTW_dadosBrutos <- DTW.dist(lat_ref, lon_ref, rota.dados)
      print(paste("DTW:",Test_dist_DTW_dadosBrutos))
      
      
      #### LINEAR METHOD Ren and Karimi 2009 ####
      print("LINEAR METHOD Ren and Karimi 2009")
      
      # Parameters
      ErroGPS <- 20
      plot_box <- FALSE
      printparam <- FALSE
      mode<-1  # 1-linear, 2-Not linear
      
      rota_true_Ren_and_Karimi_1 <- Map_matching_Ren_and_Karimi(geodf,ErroGPS,way.map,nodes.map,printparam,plot_box,mode)
      
      ## Test matching ID's
      Test_road_Ren_and_Karimi_1 <- test_dist_road(ID_true$lat, ID_true$lon, ID_true$way, rota_true_Ren_and_Karimi_1) 
      print(paste("Matching ID's:",Test_road_Ren_and_Karimi_1))
      ## Teste Average distance
      Test_dist_media_Ren_and_Karimi_1 <- test_dist(lat_ref, lon_ref, rota_true_Ren_and_Karimi_1,0)
      print(paste("Average distance:",Test_dist_media_Ren_and_Karimi_1))
      ## Teste Medium distance
      Test_dist_mediana_Ren_and_Karimi_1 <- test_dist(lat_ref, lon_ref, rota_true_Ren_and_Karimi_1,1)
      print(paste("Medium distance:",Test_dist_mediana_Ren_and_Karimi_1))
      ## Teste Summation distance
      Test_dist_sum_Ren_and_Karimi_1 <- test_dist(lat_ref, lon_ref, rota_true_Ren_and_Karimi_1,3)
      print(paste("Summation distance:",Test_dist_sum_Ren_and_Karimi_1))
      ## Teste Exponential distance
      Test_dist_exp_Ren_and_Karimi_1 <- test_dist(lat_ref, lon_ref, rota_true_Ren_and_Karimi_1,4)
      print(paste("Exponential distance:",Test_dist_exp_Ren_and_Karimi_1))
      ## Teste DTW
      Test_dist_DTW_Ren_and_Karimi_1 <- DTW.dist(lat_ref, lon_ref, rota_true_Ren_and_Karimi_1)
      print(paste("DTW:",Test_dist_DTW_Ren_and_Karimi_1))
      
      #### NOT LINEAR METHOD Ren and Karimi 2009 ####
      print("NOT INEAR METHOD Ren and Karimi 2009")
      
      # Parameters
      ErroGPS <- 20
      plot_box <- FALSE
      printparam <- FALSE
      mode<-1  # 1-linear, 2-Not linear
      
      rota_true_Ren_and_Karimi_2 <- Map_matching_Ren_and_Karimi(geodf,ErroGPS,way.map,nodes.map,printparam,plot_box,mode)
      
      ## Test matching ID's
      Test_road_Ren_and_Karimi_2 <- test_dist_road(ID_true$lat, ID_true$lon, ID_true$way, rota_true_Ren_and_Karimi_2) 
      print(paste("Matching ID's:",Test_road_Ren_and_Karimi_2))
      ## Teste Average distance
      Test_dist_media_Ren_and_Karimi_2 <- test_dist(lat_ref, lon_ref, rota_true_Ren_and_Karimi_2,0)
      print(paste("Average distance:",Test_dist_media_Ren_and_Karimi_2))
      ## Teste Medium distance
      Test_dist_mediana_Ren_and_Karimi_2 <- test_dist(lat_ref, lon_ref, rota_true_Ren_and_Karimi_2,1)
      print(paste("Medium distance:",Test_dist_mediana_Ren_and_Karimi_2))
      ## Teste Summation distance
      Test_dist_sum_Ren_and_Karimi_2 <- test_dist(lat_ref, lon_ref, rota_true_Ren_and_Karimi_2,3)
      print(paste("Summation distance:",Test_dist_sum_Ren_and_Karimi_2))
      ## Teste Exponential distance
      Test_dist_exp_Ren_and_Karimi_2 <- test_dist(lat_ref, lon_ref, rota_true_Ren_and_Karimi_2,4)
      print(paste("Exponential distance:",Test_dist_exp_Ren_and_Karimi_2))
      ## Teste DTW
      Test_dist_DTW_Ren_and_Karimi_2 <- DTW.dist(lat_ref, lon_ref, rota_true_Ren_and_Karimi_2)
      print(paste("DTW:",Test_dist_DTW_Ren_and_Karimi_2))
      
      
      #### METHOD Jagadeesh 2004 ####
      print('METHOD Jagadeesh 2004')
      
      # Parameters
      ErroGPS <- 20
      plot_box <- FALSE
      printparam <- FALSE
      
      rota_true_Jagadeesh <- Map_matching_Jagadeesh(geodf,ErroGPS,way.map,nodes.map,printparam,plot_box)

      ## Test matching ID's
      Test_road_Jagadeesh <- test_dist_road(ID_true$lat, ID_true$lon, ID_true$way, rota_true_Jagadeesh) 
      print(paste("Matching ID's:",Test_road_Jagadeesh))
      ## Teste Average distance
      Test_dist_media_Jagadeesh <- test_dist(lat_ref, lon_ref, rota_true_Jagadeesh,0)
      print(paste("Average distance:",Test_dist_media_Jagadeesh))
      ## Teste Medium distance
      Test_dist_mediana_Jagadeesh <- test_dist(lat_ref, lon_ref, rota_true_Jagadeesh,1)
      print(paste("Medium distance:",Test_dist_mediana_Jagadeesh))
      ## Teste Summation distance
      Test_dist_sum_Jagadeesh <- test_dist(lat_ref, lon_ref, rota_true_Jagadeesh,3)
      print(paste("Summation distance:",Test_dist_sum_Jagadeesh))
      ## Teste Exponential distance
      Test_dist_exp_Jagadeesh <- test_dist(lat_ref, lon_ref, rota_true_Jagadeesh,4)
      print(paste("Exponential distance:",Test_dist_exp_Jagadeesh))
      ## Teste DTW
      Test_dist_DTW_Jagadeesh <- DTW.dist(lat_ref, lon_ref, rota_true_Jagadeesh)
      print(paste("DTW:",Test_dist_DTW_Jagadeesh))
      
      
      #### METHOD Quddus 2006 ####
      print('METHOD Quddus 2006')
      
      mmTraj <- fuzzyMMfct(trajetoria)
      rota_true_fuzzyMMfc <- as.data.frame(mmTraj)
      
      ## Test matching ID's
      Test_road_fuzzyMMfc <- test_dist_road(ID_true$lat, ID_true$lon, ID_true$way, rota_true_fuzzyMMfc) 
      print(paste("Matching ID's:",Test_road_fuzzyMMfc))
      ## Teste Average distance
      Test_dist_media_fuzzyMMfc <- test_dist(lat_ref, lon_ref, rota_true_fuzzyMMfc,0)
      print(paste("Average distance:",Test_dist_media_fuzzyMMfc))
      ## Teste Medium distance
      Test_dist_mediana_fuzzyMMfc <- test_dist(lat_ref, lon_ref, rota_true_fuzzyMMfc,1)
      print(paste("Medium distance:",Test_dist_mediana_fuzzyMMfc))
      ## Teste Summation distance
      Test_dist_sum_fuzzyMMfc <- test_dist(lat_ref, lon_ref, rota_true_fuzzyMMfc,3)
      print(paste("Summation distance:",Test_dist_sum_fuzzyMMfc))
      ## Teste Exponential distance
      Test_dist_exp_fuzzyMMfc <- test_dist(lat_ref, lon_ref, rota_true_fuzzyMMfc,4)
      print(paste("Exponential distance:",Test_dist_exp_fuzzyMMfc))
      ## Teste DTW
      Test_dist_DTW_fuzzyMMfc <- DTW.dist(lat_ref, lon_ref, rota_true_fuzzyMMfc)
      print(paste("DTW:",Test_dist_DTW_fuzzyMMfc))
      
      
      ###################plot############

      lat.leftupper <- as.numeric(muc_bbox[4]) 
      lat.rightlower <- as.numeric(muc_bbox[2])
      lon.leftupper <- as.numeric(muc_bbox[1]) 
      lon.rightlower <- as.numeric(muc_bbox[3])
      
      x11()
      #Obtém mapa
      map <- openmap(c(lat.leftupper, lon.leftupper),
                     c(lat.rightlower, lon.rightlower), type = 'maptoolkit-topo')
      #Plota o mapa
      transmap <- openproj(map, projection = '+proj=longlat')
      plot(transmap, raster=T)
      
      plot_ways(hways_muc, add = TRUE)
      plot_nodes(hways_muc, add = TRUE, col = "black")
      
      legend("topleft",c("Nós","Segmentos"),pch = c(1,NA),
             lty = c(NA,1),cex = 1)
      

      x11()
      plot(transmap, raster=T)
      par(new = T)
      points(geodf$lon,geodf$lat)
      lines(rota_true_Ren_and_Karimi_1$lon_GPS,rota_true_Ren_and_Karimi_1$lat_GPS)
      title(main = "NOT LINEAR METHOD Ren and Karimi 2009")
      
      x11()
      plot(transmap, raster=T)
      par(new = T)
      points(geodf$lon,geodf$lat)
      lines(rota_true_Ren_and_Karimi_2$lon_GPS,rota_true_Ren_and_Karimi_2$lat_GPS)
      title(main = "LINEAR METHOD Ren and Karimi 2009")

      x11()
      plot(transmap, raster=T)
      par(new = T)
      points(geodf$lon,geodf$lat)
      lines(rota_true_Jagadeesh$lon_GPS,rota_true_Jagadeesh$lat_GPS)
      title(main = "METHOD Jagadeesh 2004")
      
      x11()
      plot(transmap, raster=T)
      par(new = T)
      points(geodf$lon,geodf$lat)
      lines(rota_fuzzyMMfc$lon_GPS,rota_fuzzyMMfc$lat_GPS)
      title(main = "METHOD Quddus 2006")

    }
  }
}




