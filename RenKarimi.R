#Universidade Federal de Minas Gerais - Brasil
#LITC - Computational Intelligence and Technology Laboratory. http://litc.cpdee.ufmg.br.
#Code developed by:
#Jose Maia Neto - jmnt@ufmg.br, Juan Camilo Fonseca Galindo - juankmilofg@gmail.com
#Code based on article:
# Quddus, M. A., Noland, R. B. and Ochieng, W. Y. (2006). 
# Ren, M. and Karimi, H. A. (2009). 
#A chain-code-based map matching algorithm for wheelchair navigation.
#T. GIS 13(2): 197–214.

#### METODO Ren and Karimi 2009
#Map Matching Ren and Karimi
#Entrada: 
#   geodf: Rota do Usuario, data.fram com lat, lon, heading
#   ErroGPS: Error do GPS
#   way.map: rotas ou caminhos do mapa com  from_node_id, to_node_id, weights
#   nodes.map: nó do mapa com  id, lon, lat
#   printparam: bool, imprimir o nó da lista e a quantidade de caminhos candidatos
#   plot_box: bool, plotear quadrado no mapa
#   mode: 1-linear, 2-quadratio 
#Saida:
# rota_true: Latitude e Longitude do GPS, 
#            Latitude e Longitude calculadas, 
#            Semelhança (do ponto projetado e o ponto do GPS),
#            Probabilidade de esta na rota certa,
#            o ponto esta no Way Nertwork,   variavel binaria

# Load required packages
require(XML)
require(OpenStreetMap)
require(lubridate)
require(raster)
require(osmar)
require(dtw)

# File address
setwd('~/Dropbox/Trajectory Data Mining/repositorios')

# loading functions
source('map_maching.R')

Map_matching_Ren_and_Karimi <- function(geodf,ErroGPS,way.map,nodes.map,printparam,plot_box,mode)
{
  rota_true <- data.frame(lat_GPS = as.numeric(), lon_GPS = as.numeric(), lat_cal = as.numeric(), lon_cal = as.numeric(), id_way = as.numeric() ,matching =as.numeric() )
  
  if (mode==1)
  {
    #parametros metodo linear
    Wij  <- 0.5
    Wij1 <- 0.1
    Wij2 <- 0.2
    Wij3 <- 0.2
    
    for(ponto in c(1:length(geodf[all(TRUE),1]-2))) 
    {
      if(printparam){print(ponto)}
      lat_ponto_aux <- geodf[ponto,'lat']
      lon_ponto_aux <- geodf[ponto,'lon']  
      heading_ponto_aux1<- geodf[ponto,'heading']
      heading_ponto_aux2<- geodf[ponto+1,'heading']
      heading_ponto_aux3<- geodf[ponto+2,'heading']
      
      #Limitar a busqueda em um recuadro de 2 veces a maior comprimento das ruas
      box_way_limites <- box_way_lim(lat_ponto_aux,lon_ponto_aux,way.map,nodes.map,geodf,plot_box)
      #Lista de candidatos
      if(length(box_way_limites$way_id))
      {
        way_candidates<-list_candidatos(lat_ponto_aux,lon_ponto_aux,box_way_limites,nodes.map,ErroGPS)
        
        rota_candidatos <- data.frame(lat_cal = as.numeric(), lon_cal = as.numeric(), matching =as.numeric() )
        
        if (!is.null(way_candidates))
        {
          for(way_candidato in c(1:length(way_candidates[all(TRUE),1])))
          {
            from_node_aux<- subset(nodes.map, id == way_candidates$from_node_id[way_candidato])
            to_node_aux<- subset(nodes.map, id == way_candidates$to_node_id[way_candidato])
            id_cal_way_aux<- way_candidates$way_id
            
            D_Ren <- DistMin_way_point(from_node_aux$lat[1],from_node_aux$lon[1],to_node_aux$lat[1],to_node_aux$lon[1],lat_ponto_aux,lon_ponto_aux)
            Dcc_Ren1 <- Dcc(from_node_aux$lat[1],from_node_aux$lon[1],to_node_aux$lat[1],to_node_aux$lon[1],heading_ponto_aux1)
            Dcc_Ren2 <- Dcc(from_node_aux$lat[1],from_node_aux$lon[1],to_node_aux$lat[1],to_node_aux$lon[1],heading_ponto_aux2)
            Dcc_Ren3 <- Dcc(from_node_aux$lat[1],from_node_aux$lon[1],to_node_aux$lat[1],to_node_aux$lon[1],heading_ponto_aux3)
            
            matching_aux<- 1 / ( Wij*D_Ren[1] + sum ( Wij1*Dcc_Ren1, Wij2*Dcc_Ren2, Wij3*Dcc_Ren3) )
            lat_cal_aux<-D_Ren[2]
            lon_cal_aux<-D_Ren[3]
            rota_candidatos_aux <- data.frame(lat_cal = lat_cal_aux, lon_cal = lon_cal_aux, id_way = id_cal_way_aux ,matching =matching_aux )
            rota_candidatos<- rbind(rota_candidatos,rota_candidatos_aux)
            rm(list=c('from_node_aux', 'to_node_aux', 'D_Ren1', 'Dcc_Ren1','Dcc_Ren2', 'Dcc_Ren3','matching_aux','lat_cal_aux' ,'lon_cal_aux', 'rota_candidatos_aux'))
            
          }
          
          candidato_true<-max(rota_candidatos$matching)
          rota_candidato <- data.frame(lat_cal = as.numeric(), lon_cal = as.numeric(), id_way = as.numeric() ,matching =as.numeric() )
          rota_candidato <- rbind (rota_candidato,subset(rota_candidatos, matching == candidato_true))
          
          rota_true1 <- data.frame(lat_GPS = lat_ponto_aux, lon_GPS = lon_ponto_aux, lat_cal = rota_candidato$lat_cal[1], lon_cal = rota_candidato$lon_cal[1], id_way = rota_candidato$id_way[1], matching = rota_candidato$matching[1])
          rota_true<- rbind(rota_true,rota_true1)
          
          rm(list=c('candidato_true','rota_candidato','rota_true1'))
        }else{
          rota_true1 <- data.frame(lat_GPS = lat_ponto_aux, lon_GPS = lon_ponto_aux, lat_cal = lat_ponto_aux, lon_cal = lon_ponto_aux, id_way = NA, matching = NA)
          rota_true<- rbind(rota_true,rota_true1)
        }
      }else{
        rota_true1 <- data.frame(lat_GPS = lat_ponto_aux, lon_GPS = lon_ponto_aux, lat_cal = lat_ponto_aux, lon_cal = lon_ponto_aux, id_way = NA, matching = NA)
        rota_true<- rbind(rota_true,rota_true1)}
      
      
      
    }
  }else if (mode==2)
  {
    #parametros metodo linear
    u<-0
    sigma_D_Ren <- ErroGPS
    sigma_Dcc_Ren1 <- 2.5
    sigma_Dcc_Ren2 <- 3
    sigma_Dcc_Ren3 <- 3.5
    
    for(ponto in c(1:length(geodf[all(TRUE),1]-2))) 
    {
      if(printparam){print(ponto)}
      lat_ponto_aux <- geodf[ponto,'lat']
      lon_ponto_aux <- geodf[ponto,'lon']  
      heading_ponto_aux1<- geodf[ponto,'heading']
      heading_ponto_aux2<- geodf[ponto+1,'heading']
      heading_ponto_aux3<- geodf[ponto+2,'heading']
      
      #Limitar a busqueda em um recuadro de 2 veces a maior comprimento das ruas
      box_way_limites <- box_way_lim(lat_ponto_aux,lon_ponto_aux,way.map,nodes.map,geodf,plot_box)
      if(length(box_way_limites$way_id))
      {
        way_candidates<-list_candidatos(lat_ponto_aux,lon_ponto_aux,box_way_limites,nodes.map,ErroGPS)
        
        rota_candidatos <- data.frame(lat_cal = as.numeric(), lon_cal = as.numeric(), matching =as.numeric() )
        
        if (!is.null(way_candidates))
        {
          for(way_candidato in c(1:length(way_candidates[all(TRUE),1])))
          {
            from_node_aux<- subset(nodes.map, id == way_candidates$from_node_id[way_candidato])
            to_node_aux<- subset(nodes.map, id == way_candidates$to_node_id[way_candidato])
            id_cal_way_aux<- way_candidates$way_id
            
            D_Ren <- DistMin_way_point(from_node_aux$lat[1],from_node_aux$lon[1],to_node_aux$lat[1],to_node_aux$lon[1],lat_ponto_aux,lon_ponto_aux)
            Dcc_Ren1 <- Dcc(from_node_aux$lat[1],from_node_aux$lon[1],to_node_aux$lat[1],to_node_aux$lon[1],heading_ponto_aux1)
            Dcc_Ren2 <- Dcc(from_node_aux$lat[1],from_node_aux$lon[1],to_node_aux$lat[1],to_node_aux$lon[1],heading_ponto_aux2)
            Dcc_Ren3 <- Dcc(from_node_aux$lat[1],from_node_aux$lon[1],to_node_aux$lat[1],to_node_aux$lon[1],heading_ponto_aux3)
            
            D_Ren_fuzzy    <-  exp(- ( abs(  D_Ren[1] - u[1] )^2 ) / ( 2 * sigma_D_Ren[1]^2    ) )
            Dcc_Ren1_fuzzy <- exp(- ( abs(Dcc_Ren1[1]-u[1])^2 ) / ( 2 * sigma_Dcc_Ren1[1]^2 ) )
            Dcc_Ren2_fuzzy <- exp(- ( abs(Dcc_Ren2[1]-u[1])^2 ) / ( 2 * sigma_Dcc_Ren2[1]^2 ) )
            Dcc_Ren3_fuzzy <- exp(- ( abs(Dcc_Ren3[1]-u[1])^2 ) / ( 2 * sigma_Dcc_Ren3[1]^2 ) )
            
            matching_aux<- min(D_Ren_fuzzy,Dcc_Ren1_fuzzy,Dcc_Ren2_fuzzy,Dcc_Ren3_fuzzy)
            
            lat_cal_aux<-D_Ren[2]
            lon_cal_aux<-D_Ren[3]
            rota_candidatos_aux <- data.frame(lat_cal = lat_cal_aux, lon_cal = lon_cal_aux, id_way = id_cal_way_aux, matching =matching_aux )
            rota_candidatos<- rbind(rota_candidatos,rota_candidatos_aux)
            rm(list=c('from_node_aux', 'to_node_aux', 'D_Ren1', 'Dcc_Ren1','Dcc_Ren2', 'Dcc_Ren3','matching_aux','lat_cal_aux' ,'lon_cal_aux', 'rota_candidatos_aux'))
            
          }
          
          candidato_true<-max(rota_candidatos$matching)
          rota_candidato <- data.frame(lat_cal = as.numeric(), lon_cal = as.numeric(), id_way = as.numeric(), matching =as.numeric() )
          rota_candidato <- rbind (rota_candidato,subset(rota_candidatos, matching == candidato_true))
          rota_true1 <- data.frame(lat_GPS = lat_ponto_aux, lon_GPS = lon_ponto_aux, lat_cal = rota_candidato$lat_cal[1], lon_cal = rota_candidato$lon_cal[1], id_way = rota_candidato$id_way[1],matching = rota_candidato$matching[1])
          rota_true<- rbind(rota_true,rota_true1)
          
          rm(list=c('candidato_true','rota_candidato','rota_true1'))
        }else{
          rota_true1 <- data.frame(lat_GPS = lat_ponto_aux, lon_GPS = lon_ponto_aux, lat_cal = lat_ponto_aux, lon_cal = lon_ponto_aux, id_way = NA, matching = NA)
          rota_true<- rbind(rota_true,rota_true1)
        }
      }else{
        rota_true1 <- data.frame(lat_GPS = lat_ponto_aux, lon_GPS = lon_ponto_aux, lat_cal = lat_ponto_aux, lon_cal = lon_ponto_aux, id_way = NA, matching = NA)
        rota_true<- rbind(rota_true,rota_true1)
      }
      #Lista de candidatos
      
      
    }
  }
  
  return(rota_true)
}