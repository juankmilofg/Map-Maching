#Universidade Federal de Minas Gerais - Brasil
#LITC - Computational Intelligence and Technology Laboratory. http://litc.cpdee.ufmg.br.
#Code developed by:
#Jose Maia Neto - jmnt@ufmg.br, Juan Camilo Fonseca Galindo - juankmilofg@gmail.com
#Code based on article:
# Jagadeesh, G. R., Srikanthan, T. and Zhang, X. D.(2004). 
#A map matching method for gps based real-time vehicle location
#Journal of Navigation 57(3): 429–440.


#### METODO JAGADEESH 2004
#Map Matching Jagadeesh
#Entrada: 
#   geodf: Rota do Usuario, data.fram com lat, lon, heading
#   ErroGPS: Error do GPS
#   way.map: rotas ou caminhos do mapa com  from_node_id, to_node_id, weights
#   nodes.map: nó do mapa com  id, lon, lat
#   printparam: bool, imprimir o nó da lista e a quantidade de caminhos candidatos
#   plot_box: bool, plotear quadrado no mapa
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

Map_matching_Jagadeesh <- function(geodf,ErroGPS,way.map,nodes.map,printparam,plot_box)
{
  rota_true <- data.frame(lat_GPS = as.numeric(), lon_GPS = as.numeric(), lat_cal = as.numeric(), lon_cal = as.numeric(), id_way = as.numeric() ,semelhanza =as.numeric(), road_true =as.numeric(), off_road=as.logical() )
  for(ponto in c(1:length(geodf[all(TRUE),1]))) 
  {
    if(printparam){print(ponto)}
    lat_ponto_aux <- geodf[ponto,'lat']
    lon_ponto_aux <- geodf[ponto,'lon']  
    heading_ponto_aux<- geodf[ponto,'heading']
    
    #Limitar a busqueda em um recuadro de 2 veces a maior comprimento das ruas
    box_way_limites <- box_way_lim(lat_ponto_aux,lon_ponto_aux,way.map,nodes.map,geodf,plot_box)
    
    #Lista de candidatos
    if (length(box_way_limites$way_id)==0) {way_candidates<-NULL
    }else {
      way_candidates<-list_candidatos(lat_ponto_aux,lon_ponto_aux,box_way_limites,nodes.map,ErroGPS)
    }
    if (is.null(way_candidates))                         # si não retorna way perto
    {
      if(printparam){print('way_candidates is null')}
      off_road_aux<-TRUE
      lat_cal_aux<-lat_ponto_aux
      lon_cal_aux<-lon_ponto_aux
      id_way_aux<- NA
      semelhanza_aux<-1
    } else if(length(way_candidates[all(TRUE),1]) == 1)  # se só retorna um way candidata
    {
      if(printparam){print('way_candidates is 1')}
      off_road_aux<-FALSE
      from_node_aux<- subset(nodes.map, id == way_candidates$from_node_id)
      to_node_aux<- subset(nodes.map, id == way_candidates$to_node_id)
      id_way_aux<- way_candidates$way_id
      
      dist_aux <- DistMin_way_point(from_node_aux$lat,from_node_aux$lon,to_node_aux$lat,to_node_aux$lon,lat_ponto_aux,lon_ponto_aux)
      if (dist_aux[1]<50){dist_aux_fuzzy=dist_aux[1]*(-1/50)+1
      }else{dist_aux_fuzzy=0}
      
      heading_aux <- heading_way_point_abs(from_node_aux$lat,from_node_aux$lon,to_node_aux$lat,to_node_aux$lon,heading_ponto_aux)
      if (heading_aux<50){heading_aux_fuzzy=heading_aux*(-1/90)+1
      }else{heading_aux_fuzzy=0}
      
      semelhanza_aux<-min(dist_aux_fuzzy[1],heading_aux_fuzzy)
      lat_cal_aux<-dist_aux[2]
      lon_cal_aux<-dist_aux[3]
      rm(list=c('from_node_aux', 'to_node_aux', 'dist_aux', 'dist_aux_fuzzy','heading_aux','heading_aux_fuzzy'))
      
    } else                                               # se tem mais ways candidatasretorna um way
    {
      if(printparam){print('way_candidates is 2 or more')}
      rota_candidatos <- data.frame(lat_cal = as.numeric(), lon_cal = as.numeric(), semelhanza =as.numeric() )
      for(candidato in c(1:length(way_candidates[all(TRUE),1])))
      {
        from_node_aux<- subset(nodes.map, id == way_candidates$from_node_id[candidato])
        to_node_aux<- subset(nodes.map, id == way_candidates$to_node_id[candidato])
        
        dist_aux <- DistMin_way_point(from_node_aux$lat[1],from_node_aux$lon[1],to_node_aux$lat[1],to_node_aux$lon[1],lat_ponto_aux,lon_ponto_aux)
        if (dist_aux[1]<50){dist_aux_fuzzy=dist_aux[1]*(-1/50)+1
        }else{dist_aux_fuzzy<-0}
        
        heading_aux <- heading_way_point_abs(from_node_aux$lat[1],from_node_aux$lon[1],to_node_aux$lat[1],to_node_aux$lon[1],heading_ponto_aux)
        if (heading_aux<90){heading_aux_fuzzy=heading_aux*(-1/180)+1
        }else{heading_aux_fuzzy=0}
        
        semelhanza_aux<-min(dist_aux_fuzzy[1],heading_aux_fuzzy)
        lat_cal_aux<-dist_aux[2]
        lon_cal_aux<-dist_aux[3]
        id_cal_way_aux<- way_candidates$way_id
        rota_candidatos_aux <- data.frame(lat_cal = lat_cal_aux, lon_cal = lon_cal_aux, id_way= id_cal_way_aux, semelhanza =semelhanza_aux)
        rota_candidatos<- rbind(rota_candidatos,rota_candidatos_aux)
        rm(list=c('from_node_aux', 'to_node_aux', 'dist_aux', 'dist_aux_fuzzy','heading_aux','heading_aux_fuzzy'))
      }
      
      candidato_true<-max(rota_candidatos$semelhanza)
      rota_candidato <- data.frame(lat_cal = as.numeric(), lon_cal = as.numeric(), semelhanza =as.numeric() , way_id= as.numeric())
      rota_candidato <- rbind (rota_candidato,subset(rota_candidatos, semelhanza == candidato_true))
      
      lat_cal_aux<-rota_candidato$lat_cal[1]
      lon_cal_aux<-rota_candidato$lon_cal[1]
      id_way_aux<- rota_candidato$id_way[1]
      semelhanza_aux<-rota_candidato$semelhanza[1]
      off_road_aux<-FALSE
      
      rm(list=c('candidato_true','rota_candidato'))
    }
    
    road_true_aux <- 0
    if (ponto==1){road_true_aux <- 0
    }else if (rota_true$semelhanza[ponto-1]<0.01){road_true_aux <- semelhanza_aux
    }else{road_true_aux <- min (semelhanza_aux,rota_true$semelhanza[ponto-1])}
    
    if (road_true_aux < 0.3){off_road_aux <- TRUE}
    
    rota_true1 <- data.frame(lat_GPS = lat_ponto_aux, lon_GPS = lon_ponto_aux, lat_cal = lat_cal_aux, lon_cal = lon_cal_aux, id_way=id_way_aux ,semelhanza = semelhanza_aux, road_true=road_true_aux, off_road=off_road_aux)
    rota_true<- rbind(rota_true,rota_true1)
    rm(list=c('rota_true1','lat_ponto_aux','lon_ponto_aux','id_way_aux','lat_cal_aux','lon_cal_aux','semelhanza_aux','road_true_aux','off_road_aux'))
    
  }
  return(rota_true)
}