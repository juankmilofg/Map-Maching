#Universidade Federal de Minas Gerais - Brasil
#LITC - Computational Intelligence and Technology Laboratory. http://litc.cpdee.ufmg.br.
#Code developed by:
#Jose Maia Neto - jmnt@ufmg.br, Juan Camilo Fonseca Galindo - juankmilofg@gmail.com

#functions

DTW.dist <- function(lat_ref,lon_ref, traj){
  
  library(dtw)

  lon_ref <- as.numeric(lon_ref)
  lat_ref <- as.numeric(lat_ref)
  TrajRef <- t(rbind(lat_ref,lon_ref))
  DTWdist <- dtw(TrajRef, as.matrix(traj[,3:4]))$distance
  
  return(DTWdist)
  
}

#return distance coordinates (m)
geodist <- function(from_lat, from_lon, to_lat, to_lon, units="m")
{
  units <- match.arg(units, c("m","nm"))
  
  rad <- 180 / pi
  
  N1 <- from_lat / rad
  E1 <- from_lon / rad
  N2 <- to_lat   / rad
  E2 <- to_lon   / rad
  
  duplicates <- N1==N2 & E1==E2
  N1[duplicates] <- 0            # When origin and destination are the same,
  E1[duplicates] <- 0            #   set them both to 0, 0
  N2[duplicates] <- 0            # Without this, geodist(48.535, 124, 48.535, 124) returns NaN,
  E2[duplicates] <- 0            #   but geodist(0, 0, 0, 0) seems to return 0 on all machines
  
  radians <- acos(sin(N1)*sin(N2)+cos(N1)*cos(N2)*cos(E1-E2))
  if (is.na(radians)) {radians<-0}
  
  distance <- if(units=="m") 60*rad*radians*1.852 else 60*rad*radians
  
  return(distance*1000)
}

#convertor de degree para radianes
deg2rad <- function(deg) {(deg * pi) / (180)}
#convertor de radianes para degree 
rad2deg <- function(deg) {(deg * 180) / (pi)}

#retorna o heading angle de dois pontos
heading_angle <- function(latitute_a, longitude_a,latitute_b , longitude_b)
{
  
  lon1=deg2rad(longitude_a)
  lat1=deg2rad(latitute_a)
  lon2=deg2rad(longitude_b)
  lat2=deg2rad(latitute_b)
  
  X <- cos(lat2)*sin(lon2-lon1)
  Y <- cos(lat1)*sin(lat2) - sin(lat1)*cos(lat2)*cos(lon2-lon1)
  angle <- rad2deg ( atan2(X,Y) )
  
  return(angle)
}

delete_frame <- function(frame1)
{
  frame2<-geodf[c(1:length(frame1$lat)-1),]
  return(frame2)
}

#limita o espaço de busqueda, forma um recuadro  de lados (2 veces o maior comprimento do mapa) centrado no nó
box_way_lim <- function(lat_node,lon_node,way,nodes,geodf,plot_box)
{
  way_max<-max(way$weights) #maxíma distancia do mapa
  
  muc_bbox <- center_bbox(lon_node, lat_node, 2*way_max, 2*way_max)
  
  lat_nova_sup <- as.numeric(muc_bbox[4]) 
  lat_nova_inf <- as.numeric(muc_bbox[2])
  lon_nova_sup <- as.numeric(muc_bbox[3]) 
  lon_nova_inf <- as.numeric(muc_bbox[1])
  
  if (  as.numeric(min(nodes$lat)) >  lat_nova_inf ) 
  {lat_nova_inf <- min(nodes$lat)}
  if (  as.numeric(max(nodes$lat)) <  lat_nova_sup ) 
  {lat_nova_sup <- max(nodes$lat)}
  
  if (  as.numeric(min(nodes$lon)) >  lon_nova_inf ) 
  {lon_nova_inf <- min(nodes$lon)}
  if (  as.numeric(max(nodes$lon)) <  lon_nova_sup ) 
  {lon_nova_sup <- max(nodes$lon)}
  
  x=c(lat_nova_sup,lat_nova_sup,lat_nova_inf,lat_nova_inf,lat_nova_sup)
  y=c(lon_nova_inf,lon_nova_sup,lon_nova_sup,lon_nova_inf,lon_nova_inf)
  
  if (plot_box){
    lines(y,x, type = 'l', col = "blue", lwd = 2)} #graficar recuadro no mapa
  
  nodes.aux <- subset(nodes, lat <= lat_nova_sup & lat >= lat_nova_inf & lon <= lon_nova_sup & lon >= lon_nova_inf)
  aux <-subset(way, from_node_id == nodes.aux$id[1])
  for(i in c(2:length(nodes.aux$id))){
    aux <- rbind (aux,subset(way, from_node_id == nodes.aux$id[i]))
  }
  for(i in c(2:length(nodes.aux$id))){
    aux <- rbind (aux,subset(way, to_node_id == nodes.aux$id[i]))
  }
  return(aux)
}

#Some function to shift vectors(will be used later)
shift.vec <- function (vec, shift) {
  if(length(vec) <= abs(shift)) {
    rep(NA ,length(vec))
  }else{
    if (shift >= 0) {
      c(rep(NA, shift), vec[1:(length(vec)-shift)]) }
    else {
      c(vec[(abs(shift)+1):length(vec)], rep(NA, abs(shift))) } } }

#Distância minima entre um caminho e um ponto, o caminho é dado pela união de dois nós
#entradas:
# Nó 1: (x_1,y_1) (lat,lon)
# Nó 2: (x_2,y_2)
# Ponto: (x,y)
#Saida: 
# distância minima
DistMin_way_point <- function(x_1,y_1,x_2,y_2,lat,lon)
{
  if(x_2 == x_1)
  {
    x_0 <- x_1       # ponto de interseção da perpendicular e o way (x_0,y_0)
    y_0 <- lon
  }else if (y_2 == y_1)
  {
    x_0 <- lat       # ponto de interseção da perpendicular e o way (x_0,y_0)
    y_0 <- y_1
  }else
  {
    m_1 <- (y_2-y_1)/(x_2-x_1)                     # pendiente way
    m_2 <- (-1)/m_1                                # pendiente reta perpendicular entre o ponto e o way
    x_0 <- (m_2*lat-m_1*x_1+y_1-lon)/(m_2-m_1)         # ponto de interseção da perpendicular e o way (x_0,y_0)
    y_0 <- m_1*(x_0-x_1)+y_1
  }
  distamin <- 0
  #if (  (( (x_0<=x_1)&(x_0>=x_2) )|( (x_0>=x_1)&(x_0<=x_2 ))) &  (( (y_0<=y_1)&(y_0>=y_2) )|( (y_0>=x_1)&(y_0<=y_2 )))  ){
  if (  ((x_0<=x_1)&(x_0>=x_2) )|( (x_0>=x_1)&(x_0<=x_2 ) ) |  ( (y_0<=y_1)&(y_0>=y_2) )|( (y_0>=x_1)&(y_0<=y_2 ))  ){
    distamin <- geodist(x_0,y_0,lat,lon)
    a<-1
  }  else{
    dist1<-geodist(x_1,y_1,lat,lon)
    dist2<-geodist(x_2,y_2,lat,lon)
    if (dist1<dist2)  {  distamin <- dist1  
    }else             {  distamin <- dist2  }
  }
  return(c(distamin,x_0,y_0))
}

#Ângulo entre um caminho e um ponto, o caminho é dado pela união de dois nós
#entradas:
# Nó from: (x_1,y_1) (lat,lon)
# Nó to: (x_2,y_2)
# heading Ponto: head
#Saida: 
# ângulo entre o way e o ponto 
heading_way_point <- function(x_1,y_1,x_2,y_2,head_point)
{
  head_way <- heading_angle(x_1,y_1,x_2,y_2)
  diff_head<- head_way-head_point
  return(diff_head)
}

#Ângulo entre um caminho e um ponto, o caminho é dado pela união de dois nós
#entradas:
# Nó from: (x_1,y_1) (lat,lon)
# Nó to: (x_2,y_2)
# heading Ponto: head
#Saida: 
# ângulo entre o way e o ponto (abs)
heading_way_point_abs <- function(x_1,y_1,x_2,y_2,head_point)
{
  head_way <- heading_angle(x_1,y_1,x_2,y_2)
  diff_head<- abs(head_way)-abs(head_point)
  return(diff_head)
}

#Dcc (Direction Chain Code)
#entradas:
# Nó from: (x_1,y_1) (lat,lon)
# Nó to: (x_2,y_2)
# heading Ponto: head
#Saida: 
# dcc
Dcc <- function(x_1,y_1,x_2,y_2,head_point)
{
  angle <- heading_way_point(x_1,y_1,x_2,y_2,head_point)
  
  angle <- 180
  
  if(       (angle>=-22.5)  & (angle<=22.5))    {Delta<-0
  }else if( (angle>22.5)    & (angle<=67.5))    {Delta<-1
  }else if( (angle>67.5)    & (angle<=112.5))   {Delta<-2
  }else if( (angle>112.5)   & (angle<=157.5))   {Delta<-3
  }else if( (angle< -112.5) & (angle>= -157.5)) {Delta<-7
  }else if( (angle< -67.5)  & (angle>= -112.5)) {Delta<-6
  }else if( (angle< -22.5)  & (angle>= -67.5))  {Delta<-5
  }else                                         {Delta<-4}
  
  if (Delta<=4) { dcc<-Delta
  }else         { dcc<- (Delta - 8) %% 4}
  
  return(dcc)
}

#lista de candidatos
list_candidatos <- function(lat_ponto,lon_ponto,box_way_limites,nodes.map,ErroGPS)
{
  
  count_candidatos <- FALSE
  way_candidates <- NULL
  for (i in c(1:length(box_way_limites[all(TRUE),1])))
  {
    node_from_aux=subset(nodes.map, id == box_way_limites$from_node_id[i])
    node_to_aux=subset(nodes.map, id == box_way_limites$to_node_id[i])
    dist_min_node <- DistMin_way_point(node_from_aux$lat[1],node_from_aux$lon[1],node_to_aux$lat[1],node_to_aux$lon[1],lat_ponto,lon_ponto)
    
    if (dist_min_node[1]<=ErroGPS*2)
    {
      from_node_id_aux=box_way_limites$from_node_id[i]
      to_node_id_aux=box_way_limites$to_node_id[i]
      way_id_aux=box_way_limites$way_id[i]
      weights_aux=box_way_limites$weights[i]
      way_candidates1 <- data.frame(from_node_id=from_node_id_aux, to_node_id=to_node_id_aux, way_id=way_id_aux, weights=weights_aux)
      way_candidates <- rbind(way_candidates,way_candidates1)  
      count_candidatos <- TRUE
    }
  }
  return(way_candidates)
}

#Funções
#A função as_igraph2 transforma um objeto .osm do mapa digital de uma região 
#(digital road network map) em uma tabela que representa o modelo da rede de 
#nós e ruas. A saída desta função é um data frame contento as seguintes informações
#sobre todos os segmentos de ruas que estão na região do mapa da qual o objeto
#.osm (obj) diz respeito:
#1a coluna : from_node_id -> representa o id do nó de partida
#2a coluna : to_node_id  -> representa o id do nó de destino
#3a coluna : way_id -> representa o id do rua na qual o segemento (from_node_id -> to_node_id) está contido.
#4a coluna : weights -> representa o comprimento deste segmento (a distância entre os nós de partida e chegada) baseado na distância Haversine.(Modelo esférico da Terra e raio da Terra igual a 6378137 metros) em metros.

as_igraph2 <- function (obj) 
{
  stopifnot(is_osmar(obj))
  stopifnot(require("igraph"))
  dat <- merge_ways_nodes(obj$ways[[3]], obj$nodes[[1]])
  dat <- split(dat, dat$id)
  dat <- dat[sapply(dat, nrow) >= 2]
  edges <- lapply(dat, function(x) {
    n <- nrow(x)
    from <- 1:(n - 1)
    to <- 2:n
    #weightsold <- distHaversine(x[from, c("lon", "lat")], x[to, c("lon", "lat")])
    weights <- geodist( x[from, c("lat")],x[from, c("lon")], x[to, c("lat")],x[to, c("lon")])
    cbind(from_node_id = x[from, "ref"], to_node_id = x[to,"ref"], way_id = x[1, "id"], weights = weights)#,weightsold = weightsold)
  })
  
  edges <- do.call(rbind, edges)
  edges1 <- edges
  weights <- edges[, "weights"]
  #weightsold <- edges[, "weightsold"]
  names <- edges[, "way_id"]
  edges <- cbind(as.character(edges[, "from_node_id"]), as.character(edges[,"to_node_id"]))
  graph <- graph.edgelist(edges)
  E(graph)$weight <- weights
  #E(graph)$weightsold <- weightsold
  E(graph)$name <- names
  graph
  edges1
}

merge_ways_nodes <- function(ways, nodes) {
  colnames(ways) <- sprintf("w%s", colnames(ways))
  colnames(nodes) <- sprintf("n%s", colnames(nodes))
  
  m <- match(ways$wref, nodes$nid)
  
  dat <- cbind(ways, nodes[m, ])
  # dat <- na.omit(dat)
  dat <- dat[!is.na(dat$nlat), ]
  
  dat$nid <- NULL
  colnames(dat) <- substring(colnames(dat), 2)
  
  dat
}

is_osmar <- function(obj) {
  OSMAR_CLASS %in% class(obj)
}

NODES_CLASS <- "nodes"
WAYS_CLASS <- "ways"
RELATIONS_CLASS <- "relations"
OSMAR_CLASS <- "osmar"
osmar_class <- function(obj) {
  stopifnot(length(obj) == 3)
  #stopifnot(sapply(obj,
  # function(k) class(k)[1])==c("NODE", "WAY", "RELATION"))
  subclass(obj, OSMAR_CLASS)
}

#teste de distância entre a rota calculada e a rota referencia.
#Entrada:
#lat_ref = longitudes da rota de referencia
#lon_ref = latitude da rota de referencia
#rota_true = rota calculada
#saida
#distancia_min = vetor de da distância minima entre a rota calculada e a roa referencia
test_dist <- function(lat_ref, lon_ref, rota_true,param)
{
  #param  0 media, 1 mediana ,2 vetor de distancia ,3 sumatoria das distancias, 4 sum (exp(dismin))
  distancia_min<-NULL
  for(ponto in c(1:length(rota_true$lat_cal)))
  {
    distancias<-NULL
    for(candidato in c(1:(length(lat_ref)-1)))
    {
      lat1 <- as.numeric(lat_ref[candidato])
      lon1 <- as.numeric(lon_ref[candidato])
      lat2 <- as.numeric(lat_ref[candidato+1])
      lon2 <- as.numeric(lon_ref[candidato+1])
      lat3 <- as.numeric(rota_true$lat_cal[ponto])
      lon3 <- as.numeric(rota_true$lon_cal[ponto])
      distancias[candidato]<-DistMin_way_point(lat1,lon1,lat2,lon2,lat3,lon3)[1]
    }
    distancia_min[ponto]<-min(distancias)
  }
  if (param==1){distancia_total<-median(distancia_min)
  }else if(param==3){distancia_total<-sum(distancia_min)
  }else if(param==4){distancia_total<-sum(1/exp(distancia_min))
  }else if(param==0){distancia_total<-mean(distancia_min)
  }else{distancia_total<-distancia_min}
    
  return(distancia_total)
}
test_road <- function(WaysID, rota_true)
{
  roads_true<-0
  for(ponto in c(1:length(rota_true$lat_cal)))
  {
    match_true <- match(rota_true$id_way[ponto],WaysID)
    if (is.na(match_true))
    {
    }else
    {roads_true <- roads_true+1}
    road_true <- roads_true / length(rota_true$lat_cal)
  }
  return(road_true)
}

test_dist_road = function(lat_ref, lon_ref, IDway_ref, rota_true)
{
  distancia_min<-NULL
  roads_true<-0
  for(ponto in c(1:length(rota_true$lat_cal)))
  {
    distancias<-NULL
    for(candidato in c(1:(length(lat_ref)-1)))
    {
      lat1 <- as.numeric(lat_ref[candidato])
      lon1 <- as.numeric(lon_ref[candidato])
      lat2 <- as.numeric(lat_ref[candidato+1])
      lon2 <- as.numeric(lon_ref[candidato+1])
      lat3 <- as.numeric(rota_true$lat_cal[ponto])
      lon3 <- as.numeric(rota_true$lon_cal[ponto])
      distancias[candidato]<-DistMin_way_point(lat1,lon1,lat2,lon2,lat3,lon3)[1]
    }
    distancia_min[ponto]<-min(distancias)[1]
    way_aux <- which(distancias==distancia_min[ponto])
    
    id_way_rota <- rota_true$id_way[ponto]
    id_way_ref <- IDway_ref[way_aux]  
    if( !is.na(id_way_rota) & !is.na(id_way_ref) ) {
      if (id_way_rota==id_way_ref){roads_true<-roads_true +1}
    }
  }
  return(roads_true/length(rota_true$lat_cal))
}


way_ref <- function(lat_ref,lon_ref,nodes.map,way.map,geodf)
{
  lat_ref_from=as.numeric(lat_ref[1])
  lon_ref_from=as.numeric(lon_ref[1])
  
  dist.1 <- NULL
  dist.match1 <- NULL
  dist.match1[1]<-Inf
  
  for (di in c(1:length(nodes.map$id)))
  {
    #di<-1
    dist.1[di] <- geodist(lat_ref_from,lon_ref_from,nodes.map$lat[di],nodes.map$lon[di])
  }
  dist_aux.1 <- which(dist.1==min(dist.1))[1]
  #print(min(dist.1))
  
  nodes_match  <- data.frame(id=nodes.map$id[dist_aux.1], lat=nodes.map$lat[dist_aux.1], lon=nodes.map$lon[dist_aux.1],way=NA)
  count=1
  
  for (point in c(1:(length(lat_ref)-1) ))
  {
    lat_ref_from=as.numeric(lat_ref[point])
    lon_ref_from=as.numeric(lon_ref[point])
    lat_ref_to=as.numeric(lat_ref[point+1])
    lon_ref_to=as.numeric(lon_ref[point+1])
    count2<-1
    
    if( (lat_ref_to>min(geodf$lat)) & (lat_ref_to<max(geodf$lat)) & (lon_ref_to>min(geodf$lon)) & (lon_ref_to<max(geodf$lon))  )
    {
      end <- FALSE
      #print("entro")
      while(end ==FALSE) {
        
        aux_to.a <- subset(way.map, from_node_id == nodes_match$id[count])
        if (length(aux_to.a$way_id)==0){dist.a<- Inf
        }else{
          dist.a<-NULL
          for (di in c(1:(length(aux_to.a$way_id))))
          {
            node.a <- subset(nodes.map, id == aux_to.a$to_node_id[di])
            if (count==1){
              lat_to <- node.a$lat[1]
              lon_to <- node.a$lon[1]
              dist.a[di] <- DistMin_way_point(lat_ref_from,lon_ref_from,lat_ref_to,lon_ref_to,lat_to,lon_to)
            }else{
              if (node.a$id[1]==nodes_match$id[count-1]){dist.a[di] <- Inf
              }else{
                lat_to <- node.a$lat[1]
                lon_to <- node.a$lon[1]
                dist.a[di] <- DistMin_way_point(lat_ref_from,lon_ref_from,lat_ref_to,lon_ref_to,lat_to,lon_to)
              }
            }
          }
        }
        dist_aux.a <- which(dist.a==min(dist.a))[1]
        id.a <- aux_to.a$to_node_id[dist_aux.a]
        
        aux_to.b <- subset(way.map, to_node_id   == nodes_match$id[count])
        if (length(aux_to.b$way_id)==0){dist.b<- Inf
        }else{
          dist.b<-NULL
          for (di in c(1:(length(aux_to.b$way_id))))
          {
            node.b <- subset(nodes.map, id == aux_to.b$from_node_id[di])
            if (count==1){
              lat_to <- node.b$lat[1]
              lon_to <- node.b$lon[1]
              dist.b[di] <- DistMin_way_point(lat_ref_from,lon_ref_from,lat_ref_to,lon_ref_to,lat_to,lon_to)
            }else{
              if (node.b$id[1]==nodes_match$id[count-1]){dist.b[di] <- Inf
              }else{
                lat_to <- node.b$lat[1]
                lon_to <- node.b$lon[1]
                dist.b[di] <- DistMin_way_point(lat_ref_from,lon_ref_from,lat_ref_to,lon_ref_to,lat_to,lon_to)
              }
            }
          }
        }
        dist_aux.b <- which(dist.b==min(dist.b))[1]
        id.b <- aux_to.b$from_node_id[dist_aux.b]
        
        if ((min(dist.a)) <= (min(dist.b)) ) {id_to=id.a
        }else{id_to=id.b}
        
        if(!((min(dist.a)==Inf) & (min(dist.b)==Inf))){
          node_to <- subset(nodes.map, id == id_to)
          
          way_1 <- subset(way.map, from_node_id == nodes_match$id[count] & to_node_id == node_to$id[1])
          if (length(way_1$way_id)==0) {
            way_1 <- subset(way.map, to_node_id == nodes_match$id[count] & from_node_id == node_to$id[1])
          }
          
          nodes_match_aux  <- data.frame(id=node_to$id[1], lat=node_to$lat[1], lon=node_to$lon[1] ,way=way_1$way_id)
          
          nodes_match  <- rbind(nodes_match,nodes_match_aux)
          
          count <- count +1
          count2 <- count2 +1
          
          dist.match <- geodist(nodes_match$lat[count],nodes_match$lon[count],lat_ref_to,lon_ref_to)
          #print(dist.match)
          
          dist.match1[count2]<-dist.match
          if (dist.match1[count2]>dist.match1[count2-1]){end<-TRUE}
          if (dist.match<1){end<-TRUE}
        }else{end<-TRUE}
        
      }
    }
  }
  return(nodes_match)
  
}   
