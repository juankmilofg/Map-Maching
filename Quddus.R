#Universidade Federal de Minas Gerais - Brasil
#LITC - Computational Intelligence and Technology Laboratory. http://litc.cpdee.ufmg.br.
#Code developed by:
#Jose Maia Neto - jmnt@ufmg.br, Juan Camilo Fonseca Galindo - juankmilofg@gmail.com
#Code based on article:
# Quddus, M. A., Noland, R. B. and Ochieng, W. Y. (2006). 
# A high accuracy fuzzy logic based map matching algorithm for road transport, 
# Journal of Intelligent Transportation Systems 10(3): 103–115
# The developed code uses parts of the work done by of the work done by Nikolai Gorte, http://github.com/ngort01/fuzzyMM

#### METHOD Quddus 2006 ####

fuzzyMMfct <- function(traj){
  
  require(rgeos)
  require(frbs)
  require(osmar)
  require(rjson)
  require(RCurl)
  require(maptools)
  require(stringr)
  require(rgdal)
  require(lubridate)
  
  #Fun??es auxiliares
  create_drn <- function(bbox) {
    x1 <- bbox[[1]] 
    y1 <- bbox[[2]]
    x2 <- bbox[[3]] 
    y2 <- bbox[[4]] 
    
    if (!requireNamespace("RCurl", quietly = TRUE))
      stop("package RCurl required")
    # Using overpass API, because it offers better possibilties for 
    # filtering data than the the OSM API
    url <- paste0("http://www.overpass-api.de/api/xapi?way[bbox=",x1,",",y1,",",x2,",",y2,"][highway=*]")  
    response <- RCurl::getURL(url, .encoding = "UTF-8")
    
    if (!requireNamespace("XML", quietly = TRUE))
      stop("package XML required")
    # Parse Data
    resp <- XML::xmlParse(response)
    
    # Transform parsed data to osmar object
    roads <- as_osmar(resp)
    
    v <- k <- NULL # make visible bindings and R CMD check happy
    
    # Get ID's of all streets used by cars
    id <- find(roads, way(tags(k == "highway" & !(v %in% c("cycleway", 
                                                           "footway", "bridleway", "steps", "path")))))
    roads <- subset(roads, ids = find_down(roads, way(id)))
    
    # Get coordinates of each node
    nodes <- roads$nodes
    coords <- nodes$attrs[c("id", "lat", "lon")]
    
    # Create an igraph from the osmar object
    graph <- as_igraph(roads)
    graph <- as.undirected(graph, mode = "each")
    V(graph)$id <- as.numeric(V(graph)$name)
    # Bring coordinates in the right order
    coords <- coords[match(V(graph)$id, coords$id),]
    V(graph)$lon <- coords$lon
    V(graph)$lat <- coords$lat
    
    # Convert the osmar object to spatial lines and split each line into segments
    roads <- as_sp(roads, "lines")
    roads <- lines2segments(roads)
    #roads <- SpatialLinesDataFrame(roads, as.data.frame(get.edgelist(graph)), match.ID = FALSE)
    roads <- new("DigitalRoadNetwork", sl = roads, g = graph)
    roads
  }
  
  setClass("igraph")
  setClass("DigitalRoadNetwork", representation(sl = "SpatialLinesDataFrame", g = "igraph"), 
           validity = function(object) {stopifnot(length(object@sl) == length(E(object@g)))})
  
  lines2segments <- function(sl){
    coords <- coordinates(sl)
    osm_ids <- sl@data$id # osm edge ids
    in_nrows <- lapply(coords, function(x) sapply(x, nrow))
    outn <- sapply(in_nrows, function(y) sum(y-1))
    osm_ids <- rep(osm_ids, outn)
    res <- vector(mode = "list", length = sum(outn))
    i <- 1
    for (j in seq(along = coords)) {
      for (k in seq(along = coords[[j]])) {
        for (l in 1:(nrow(coords[[j]][[k]]) - 1)) {
          res[[i]] <- coords[[j]][[k]][l:(l + 1),]
          i <- i + 1
        }
      }
    }
    res1 <- vector(mode = "list", length = sum(outn))
    for (i in seq(along = res))
      res1[[i]] <- Lines(list(Line(res[[i]])), as.character(i))
    outSL <- SpatialLines(res1, osm_crs())
    outSL <- SpatialLinesDataFrame(outSL, data.frame(osm_ids))
    outSL
  }
  
  
  ## Cache environment for package variables
  cacheEnv <- new.env()
  
  ## Internal function used to create the data.frame containing the bounds of the fuzzy subsets.
  ##
  ## Returns:
  ##  data.frame containing the the default range of the fuzzy subsets 
  init_vars <- function(){
    left_bounds <- c(3, 2, 0, 20, 25, 10, 20, 3, 4, 85, 90, 85, 90, -5, -5, 10, 15, 150, 0, 0, 5, 10)
    right_bounds <- c(6, 4, 2, 45, 60, 40, 50, 5, 6, 100, 120, 100, 120, 5, 10, 20, 30, 200, 1, 1, 15, 25)
    df <- data.frame(left_bounds, right_bounds, row.names = c("speed_high", "speed_low", "speed_zero", "HE_small",
                                                              "HE_large", "PD_short", "PD_long", "HDOP_good",
                                                              "HDOP_bad", "alpha_low", "alpha_high", "beta_low", "beta_high", 
                                                              "delta_dist_neg", "delta_dist_pos", 
                                                              "HI_small", "HI_large", "HI_180", "connectivity_direct", 
                                                              "connectivity_indirect", "dist_err_small", "dist_err_large"))
    df$ID <- 1:nrow(df)
    df
  }
  
  assign("var_bounds", init_vars(), envir = cacheEnv)
  
  
  
  #' Fuzzy subset range
 
  get_var_bounds <- function() get("var_bounds", envir = cacheEnv)
  
  
  
  #' Set bounds of fuzzy subsets
  set_var_bounds <- function(name = c("speed_high", "speed_low", "speed_zero", "HE_small",
                                      "HE_large", "PD_short", "PD_long", "HDOP_good",
                                      "HDOP_bad", "alpha_low, alpha_high", "beta_low", 
                                      "beta_high", 
                                      "delta_dist_neg", "delta_dist_pos", "HI_small", 
                                      "HI_large", "HI_180", "connectivity_direct", 
                                      "connectivity_indirect", "dist_err_small", "dist_err_large"),
                             bounds = "numeric", default = FALSE) {
    if(default) {
      assign("var_bounds", init_vars(), envir = cacheEnv) 
    } else {
      name <- match.arg(name)
      if (is.null(bounds))
        stop ("No bounds specified!")
      if (!is(bounds, "numeric"))
        stop ("Bounds must be numeric!")
      if (!length(bounds) == 2)
        stop ("Bound must be of length 2!")
      
      var_bounds <- get_var_bounds()
      var_bounds[rownames(var_bounds) == name, 1] <- bounds[1]
      var_bounds[rownames(var_bounds) == name, 2] <- bounds[2]
      assign("var_bounds", var_bounds, envir = cacheEnv) 
    }
  }
  
  
  
  #' Update Membership Functions
  update_mf <- function() {
    create_fis1()
    create_fis2()
    create_fis3()
  } 
  
  
  #' Get Fuzzy Inference System
  #' Get the Fuzzy Inference System for IMP, SMP1 or SMP2.
  get_fis <- function(name = c("IMP", "SMP1", "SMP2")) {
    name <- match.arg(name)
    if (name == "IMP")
      get("fis1", envir = cacheEnv)
    else if(name == "SMP1")
      get("fis2", envir = cacheEnv)
    else if (name == "SMP2") 
      get("fis3", envir = cacheEnv)
  }
  
  ## Set model: Takagi Sugeno Kang (TSK).
  type.model <- "TSK"
  
  ## Use standard t-norm and s-norm.
  type.tnorm <- "MIN"
  type.snorm <- "MAX"
  
  ## Type of implication function
  type.implication.func <- "MIN"
  
  ## Gets the coefficients for the sigmoidal membership functions.
  ##
  ## Args:
  ##  left: x values at which the sigmoidal membership function reachs ~0
  ##  right: x values at which the sigmoidal membership function reachs ~1
  ##  shape: shape of the sigmoidal membership function. s = ascending, z = descending
  ##
  ## Returns:
  ##  slope of the function at the crossover point
  get_params <- function(l, r, shape = c("s", "z")) {
    shape <- match.arg(shape)
    if (shape == "s") 
      y <- c(0.01, 0.5, 0.99)
    else
      y <- c(0.99, 0.5, 0.01)
    x <- c(l, (l + r)/2, r)
    slope <- ifelse(shape == "s", 1/(r - l),  1/(r - l))
    data <- list(x = x, y = y)
    fitModel <- nls(y ~ a/(1 + exp(-b * (x - ((l + r)/2)))), 
                    data, start = c(a = 1, b = slope), 
                    algorithm = "port")  
    # get the coefficients using the coef function
    params <- coef(fitModel)
    params[2]
  }
  
  ## gets parameter c 
  get_mid <- function(bounds, row) {
    (bounds[row, 1] + bounds[row, 2])/2
  }
  
  
  #' Fuzzy Inference System 1 (FIS1)
  #' 
  #' Fuzzy Inference System used in the Initial Map Matching Process (IMP).
  ## Function that creates FIS1 used in IMP
  ##
  ## Returns:
  ##  frbs object
  create_fis1 <- function() {
    
    # Types used: 6 = sigmoidal
    var_bounds <- get_var_bounds()
    m <- matrix(c(6, get_params(var_bounds[1, 1], var_bounds[1, 2], "s"), #speed_high
                  get_mid(var_bounds, 1), NA, NA,   
                  6, get_params(var_bounds[2, 1], var_bounds[2, 2], "z"), #sped_low 
                  get_mid(var_bounds, 2), NA, NA, 
                  6, get_params(var_bounds[3, 1], var_bounds[3, 2], "z"), #speed_zero
                  get_mid(var_bounds, 3), NA, NA, 
                  6, get_params(var_bounds[4, 1], var_bounds[4, 2], "z"), #HE_small 
                  get_mid(var_bounds, 4), NA, NA,  
                  6, get_params(var_bounds[5, 1], var_bounds[5, 2], "s"), #HE_large
                  get_mid(var_bounds, 5), NA, NA, 
                  6, get_params(var_bounds[6, 1], var_bounds[6, 2], "z"), #PD_short
                  get_mid(var_bounds, 6), NA, NA, 
                  6, get_params(var_bounds[7, 1], var_bounds[7, 2], "s"), #PD_long
                  get_mid(var_bounds, 7), NA, NA, 
                  6, get_params(var_bounds[8, 1], var_bounds[8, 2], "z"), #HDOP_good
                  get_mid(var_bounds, 8), NA, NA, 
                  6, get_params(var_bounds[9, 1], var_bounds[9, 2], "s"), #HDOP_bad
                  get_mid(var_bounds, 9), NA, NA), 
                nrow = 5, byrow = FALSE)
    assign("varinp.mf1", m, envir = cacheEnv)
    
    
    
    
    
    # Names of the variables
    colnames.var1 <- c("v", "HE", "PD", "HDOP", "output")
    
    # Give the names of the fuzzy terms of each input variable
    # Names of the fuzzy terms must be unique
    varinput.1 <- c("high", "low", "zero")
    varinput.2 <- c("small", "large")
    varinput.3 <- c("short", "long")
    varinput.4 <- c("good", "bad")
    names.varinput1 <- c(varinput.1, varinput.2, varinput.3, varinput.4)
    
    # Set interval of data. "v", "HE", "PD", "HDOP", "output"
    range.data1 <- matrix(c(0, 50, 0, 360, 0, 60, 0, 20, 0, 100), nrow = 2)
    
    # Define number of fuzzy terms of input variables.
    num.fvalinput1 <- matrix(c(3, 2, 2, 2), nrow = 1)
    
    
    # Set the name of the simulation
    name1 <- "Sim-1"
    
    
    # Define the fuzzy IF-THEN rules 
    # For TSK model, it isn't necessary to put linguistic terms in consequent parts
    r1 <- c("high","and","small","and","dont_care","and","dont_care", "->")
    r2 <- c("high","and","large", "and", "dont_care", "and", "dont_care", "->")
    r3 <- c("dont_care", "and", "dont_care", "and", "short","and","good", "->")
    r4 <- c("dont_care", "and", "dont_care", "and", "long","and","good", "->")
    r5 <- c("dont_care", "and", "small","and","short", "and", "dont_care", "->")
    r6 <- c("dont_care", "and", "large","and","long", "and", "dont_care", "->")
    rule1 <- list(r1, r2, r3, r4, r5, r6)
    rule1 <- do.call(rbind, rule1)
    
    # Define linear functions of TSK 
    func.tsk1 <- matrix(c(50, 10, 50, 10, 100, 10), nrow = 6, byrow = TRUE)
    
    
    
    
    # Generate a fuzzy model with frbs.gen
    # For TSK model we do not need to input: 
    # num.fvaloutput, varout.mf, names.varoutput, type.defuz
    varinp.mf1 <- get("varinp.mf1", envir = cacheEnv)
    fis1 <- frbs.gen(range.data1, num.fvalinput1, names.varinput1, num.fvaloutput = NULL, 
                     varout.mf = NULL, names.varoutput = NULL, rule1, 
                     varinp.mf1, type.model, type.defuz = NULL, type.tnorm, type.snorm, func.tsk1, colnames.var1, 
                     type.implication.func, name1)
    
    assign("fis1", fis1, envir = cacheEnv)
  }
  
  create_fis1()
  
  
  #' Fuzzy Inference System 2 (FIS2)
  #' 
  #' Fuzzy Inference System used in the Subsequent Map Matching Process
  #' along a link (SMP-1).
  ## Function that creates FIS2 used in SMP-1
  ##
  ## Returns:
  ##  frbs object
  create_fis2 <- function() {
    # Types used: 6 = sigmoidal
    var_bounds <- get_var_bounds()
    m <- matrix(c(6, get_params(var_bounds[1, 1], var_bounds[1, 2], "s"), #speed_high
                  get_mid(var_bounds, 1), NA, NA,   
                  6, get_params(var_bounds[2, 1], var_bounds[2, 2], "z"), #speed_low 
                  get_mid(var_bounds, 2), NA, NA, 
                  6, get_params(var_bounds[3, 1], var_bounds[3, 2], "z"), #speed_zero
                  get_mid(var_bounds, 3), NA, NA, 
                  6, get_params(var_bounds[8, 1], var_bounds[8, 2], "z"), #HDOP_good
                  get_mid(var_bounds, 8), NA, NA, 
                  6, get_params(var_bounds[9, 1], var_bounds[9, 2], "s"), #HDOP_bad
                  get_mid(var_bounds, 9), NA, NA,
                  6, get_params(var_bounds[10, 1], var_bounds[10, 2], "z"), #alpha_below90
                  get_mid(var_bounds, 10), NA, NA,   
                  6, get_params(var_bounds[11, 1], var_bounds[11, 2], "s"), #alpha_above90
                  get_mid(var_bounds, 11), NA, NA, 
                  6, get_params(var_bounds[12, 1], var_bounds[12, 2], "z"), #beta_below90
                  get_mid(var_bounds, 12), NA, NA,   
                  6, get_params(var_bounds[13, 1], var_bounds[13, 2], "s"), #beta_above90
                  get_mid(var_bounds, 13), NA, NA,  
                  6, get_params(var_bounds[14, 1], var_bounds[14, 2], "z"), #delta_dist_neg
                  get_mid(var_bounds, 14), NA, NA, 
                  6, get_params(var_bounds[15, 1], var_bounds[15, 2], "s"), #delta_dist_pos
                  get_mid(var_bounds, 15), NA, NA, 
                  6, get_params(var_bounds[16, 1], var_bounds[16, 2], "z"), #HI_small
                  get_mid(var_bounds, 16), NA, NA, 
                  6, get_params(var_bounds[17, 1], var_bounds[17, 2], "s"), #HI_large
                  get_mid(var_bounds, 17), NA, NA, 
                  4, 125, 150, 200, 225),#HI_180
                nrow = 5, byrow = FALSE)
    
    assign("varinp.mf2", m, envir = cacheEnv)
    
    
    # Names of the variables
    colnames.var2 <- c("speed", "HDOP", "alpha", "beta", "delta_dist", 
                       "HI", "HI1", "output")
    
    # Give the names of the fuzzy terms of each input variable
    # Names of the fuzzy terms must be unique
    varinput2.1 <- c("fast", "slow", "zero")
    varinput2.2 <- c("good", "bad")
    varinput2.3 <- c("below90", "above90")
    varinput2.4 <- c("below90b", "above90b")
    varinput2.6 <- c("pos", "neg")
    varinput2.7 <- c("small", "large")
    varinput2.8 <- c("180")
    names.varinput2 <- c(varinput2.1, varinput2.2, varinput2.3, varinput2.4,
                         varinput2.6, varinput2.7, varinput2.8)
    
    # Set interval of data."speed", "HDOP", "alpha", "beta", "delta_dist", 
    # "HI", "HI1", "output"
    range.data2 <- matrix(c(0, 50, 0, 20, 0, 360, 0, 360, -500, 500, 0, 360, 0, 360, 0, 100), nrow = 2)
    
    # Define number of fuzzy terms of input variables
    num.fvalinput2 <- matrix(c(3, 2, 2, 2, 2, 2, 1), nrow = 1)
    
    
    # Define the fuzzy IF-THEN rules; 
    # For TSK model, it isn't necessary to put linguistic terms in consequent parts.
    r1 <- c("dont_care", "and", "dont_care", "and","below90","and","below90b","and",
            "dont_care", "and", "dont_care", "and", "dont_care","->")
    r2 <- c("dont_care", "and","dont_care", "and","above90","and", "dont_care", "and",
            "pos","and", "dont_care", "and","dont_care", "->")
    r3 <- c("dont_care", "and","dont_care", "and","dont_care","and", "above90b", "and",
            "pos","and", "dont_care", "and","dont_care", "->")
    r4 <- c("dont_care", "and","dont_care", "and", "below90", "and","below90b","and",
            "dont_care", "and","small","and","dont_care", "->")
    r5 <- c("dont_care", "and","dont_care", "and", "above90", "and","dont_care","and",
            "pos", "and","small","and","dont_care", "->")
    r6 <- c("dont_care", "and","dont_care", "and", "dont_care", "and","above90b","and",
            "pos", "and","small","and","dont_care", "->")
    r8 <- c("dont_care", "and","dont_care", "and", "below90", "and","below90b","and",
            "dont_care", "and","large","and","dont_care", "->")
    r9 <- c("zero","and", "good", "and","dont_care", "and","dont_care", "and",
            "dont_care", "and","dont_care", "and","dont_care", "->")
    r10 <- c("dont_care", "and","good","and","dont_care", "and","dont_care", "and","neg", 
             "and","dont_care", "and","dont_care", "->")
    r11 <- c("dont_care", "and","good","and","dont_care", "and","dont_care", "and","pos", 
             "and","dont_care", "and","dont_care", "->")
    r12 <- c("fast","and","dont_care", "and","dont_care", "and","dont_care", "and",
             "dont_care", "and","small", "and", "dont_care", "->")
    r13 <- c("fast", "and","good","and","dont_care", "and","dont_care", "and",
             "dont_care", "and","dont_care", "and","180", "->")
    rule2 <- list(r1, r2, r3, r4, r5, r6, r8, r9, r10, r11, r12, r13)
    rule2 <- do.call(rbind, rule2)
    
    
    # Set the name of the simulation
    name2 <- "Sim-2"
    
    
    # Define linear functions of TSK  
    func.tsk2 <- matrix(c(100, 10, 10, 100, 10, 10, 10, 100, 50, 10, 50, 100), 
                        nrow = 12, byrow = TRUE)
    
    
    
    # Generate a fuzzy model with frbs.gen
    # For TSK model we do not need to input: 
    # num.fvaloutput, varout.mf, names.varoutput, type.defuz
    varinp.mf2 <- get("varinp.mf2", envir = cacheEnv)
    fis2 <- frbs.gen(range.data2, num.fvalinput2, names.varinput2, num.fvaloutput = NULL, 
                     varout.mf = NULL, names.varoutput = NULL, rule2, 
                     varinp.mf2, type.model, type.defuz = NULL, type.tnorm, type.snorm, func.tsk2, 
                     colnames.var2, type.implication.func, name2)
    
    assign("fis2", fis2, envir = cacheEnv)
  }
  
  create_fis2()
  
  #' Fuzzy Inference System 3 (FIS3)
  #' 
  #' Fuzzy Inference System used in the Subsequent Map Matching Process
  #' at a junction (SMP-2).
  ## Function that creates FIS3 used in SMP-2
  ##
  ## Returns:
  ##  frbs object
  create_fis3 <- function() {
    # Define the shape and parameters of the membership functions
    # Parameters are aquired through the get_params() function
    # Please see fuzzifier function to contruct the matrix
    # First row is the type of membership function
    # Types used: 1 = trinagle, 6 = sigmoidal
    var_bounds <- get_var_bounds()
    m <- matrix(c(6, get_params(var_bounds[1, 1], var_bounds[1, 2], "s"), #speed_high
                  get_mid(var_bounds, 1), NA, NA,   
                  6, get_params(var_bounds[2, 1], var_bounds[2, 2], "z"), #sped_low 
                  get_mid(var_bounds, 2), NA, NA, 
                  6, get_params(var_bounds[3, 1], var_bounds[3, 2], "z"), #speed_zero
                  get_mid(var_bounds, 3), NA, NA, 
                  6, get_params(var_bounds[4, 1], var_bounds[4, 2], "z"), #HE_small 
                  get_mid(var_bounds, 4), NA, NA,  
                  6, get_params(var_bounds[5, 1], var_bounds[5, 2], "s"), #HE_large
                  get_mid(var_bounds, 5), NA, NA, 
                  6, get_params(var_bounds[6, 1], var_bounds[6, 2], "z"), #PD_short
                  get_mid(var_bounds, 6), NA, NA, 
                  6, get_params(var_bounds[7, 1], var_bounds[7, 2], "s"), #PD_long
                  get_mid(var_bounds, 7), NA, NA, 
                  6, get_params(var_bounds[8, 1], var_bounds[8, 2], "z"), #HDOP_good
                  get_mid(var_bounds, 8), NA, NA, 
                  6, get_params(var_bounds[9, 1], var_bounds[9, 2], "s"), #HDOP_bad
                  get_mid(var_bounds, 9), NA, NA,
                  1, -1, 0, 1, NA, #indirect
                  1, 0, 1, 2, NA, #direct
                  6, get_params(var_bounds[21, 1], var_bounds[21, 2], "z"), #dist_err_small
                  get_mid(var_bounds, 21), NA, NA,
                  6, get_params(var_bounds[22, 1], var_bounds[22, 2], "s"), #dist_err_large
                  get_mid(var_bounds, 22), NA, NA),
                nrow = 5, byrow = FALSE)
    
    assign("varinp.mf3", m, envir = cacheEnv)
    
    
    # Names of the variables
    colnames.var3 <- c("v", "HE", "PD", "HDOP", "Connectivity", "dist_err", "output")
    
    # Give the names of the fuzzy terms of each input variable.
    # Names of the fuzzy terms must be unique
    varinput3.1 <- c("high", "low", "zero")
    varinput3.2 <- c("small", "large")
    varinput3.3 <- c("short", "long")
    varinput3.4 <- c("good", "bad")
    varinput3.5 <- c("indirect", "direct")
    varinput3.6 <- c("small2", "large2")
    names.varinput3 <- c(varinput3.1, varinput3.2, varinput3.3, varinput3.4,
                         varinput3.5, varinput3.6)
    
    # Set interval of data, "v", "HE", "PD", "HDOP", "Connectivity", "dist_err", "output".
    range.data3 <- matrix(c(0, 50, 0, 360, 0, 60, 0, 20, 0, 1, 0, 1000, 0, 100), nrow=2)
    
    # Define number of fuzzy terms of input variables
    num.fvalinput3 <- matrix(c(3, 2, 2, 2, 2, 2), nrow=1)
    
    
    # Set the name of the simulation.
    name3 <- "Sim-3"
    
    
    
    # Define linear functions of TSK 
    # The dimension of this matrix is: 
    # [<number_of_rules>, <number_of_variables> + 1]
    func.tsk3 <- matrix(c(50, 10, 50, 10, 100, 
                          10, 10, 100, 100, 10), nrow = 10, byrow = TRUE)
    
    
    
    # Define the fuzzy IF-THEN rules; 
    # For TSK model, it isn't necessary to put linguistic terms in consequent parts.
    r1 <- c("high","and","small","and","dont_care","and","dont_care", 
            "and","dont_care","and","dont_care","->")
    r2 <- c("high","and","large", "and", "dont_care", "and", "dont_care", 
            "and","dont_care","and","dont_care","->")
    r3 <- c("dont_care", "and", "dont_care", "and", "short","and","good", 
            "and","dont_care","and","dont_care","->")
    r4 <- c("dont_care", "and", "dont_care", "and", "long","and","good", 
            "and","dont_care","and","dont_care","->")
    r5 <- c("dont_care", "and", "small","and","short", "and", "dont_care", 
            "and","dont_care","and","dont_care","->")
    r6 <- c("dont_care", "and", "large","and","long", "and", "dont_care", 
            "and","dont_care","and","dont_care","->")
    r7 <- c("dont_care","and","dont_care","and","dont_care","and","dont_care",
            "and","indirect","and","dont_care", "->")
    r8 <- c("dont_care","and","dont_care","and","dont_care","and","dont_care",
            "and","direct","and","dont_care", "->")
    r9 <- c("dont_care","and","dont_care","and","dont_care","and","dont_care","and",
            "dont_care","and","small2","->")
    r10 <- c("dont_care","and","dont_care","and","dont_care","and","dont_care","and",
             "dont_care","and","large2","->")
    rule3 <- list(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10)
    rule3 <- do.call(rbind, rule3)
    
    # Generate a fuzzy model with frbs.gen
    # For TSK model we do not need to input: 
    # num.fvaloutput, varout.mf, names.varoutput, type.defuz
    varinp.mf3 <- get("varinp.mf3", envir = cacheEnv)
    fis3 <- frbs.gen(range.data3, num.fvalinput3, names.varinput3, num.fvaloutput = NULL, 
                     varout.mf = NULL, names.varoutput = NULL, rule3, 
                     varinp.mf3, type.model, type.defuz = NULL, type.tnorm, type.snorm, func.tsk3, colnames.var3, 
                     type.implication.func, name3)
    
    assign("fis3", fis3, envir = cacheEnv)
  }
  
  create_fis3()
  
  
  
  ## Initial MapMatching Process: identifies the intial road link 
  imp <- function (traj, roads = "DigitalRoadNetwork", err_region) { 
    
    i <- 1
    count <- 0
    found <- FALSE
    initial_links <- data.frame(V1 = numeric(0), V2 = numeric(0), edge_id = numeric(0), 
                                direction = numeric(0), NP_x = numeric(0), NP_y = numeric(0))
    
    while (!found){#Enquanto o primeiro link n?o for encontrado execute:
      # Get coordinates of the current position and create SpatialPoints
      lon <- traj$coords.x1[i] #lon
      lat <- traj$coords.x2[i] #lat
      current_pt <- cbind(lon, lat)
      rec <- err_region(lon, lat)
      
      if (!requireNamespace("rgeos", quietly = TRUE))
        stop("package rgeos required")
      
      # Get edges inside the error region(contidos ou que interseptam o rec) 
      candidate_links <- data.frame(edge_id = unique(c(which(rgeos::gIntersects(rec, roads@sl, byid = TRUE)), 
                                                       which(rgeos::gContains(rec, roads@sl, byid = TRUE)))))
      
      # Nodes of the candidate links
      candidate_links$V1 <- get.edgelist(roads@g)[candidate_links$edge_id, 1]
      candidate_links$V2 <- get.edgelist(roads@g)[candidate_links$edge_id, 2]
      #V1 ? o ID do from.node.id e V2 to.node.id
      
      if (!requireNamespace("geosphere", quietly = TRUE))
        stop("package geosphere required")
      # Calculate the perpendicular distance from the current point to all 
      # segments inside the error region and the closest point on the segments
      PD <- sapply(candidate_links[,c("edge_id")], 
                   function(x) geosphere::dist2Line(current_pt, 
                                                    roads@sl@lines[[x]]@Lines[[1]]@coords))
      #PD ? a dist?ncia perpendicular do ponto para os segmentos que
      #est?o dentro da error region
      
      # Perpendicular distance
      candidate_links$PD <- PD[1,]
      # Nearest point
      candidate_links$NP_x <- PD[2,] #Longitude do Ponto mais pr?ximo no segmento 
      candidate_links$NP_y <- PD[3,] #Latitude do ponto mais pr?ximo no segmento
      
      # Calculate the beraing of the segments
      # If a segment is defined by the points a and b, bearing can be:
      # bearing(a,b) or bearing(b,a)
      # Which one is chosen depends on the differnce between the bearing and the GPS.Bearing
      gps_bearing <- traj$GPS.Bearing[i]
      candidate_links$direction <- sapply(candidate_links$edge_id, 
                                          function(x) {
                                            bearing <- geosphere::bearing(roads@sl@lines[[x]]@Lines[[1]]@coords[1,],
                                                                          roads@sl@lines[[x]]@Lines[[1]]@coords[2,])
                                            if (bearing - gps_bearing <= -90) {
                                              bearing <- bearing + 180
                                              if (bearing > 360) bearing <- bearing - 360
                                              bearing
                                            } else if (bearing - gps_bearing > 90) {
                                              bearing <- bearing - 180
                                              if (bearing < 0) bearing <- bearing + 360
                                              bearing
                                            } else bearing #Verificar como funciona o bearing
                                          }) 
      
      # Calculate the heading error
      candidate_links$HE <- abs(candidate_links$direction - traj$GPS.Bearing[i])
      
      speed <- traj$GPS.Speed[i]/3.6
      hdop <- traj$GPS.HDOP[i]
      
      # Prepare data for FIS1: speed, HE, PD, HDOP
      newdata <- cbind(rep(speed, nrow(candidate_links)), 
                       candidate_links$HE,
                       candidate_links$PD, 
                       rep(hdop, nrow(candidate_links)))
      
      
      
      # Get FIS1 from the cache environment
      fis1 <- get("fis1", envir = cacheEnv)
      
      # Probability of being the correct link
      # Sometimes warnings are produced saying the data is out of the specified range,
      # but these have no negativ impact on the link identification
      candidate_links$pred <- predict(fis1, newdata)$predicted.val
      initial_links[i,] <- candidate_links[candidate_links$pred == 
                                             max(candidate_links$pred),][,c("V1", "V2", "edge_id", "direction", "NP_x", "NP_y")]
      
      #plot(newdata[,2],newdata[,3])
      #points(newdata[as.numeric(rownames(initial_links)),2],newdata[as.numeric(rownames(initial_links)),3],col = "red")
      
      # Decide if  a link is the initial link depending on the number
      # of points matched to the link(quando o mesmo link ? assinalado duas
      #vezes consecutivas, ent?o o IMP aceita este link como solu??o inicial,
      #caso contr?rio, o contador ? zerado)
      
      if (i > 1) {
        if (identical(E(roads@g)[initial_links$edge_id[i]]$name, 
                      E(roads@g)[initial_links$edge_id[i - 1]]$name)) {
          count <- count +1
        } else {
          count <- 0
        } 
        if (count == 2) {
          initial_links <- initial_links[(i - 2):i,]
          found <- TRUE
        }
      }
      i <- i + 1
    }
    
    #Se encontrou a solu??o inicial, ent?o sai do while
    
    # Matched coordinates 
    traj$coords.x1[(i - 3):(i - 1)] <- initial_links$NP_x
    traj$coords.x2[(i - 3):(i - 1)] <- initial_links$NP_y
    #OSM ID of the roads 
    traj$OSM_ID[(i - 3):(i - 1)] <- E(roads@g)[initial_links$edge_id]$name
    
    list(traj = traj, index = i, current_link = initial_links[nrow(initial_links),])
  }
  
  
  
  ## Create the error region with fixed size using UTM coordinates
  ## 
  ## Args:
  ##  x: lon
  ##  y: lat
  ##
  ## Returns:
  ##  Error region as SpatialPolygon
  err_region <- function(x, y, size = 38) {  
    current_pt <- data.frame(x, y)
    coordinates(current_pt) <- ~x + y 
    proj4string(current_pt) <- osm_crs()
    
    # Transform to UTM
    #UTMzone <- trunc((180 + 5) / 6) + 1
    UTMzone <- long2UTM(x)
    UTM <- CRS(paste0("+proj=utm +zone=",UTMzone," +ellps=WGS84 +datum=WGS84"))
    current_pt <- spTransform(current_pt, UTM)
    x <- coordinates(current_pt )[1]
    y <- coordinates(current_pt )[2]
    
    # Create the error region with fixed size and transform back to osm_crs()
    #Cria as coordenadas do quadrado
    x_rec <- c(x - size, x + size, x + size, x - size, x - size)
    y_rec <- c(y - size, y - size, y + size, y + size, y - size)
    rec <- cbind(x_rec, y_rec)
    #Cria objeto da classe Polygons (sp)
    rec <- Polygons(list(Polygon(rec, hole = FALSE)), 1)
    rec <- SpatialPolygons(list(rec), proj4string = UTM)
    rec <- spTransform(rec, osm_crs())
  }
  
  
  long2UTM <- function(long) {
    (floor((long + 180)/6) %% 60) + 1
  }
  
  
  
  
  ## Subsequent MapMatching Process (SMP-1) along a link 
  smp1 <- function (traj, roads = "DigitalRoadNetwork", current_link, pt_index = "numeric") {
    
    last_fix <- cbind(traj$coords.x1[pt_index - 1], traj$coords.x2[pt_index - 1])
    current_pt <- cbind(traj$coords.x1[pt_index], traj$coords.x2[pt_index])
    edge_id <- current_link$edge_id
    name <- NULL # make R CMD check happy
    
    #Encontra o final do link atual(downstream junction)
    if (0 <= current_link$direction && current_link$direction <= 180) {
      
      current_link_end <- ifelse(V(roads@g)[name == current_link$V1]$lon 
                                 >= V(roads@g)[name == current_link$V2]$lon
                                 ,current_link$V1, current_link$V2)
    } else {
      current_link_end <- ifelse(V(roads@g)[name == current_link$V1]$lon 
                                 < V(roads@g)[name == current_link$V2]$lon
                                 ,current_link$V1, current_link$V2)
    }
    
    if (!requireNamespace("geosphere", quietly = TRUE))
      stop("package geosphere required")
    
    
    # location of the current position fix, relative to the link, as seen from
    # the last matched position
    
    alpha <- abs(geosphere::bearing(last_fix, current_pt) - current_link$direction)
    
    # location of the current position fix, relative to the link, as seen from
    # the downstream junction
    beta <- abs(geosphere::bearing(current_pt, cbind(V(roads@g)[name == current_link_end]$lon,
                                                     V(roads@g)[name == current_link_end]$lat))
                - current_link$direction)
    
    # Heading increment
    HI <- abs(traj$GPS.Bearing[pt_index] - traj$GPS.Bearing[pt_index - 1])
    
    # Distance from the last fix to the downstream junction (m)
    #Verificar esta parte
    d1 <- spDists(last_fix, cbind(V(roads@g)[name == current_link_end]$lon,
                                  V(roads@g)[name == current_link_end]$lat), 
                  longlat = TRUE) * 1000
    
    # Distance travelled since last position fix
    t <- as.double(traj$time[pt_index] - traj$time[pt_index-1])
    d2 <- (traj$GPS.Speed[pt_index]/3.6) * t
    
    # the difference between the distance from the last matched position to
    # the downstream junction and the distance travelled by the vehicle since
    # the last position fix
    delta_d <- d1 - d2
    
    speed <- traj$GPS.Speed[pt_index] / 3.6
    hdop <- traj$GPS.HDOP[pt_index]
    newdata <- cbind(speed, hdop, alpha, beta, delta_d, HI, HI)
    fis2 <- get("fis2", envir = cacheEnv)
    
    # Sometimes warnings are produced saying the data is out of the specified range,
    # but these have no negativ impact on the link identification
    pred_val <- predict(fis2, newdata)
    
    # Probability of matching the psotion fix on the current link
    pred_val
  }
  
  
  
  
  ## Subsequent MapMatching Process (SMP-2) at a junction 
  smp2 <- function(traj, roads = "DigitalRoadNetwork", current_link, pt_index = "numeric", err_region) {
    
    lon <- traj$coords.x1[pt_index] 
    lat <- traj$coords.x2[pt_index] 
    rec <- err_region(lon, lat, 38)
    
    current_pt <- cbind(lon, lat)
    last_fix <- cbind(traj$coords.x1[pt_index - 1], traj$coords.x2[pt_index - 1])
    edge_id <- current_link$edge_id
    
    # Current selected link becomes the previous link because SMP-2 chooses a new link
    prev_link <- current_link
    
    # make R CMD check happy
    name <- NULL
    #from <- NA 
    #rm(from)
    
    # Check which node of the prev_link is the end node
    if (0 <= prev_link$direction && prev_link$direction <= 180) {
      prev_link_end <- ifelse(V(roads@g)[name == prev_link$V1]$lon >= V(roads@g)[name == prev_link$V2]$lon
                              ,prev_link$V1, prev_link$V2)
    } else {
      prev_link_end <- ifelse(V(roads@g)[name == prev_link$V1]$lon < V(roads@g)[name == prev_link$V2]$lon
                              ,prev_link$V1, prev_link$V2)
    }
    
    if (!requireNamespace("rgeos", quietly = TRUE))
      stop("package rgeos required")
    
    # Get edges inside the error region 
    candidate_links <- data.frame(edge_id = unique(c(which(rgeos::gIntersects(rec, roads@sl, byid = TRUE)), 
                                                     which(rgeos::gContains(rec, roads@sl, byid = TRUE)))))
    
    
    
    # Nodes of the candidate links
    candidate_links$V1 <- get.edgelist(roads@g)[candidate_links$edge_id, 1]
    candidate_links$V2 <- get.edgelist(roads@g)[candidate_links$edge_id, 2]
    candidate_links <- candidate_links[!candidate_links$edge_id == edge_id,]
    
    
    if(length(candidate_links$edge_id) != 0){
    
    # Check if the line segments are connected to the prev_link
    candidate_links$conn <- sapply(candidate_links[,c("edge_id")],                         
                                   function(x) {
                                     # edges connected to the previously selected link
                                     conn_edges <- E(roads@g)[from(prev_link_end)] 
                                     if (isTRUE(any(as.vector(conn_edges) == x))) 1 else 0
                                   })
    
    if (!requireNamespace("geosphere", quietly = TRUE))
      stop("package geosphere required")
    
    # Calculate the perpendicular distance from the current point to all 
    # segments inside the error region and the closest point on the segments
    #str(candidate_links[,c("edge_id")])
    PD <- sapply(candidate_links[,c("edge_id")], 
                 function(x) geosphere::dist2Line(current_pt, roads@sl@lines[[x]]@Lines[[1]]@coords)) 
    
    #str(PD) -- EP: might be list(), which then breaks;
    #str(class(PD))
    #str(candidate_links) #mostra os links candidatos
    
    if (length(PD) == 0) {
      # Perpendicular distance
      #candidate_links$PD <- 1e9 # large
      # Nearest point
      #candidate_links$NP_x <- 0
      #candidate_links$NP_y <- 0
    } else {
      # Perpendicular distance
      candidate_links$PD <- PD[1,]
      # Nearest point
      candidate_links$NP_x <- PD[2,]
      candidate_links$NP_y <- PD[3,]
    }
    
    # Calculate the beraing of the segments
    # If a segment is defined by the points a and b, bearing can be:
    # bearing(a,b) or bearing(b,a)
    # Which one is chosen depends on the differnce between the bearing and the GPS.Bearing
    gps_bearing <- traj$GPS.Bearing[pt_index]
    candidate_links$direction <- sapply(candidate_links$edge_id, 
                                        function(x) {
                                          bearing <- geosphere::bearing(roads@sl@lines[[x]]@Lines[[1]]@coords[1,],
                                                                        roads@sl@lines[[x]]@Lines[[1]]@coords[2,])
                                          if (bearing - gps_bearing <= -90) {
                                            bearing <- bearing + 180
                                            if (bearing > 360) bearing <- bearing - 360
                                            bearing
                                          } else if (bearing - gps_bearing > 90) {
                                            bearing <- bearing - 180
                                            if (bearing < 0) bearing <- bearing + 360
                                            bearing
                                          } else bearing
                                        }) 
    
    # Calculate the heading error
    candidate_links$HE <- abs(candidate_links$direction - traj$GPS.Bearing[pt_index])
    
    # Distance (m) from last fix to the end node of the prev_link
    end_node <- cbind(V(roads@g)[name == prev_link_end]$lon, V(roads@g)[name == prev_link_end]$lat)
    d1 <- spDists(end_node, last_fix, longlat = TRUE) * 1000
    
    
    # Shortest path from prev_link_end to the segments and closest vertex of the segments
    sp <- as.data.frame(do.call(rbind, 
                                lapply(1:nrow(candidate_links), 
                                       function(x) {
                                         spV1 <- shortest.paths(roads@g, prev_link_end, candidate_links$V1[x])
                                         spV2 <- shortest.paths(roads@g, prev_link_end, candidate_links$V2[x])
                                         if (spV1 < spV2) {
                                           c(candidate_links$V1[x], spV1)
                                         } else {
                                           c(candidate_links$V2[x], spV2)}})))
    
    candidate_links$cl_vertex <- as.character(as.vector(sp[,1]))#Vertice mais proximo
    
    # length of the shortest path (m)
    candidate_links$sp <- as.numeric(as.vector(sp[,2]))
    candidate_links$sp[is.infinite(candidate_links$sp)] <- 500
    
    # Distance on the candidate links: from the start node of the link to
    # the nearest point on the link from the current position fix
    candidate_links$d_link <- apply(candidate_links[,c("NP_x", "NP_y", "cl_vertex")], 1,
                                    function(z) 
                                      spDists(cbind(V(roads@g)[name == z[3]]$lon,V(roads@g)[name == z[3]]$lat),
                                              cbind(as.numeric(z[1]), as.numeric(z[2])), 
                                              longlat = TRUE) * 1000)
    
    # Distance travelled since last position fix
    t <- as.double(traj$time[pt_index] - traj$time[pt_index-1])
    d <- (traj$GPS.Speed[pt_index]/3.6) * t
    
    # Distance error
    candidate_links$dist_err <- apply(candidate_links[,c("sp", "d_link")], 1,
                                      function(x) abs(d - (d1 + x[1] + x[2]))) 
    
    
    speed <- traj$GPS.Speed[pt_index] / 3.6
    hdop <- traj$GPS.HDOP[pt_index]
    
    # Prepare data for FIS3: speed, HE, PD, HDOP, connectivity, dist_err
    newdata <- cbind(rep(speed, nrow(candidate_links)), 
                     candidate_links$HE,
                     candidate_links$PD, 
                     rep(hdop, nrow(candidate_links)), 
                     candidate_links$conn, 
                     candidate_links$dist_err)
    
    fis3 <- get("fis3", envir = cacheEnv)
    # Probability of being the correct link
    # Sometimes warnings are produced saying the data is out of the specified range,
    # but these have no negativ impact on the link identification
    candidate_links$pred <- predict(fis3, newdata)$predicted.val
    
    # Current link the vehicle is traveling on
    current_link <- candidate_links[which.max(candidate_links$pred),c("V1", "V2", "edge_id", "direction", "NP_x", "NP_y")] 
    
    current_link 
    }else{
      
      print("no Candidate Links")
      current_link <- prev_link
      PD <- sapply(current_link[,c("edge_id")], 
                   function(x) geosphere::dist2Line(current_pt, roads@sl@lines[[x]]@Lines[[1]]@coords))
      current_link$NP_x <- PD[2,]
      current_link$NP_y <- PD[3,]
     
      current_link
      
    }
    
  }
  
  #=====================================================================
  #Executa o MapMatching
 
  DRN <- NULL
  plot <- FALSE
  
  if (!is(traj, "SpatialPointsDataFrame")) 
    stop ("Not a SpatialPointsDataFrame object!")
  if (is.null(proj4string(traj)))
    stop ("No projection specified!")
  if (!all(c("GPS.Bearing", "GPS.HDOP", "GPS.Speed") %in% names(traj))) 
    stop("Trajectory does not contain all the necessary data (GPS.Bearing, GPS.HDOP, GPS.Speed)!")
  if (!is(traj$time, "POSIXct") && !is(traj$time, "POSIXlt"))
    stop ("time must be of class POSIXct or POSIXlt!")
  
  traj <- spTransform(traj, osm_crs())
  coordnames(traj) <- c("coords.x1", "coords.x2")
  traj@data[is.na(traj@data)] <- 0
  bbox <- bbox(traj)#verificar se precisa mudar este bbox(usar cornerbox)
  
  # Create digital road network
  if (!is(DRN, "DigitalRoadNetwork")) {
    roads <- create_drn(bbox)#Aqui podemos delimitar area em que 
  } else {
    roads <- DRN
  }
  
  traj <- as(traj, "data.frame")
  traj$OSM_ID <- 0
  
  # Execute the Initial Map-Matching Process (IMP)
  list <- imp(traj, roads, err_region=38)
  edit_traj <- list$traj
  pt_index <- list$index
  current_link <- list$current_link
  
  # Map match remainig points using SMP1 and SMP2
  for (j in pt_index:nrow(edit_traj)) {
    # Check the possibility of matching the next point on the current link
    pred_val <- smp1(edit_traj, roads, current_link, j)$predicted.val
    
    #print(j)
    #Caso pred_val(SMP1) for maior ou igual 60, ent?o assinala o ponto ao
    #link encontrado por ela
    if (pred_val >= 60) {
      if (!requireNamespace("geosphere", quietly = TRUE))
        stop("package geosphere required")
      PD <- geosphere::dist2Line(edit_traj[,c("coords.x1", "coords.x2")][j,], 
                                 roads@sl@lines[[current_link$edge_id]]@Lines[[1]]@coords)
      edit_traj$coords.x1[j] <- PD[2]
      edit_traj$coords.x2[j] <- PD[3]
      edit_traj$OSM_ID[j] <- edit_traj$OSM_ID[j - 1]
      #pt_index <- pt_index + 1
    } else {
      current_link <- smp2(edit_traj, roads, current_link, j, 100)
      edit_traj$coords.x1[j] <- current_link$NP_x
      edit_traj$coords.x2[j] <- current_link$NP_y
      edit_traj$OSM_ID[j] <- E(roads@g)[current_link$edge_id]$name
      #pt_index <- pt_index + 1
    }
  }
  
  data <- edit_traj[,!names(edit_traj) %in% c("coords.x1", "coords.x2")]
  matched_coords <- edit_traj[,c("coords.x1", "coords.x2")]
  matched_traj <- SpatialPointsDataFrame(matched_coords, data, proj4string=osm_crs())
  
  if (plot) {
    plot(traj$coords.x1, traj$coords.x2, pch = 16, col = "blue")
    points(matched_traj$coords.x1, matched_traj$coords.x2,pch = 16, col = "red")
    lines(roads@sl)
  }
  #return(matched_traj)
  saida <- as.data.frame(matched_traj)
  lat_gps <- traj$coords.x2
  lon_gps <- traj$coords.x1
  
  lat_cal <- saida$coords.x2
  lon_cal <- saida$coords.x1
  
  id_way <- saida$OSM_ID
  
  saida2 <- cbind(lat_gps,lon_gps,lat_cal,lon_cal,id_way)
  colnames(saida2) <- c('lat_gps','lon_gps','lat_cal','lon_cal','id_way')
  return(saida2)
}