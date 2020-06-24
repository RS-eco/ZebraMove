#' R-Code - MSc Thesis 
#' RS-eco
#' January 2017
## ---- markdown_options ----

# Set markdown setttings
knitr::opts_chunk$set(cache=TRUE, eval=FALSE, warning=FALSE, 
                      message=FALSE, comment=NA, 
                      tidy=FALSE, results="hide")

## ---- global_options ----

# Define working directory and file directory
workdir <- "/home/mabi/GitHub/ZebraMove/"
filedir <- "/home/mabi/Documents/Wissenschaft/Data/"

# Path specification of working directory
setwd(workdir)

# Load required libraries
sapply(c("doParallel", "ellipse", "gdalUtils", "ggmap", "ggpmisc", 
         "ggplot2", "ggrepel", "grid", "lubridate",  "move", "rgdal", 
         "raster", "rgeos", "RStoolbox", "xtable"), 
       require, character.only = TRUE)

# Define spatial projections
crs.wgs84 <- CRS("+proj=longlat +datum=WGS84 +no_defs 
                 +ellps=WGS84 +towgs84=0,0,0")
crs.laea <- CRS("+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 
                +datum=WGS84 +units=m +no_defs 
                +ellps=WGS84 +towgs84=0,0,0")

# Set number of cores to use
n <- parallel::detectCores()-2

# Define maximum pixels for plotting
maxpixels <- 1000000

# Define colour theme for plotting
colourtheme <- viridis::viridis(255)

# Set plotting theme
theme_set(theme_bw() + 
            theme(panel.border = element_rect(colour = "transparent")))

# Iso codes of countries
isos <- c("BWA", "KEN")
region <- c("Ngamiland", "Laikipia")

# Create north arrow plot
library(png)
northarrow <- readPNG("Pictures/north-arrow.png")
g <- rasterGrob(northarrow, interpolate=TRUE); rm(northarrow)
h <- ggplot(data=data.frame(x = 1:5, y = 1:10), aes(x=x, y=y)) + 
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + 
  theme(panel.background = element_rect(fill = "transparent",colour = NA), 
        panel.grid.minor = element_line(colour = NA), 
        panel.grid.major = element_line(colour = NA),
        axis.ticks = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        panel.border = element_rect(fill=NA, colour = NA), 
        plot.background = element_rect(fill = "transparent",colour = NA), 
        legend.spacing=unit(-0.5, "lines")) + labs(x=NULL, y=NULL); rm(g)

## ---- functions ----
#Build Function to Return Element Text Object
rotatedAxisElementText = function(angle,position='x'){
  angle = angle[1];
  position = position[1]
  positions = list(x=0,y=90,top=180,right=270)
  if(!position %in% names(positions))
    stop(sprintf("'position' must be one of [%s]",
                 paste(names(positions),collapse=",")),call.=FALSE)
  if(!is.numeric(angle))
    stop("'angle' must be numeric",call.=FALSE)
  rads  = (angle - positions[[ position ]])*pi/180
  hjust = 0.5*(1 - sin(rads))
  vjust = 0.5*(1 + cos(rads))
  element_text(angle=angle,vjust=vjust,hjust=hjust)
}

# Get Gimms Data for certain time period and area
gimms3g <- function(extent = c(10, 30, -15, 15), 
                    path="/home/mabi/Documents/Wissenschaft/Data/GIMMS",
                    snap="near"){
  #Create empty raster stacks
  library(raster)
  ndvi_sa_15day <- stack()
  
  # specify temporal parameters
  startyear = 1982
  endyear = 2013
  year = startyear:endyear
  year.name = c("82","83","84","85","86","87","88","89","90",
                "91","92","93","94","95", "96","97","98","99",
                "00","01","02","03","04", "05","06","07","08",
                "09","10","11","12","13")
  month = c("01", "02", "03", "04", "05", "06" ,"07","08","09",
            "10","11","12")
  month.name = c("jan", "feb", "mar", "apr", "may", "jun",
                 "jul","aug", "sep","oct","nov","dec")
  
  # Define extent
  extent <- extent(extent)
  
  # set to loop through years and months within each year
  for (i in 1:length(year.name)){
    print(year[i])
    for (j in 1:length(month)){
      print(month[j])
      
      ### process "a"-lettered files
      letter = 'a'
      
      # set name format 
      file.path = list.files(path, 
                             pattern=paste0("geo", year.name[i], 
                                            month.name[j], 15, letter), 
                             full.names=TRUE)
      
      # Read binary file
      binread = readBin(con=file.path, what="integer", n=9331200, 
                        size=2, signed=T, endian="big")
      
      # covert binary object into a matrix using specifications 
      #from the metadata
      mtrx = matrix(binread, nrow=2160, ncol=4320, byrow=F)
      
      # convert from matrix to a raster object
      rstr = raster(mtrx)
      
      # specify extent and projection
      extent(rstr) <- extent(c(-180, 180, -90, 90))
      projection(rstr) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 
      +no_defs"
      
      # set NA values
      rstr[rstr == -5000] <-- NA
      
      #scaling to actual value
      rstr = floor(rstr/10)/1000
      rstr.name = paste0("NDVI", "_", year[i], "_", month[j], "_", letter)
      
      #Cut raster to extent of studyarea
      rstr_sa <- crop(rstr, extent, snap=snap); rm(rstr)
      
      #Save 15day NDVI in raster stack
      ndvi_sa_15day <- stack(ndvi_sa_15day, rstr_sa); rm(rstr_sa)
      
      # free memory
      rm(letter,file.path, mtrx, rstr.name)
      
      ### process "b"-lettered files
      letter = 'b'
      
      # set name format
      file.path = list.files(path, pattern=
                               paste0("geo", year.name[i], 
                                      month.name[j], 15, letter), 
                             full.names=TRUE)
      
      # Read binary file
      binread = readBin(con=file.path, what="integer", n=9331200, size=2, 
                        signed=T, endian="big")
      
      # covert binary object into a matrix using specifications from 
      #the metadata
      mtrx = matrix(binread, nrow=2160, ncol=4320, byrow=F)
      
      # convert to matrix to a raster object
      rstr = raster(mtrx)
      
      # specify extent and projection
      extent(rstr) <- extent(c(-180, 180, -90, 90))
      projection(rstr) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 
      +no_defs"
      
      # set NA values
      rstr[rstr == -5000] <-- NA
      
      # division by 10000 for NDVI sets the actual value
      rstr = floor(rstr/10)/1000
      rstr.name = paste0("NDVI", "_", year[i], "_", month[j], "_", letter)
      
      #Cut raster to extent of study area
      rstr_sa <- crop(rstr, extent, snap=snap)
      rm(rstr)
      
      #Save 15day NDVI in raster stack
      ndvi_sa_15day <- stack(ndvi_sa_15day, rstr_sa)
      rm(rstr_sa)
      
      # free memory
      rm(letter,file.path,mtrx,rstr.name)
    }
  }
  # Clear memory
  rm(binread, i, j, year, month, month.name, year.name)
  
  #Set time of raster values 
  ndvi_sa_15day <- setZ(ndvi_sa_15day, 
                        seq(from = as.Date(paste(startyear, "01-15", sep="-")), 
                            to = as.Date(paste(endyear, "12-30", sep="-")), 
                            length=24*(endyear-startyear+1)))
  names(ndvi_sa_15day) <- format(
    seq(from = as.Date(paste(startyear, "01-15", sep="-")), 
        to = as.Date(paste(endyear, "12-30", sep="-")), 
        length=24*(endyear-startyear+1)), '%Y%m%d')
  return(ndvi_sa_15day)
}

# Add Scale bar to ggplot2
ggScaleBar = function(ggplot_obj, extent_obj, addParams =
                        list()) {
  addParamsDefaults = list(noBins = 5, unit = "m",
                           placement = "bottomright", sbLengthPct = 0.3, 
                           sbHeightvsWidth = 1/14, size=6, 
                           margin=1, offset_unit=1)
  addParams = modifyList(addParamsDefaults, addParams)
  
  range_x = extent_obj[2]-extent_obj[1]
  range_y = extent_obj[4]-extent_obj[3]
  lengthScalebar = addParams[["sbLengthPct"]] * range_x
  ## OPTION: use pretty() instead
  
  makeNiceNumber = function(num, num.pretty = 1) {
    # Rounding provided by code from Maarten Plieger
    return((round(num/10^(round(log10(num))-1))*(10^(round(log10(num))-1))))
  }
  widthBin = makeNiceNumber(lengthScalebar / addParams[["noBins"]])
  heightBin = lengthScalebar * addParams[["sbHeightvsWidth"]]
  if(addParams$placement == "bottomleft"){
    ScaleBar = c(x = extent_obj[1]+addParams[["margin"]], 
                 y = extent_obj[3]+addParams[["margin"]])
  }else if(addParams$placement == "bottomright"){
    ScaleBar = c(x = extent_obj[2], y = extent_obj[3])
  }
  createBoxPolygon = function(llcorner, width, height) {
    relativeCoords = data.frame(c(0, 0, width, width, 0), 
                                c(0, height, height, 0, 0))
    names(relativeCoords) = names(llcorner)
    return(t(apply(relativeCoords, 1, function(x) llcorner + x)))
  }
  scaleBarPolygon = do.call("rbind", lapply(0:(addParams[["noBins"]] - 1), 
                                            function(n) {
    dum = data.frame(createBoxPolygon(ScaleBar + c((n * widthBin), 0), 
                                      widthBin, heightBin))
    if(!(n + 1) %% 2 == 0) dum$cat = "odd" else dum$cat = "even"
    return(dum)
  }))
  if(addParams$unit == "m"){
    textScaleBar = data.frame(
      x = ScaleBar[["x"]] + 
        (c(c(0:(addParams[["noBins"]])) * widthBin, 
           addParams[["noBins"]]* widthBin + 
             widthBin/addParams[["offset_unit"]])), 
      y = ScaleBar[["y"]], 
      label = c(as.character(0:(addParams[["noBins"]]) * widthBin), 
                paste0("m"))) 
  } else if(addParams$unit == "km"){
    textScaleBar = data.frame(
      x = ScaleBar[["x"]] + (c(c(0:(addParams[["noBins"]])) * widthBin, 
                               addParams[["noBins"]]* widthBin + 
                                 widthBin/addParams[["offset_unit"]])), 
      y = ScaleBar[["y"]], 
      label = c(as.character(0:(addParams[["noBins"]]) * widthBin / 1000), 
                paste0("km"))) 
  }
  return(ggplot_obj +
           geom_polygon(data = subset(scaleBarPolygon, cat == "odd"), 
                        aes(x=x, y=y), fill = "black", color = "black") +
           geom_polygon(data = subset(scaleBarPolygon, cat == "even"), 
                        aes(x=x, y=y), fill = "white", color = "black") +
           geom_text(aes(x=x, y=y, label = label), color = "black", 
                     size = addParams[["size"]], data = textScaleBar, 
                     hjust = 0.5, vjust = 1.2))
}

# Obtain cleaned GBIF data
terrestrialGBIF <- 
  function(species="Equus quagga", genus=NULL, limit = 50000, 
           path = "/home/mabi/Documents/Wissenschaft/Data/", 
           overwrite = FALSE){
  # Get data of species
  species_name <- paste0(strsplit(species, split=" ")[[1]][1], ".", 
                         strsplit(species, split=" ")[[1]][2])
  if(overwrite == FALSE & file.exists(paste0(path, "GBIF/", 
                                             species_name, ".rds"))){
    sp_data <- readRDS(paste0(path, "GBIF/", species_name, ".rds"))
  } else{
    # Download species location data from gbif
    library(rgbif)
    if(is.null(species)){
      species_data <- 
        occ_search(taxonKey=name_backbone(name=species)$speciesKey, 
                   return="data",
                   fields=c("species", "year", "month", 
                            "decimalLatitude", "decimalLongitude"), 
                   limit= limit)
    } else if(is.null(genus)){
      species_data <- 
        occ_search(genusKey=name_backbone(name=species)$genusKey, 
                   return="data", 
                   fields=c("species", "year", "month", 
                            "decimalLatitude", "decimalLongitude"), 
                   limit= limit)
    } else{print("No species or genus specified!");break()}
    
    # Discard data with errors in coordinates
    species_data <- 
      species_data[complete.cases(
        species_data[,c("decimalLongitude", "decimalLatitude")]),]
    
    # Set spatial coordinates
    coordinates(species_data) <- c("decimalLongitude", "decimalLatitude")
    
    # Define spatial projection
    proj4string(species_data) <- 
      CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 
          +no_defs +towgs84=0,0,0")
    
    # Save species records as RDS-File
    saveRDS(species_data, file = paste0(path, "GBIF/", species_name, ".rds"))
    rm(species_data)
    sp_data <- readRDS(paste0(path, "GBIF/", species_name, ".rds"))
  }
  
  # Get high resolution country data
  data(countriesHigh, package="rworldxtra")
  
  # Remove species data at sea
  country_per_point <- over(sp_data, countriesHigh)
  sp_data <- sp_data[!is.na(country_per_point$SOV_A3),]
  rm(country_per_point)
  
  # Remove outliers
  library(outliers)
  sp_data <- sp_data[sp_data$decimalLongitude != 
                       outlier(sp_data$decimalLongitude),]
  sp_data <- sp_data[sp_data$decimalLatitude != 
                       outlier(sp_data$decimalLatitude),]
  
  # Return data
  return(sp_data)
}

# SSF Functions

thin_split <- function(object, target_interval, 
                       upper_threshold = 1.1 * target_interval, 
                       units = "min",
                       minimum_burst_length = 3) {
  # target interval can be specified as an integer, it is 
  # pasted with units to provide a value for the seq-algorithm: 
  # e.g. as "60 mins"
  # TODO: currently we throw away more observations than necessary
  time <- timestamps(object)
  interval_id <- cut(time, seq(min(time), max(time), 
                               by = paste(target_interval, units)), 
                     include.lowest = TRUE, labels = FALSE)
  object <- object[!duplicated(interval_id)]
  dt <- timeLag(object, units = units)
  burstId <- c(0, cumsum(dt > upper_threshold)) 
  # burstId is incremented whenever the time lag to the next observation 
  #is larger than upper_threshold
  # Remove bursts that are shorter than the minimum_burst_length
  remove <- 
    which(burstId %in% names(which(table(burstId) < minimum_burst_length)))
  if (length(remove) > 0) {
    object <- object[-remove, ]  
    burstId <- burstId[-remove]
  }
  object$burstId <- factor(burstId)
  object
}

rdist_2d <- function(n, dist = "exponential", theta, kappa = NULL, 
                     start_angle = NULL) {
  # n: number of samples to draw
  # lambda: rate parameter, mean distance is 1/lambda
  r <- switch(dist,
              exponential = rexp(n, theta),
              gamma = rgamma(n, theta[1], theta[2]),
              lognormal = rlnorm(n, theta[1], theta[2]),
              halfnormal = abs(rnorm(n, 0, sqrt(pi/2)/theta)))
  if (is.null(kappa)) {
    theta <- runif(n, 0, 2*pi)  
  } else {
    if(is.null(start_angle)) stop("A start_angle is required when 
                                  kappa is not NULL")
    theta <- rvonmises(n, circular(0), kappa) + start_angle
  }
  dx <- r * cos(theta)
  dy <- r * sin(theta)
  cbind(dx, dy)  
}

wrap_rel_angle <- function(x) {
  # wraps relative angles to the range from -pi to pi
  ifelse(abs(x) > pi, sign(x) * (abs(x) - 2 * pi), x)
}

relative_control_locations <- 
  function(K, method = c("exponential", "gamma", "lognormal", "halfnormal"), 
                                       start_angle = NULL, theta = NULL, 
           kappa = NULL) {
  dxy <- rdist_2d(K, dist = method, theta, kappa, start_angle)
  dxy
}

degrees2radians <- function(x) {
  -2*pi*(x-90)/360
}

radians2degrees <- function(x) {
  round((x*360/-2/pi)+90)
}

move2spdf <- function(object) {
  # Transform a move object to a SpatialPointsDataFrame
  # date, id, burst, abs.angle
  if (!is.null(object$individual.local.identifier)) {
    id <- factor(as.character(object$individual.local.identifier))
  } else {
    id <- "unknown"
  }
  res <- data.frame(coordinates(object), 
                    "date" = timestamps(object),
                    "id" = id,
                    "burst" = object$burstId,
                    "dist" = c(distance(object), NA),
                    "abs.angle" = c(degrees2radians(angle(object)), NA)                    
  )
  colnames(res)[1:2] <- c("x", "y")
  coordinates(res) <- ~ x + y
  proj4string(res) <- proj4string(object)
  res  
}

fast_sample_path <- function(path, angle_dist = NULL, K, 
                             method = c("exponential", "gamma", 
                                        "halfnormal", "lognormal", 
                                        "empirical"), theta = NULL, 
                             kappa = NULL) {
  # For each step in a path, create K alternative steps 
  # according to method ("exponential" or "empirical")
  xy_path <- coordinates(path)  # xy coordinates 
  # abs_angles <- atan2(diff(xy_path[, "y"]), diff(xy_path[, "x"]))
  abs_angles <- path$abs.angle
  start_angles <- rep(abs_angles[-c(length(abs_angles)-1, 
                                    length(abs_angles))], each = K)
  path$used <- 1
  path$rel.angle <- wrap_rel_angle(c(NA, diff(abs_angles)))
  path$stratum <- 1:NROW(path)
  
  # The values in the path object for distance and angle are prospective.
  # In the analysis we need retrospective values 
  #(how far did we need to move to get to the new position?)
  dxy <- 
    relative_control_locations(K * (NROW(xy_path) - 2), method, 
                               start_angle = start_angles, theta, kappa)
  start_i <- rep(2:(NROW(xy_path)-1), each = K) 
  # So we ignore the first and last position
  stop_i <- rep(3:(NROW(xy_path)), each = K)
  # the first because we do not know its relative angle, 
  #the last because we do not know where the animal moved to further
  xy <- data.frame("x" = dxy[, 1] + xy_path[start_i, 1], 
                   "y" = dxy[, 2] + xy_path[start_i, 2])
  xy$date <- path$date[start_i]
  xy$used <- 0
  xy$dist <- sqrt((xy$x - xy_path[start_i, 1])^2 + 
                    (xy$y - xy_path[start_i, 2])^2)  
  xy$id <- path$id[start_i]
  xy$burst <- path$burst[start_i]
  abs.angle <- atan2(xy$y - xy_path[start_i, 2], 
                     xy$x - xy_path[start_i, 1])
  xy$rel.angle <- wrap_rel_angle(abs.angle - path$abs.angle[start_i - 1])
  xy$stratum <- path$stratum[start_i]
  realized_path <- data.frame(xy_path[-c(1:2), , drop = FALSE])
  realized_path <- 
    cbind(realized_path, as.data.frame(path)[2:(NROW(path)-1), 
                                             c("date", "used", "id", 
                                               "stratum", "burst", "dist", 
                                               "rel.angle"), drop = FALSE])  
  cols <- c("x", "y", "date", "used", "id", "stratum", "burst", "dist", 
            "rel.angle")
  res <- rbind(realized_path[, cols], xy[, cols])
  coordinates(res) <- ~ x + y
  proj4string(res) <- proj4string(path)
  res
}

extract_angle_dist <- function(path) {
  dist <- distance(path)
  abs.angle <- angle(path)
  abs.angle <- degrees2radians(abs.angle)
  rel.angle <- wrap_rel_angle(c(NA, diff(abs.angle)))
  burstIds <- path$burstId
  # we remove "steps" across bursts
  remove <- c(1, which(diff(as.numeric(burstIds)) != 0) + 1)
  remove <- c(remove, remove + 1)
  angle_dist <- cbind("dist" = dist[-remove], 
                      "rel.angle" = rel.angle[-remove])
  angle_dist
}

prepare_angle_dist <- function(path_object) {
  if (class(path_object) == "MoveStack") {
    path_list <- split(path_object)
  }
  if (class(path_object) == "Move") {
    path_list <- list(path_object)
  }
  angle_dist <- lapply(path_list, extract_angle_dist)
  angle_dist <- do.call("rbind", angle_dist)
  na.omit(angle_dist)    
}

prepare_ssf_steps <- function(path_object, method, K, angle_dist = NULL, 
                              theta = NULL, kappa = NULL, disc_radius = NULL, 
                              crs, verbose = FALSE) {  
  # TODO: Rethink whether to do everything in lat-lon to keep projection as 
  # in move objects
  orig_proj4string <- proj4string(path_object)
  path_object <- spTransform(path_object, crs)    
  if (class(path_object) == "MoveStack") { # convert move stact to list 
    path_list <- split(path_object)
  }
  if (class(path_object) == "Move") {
    path_list <- list(path_object)
  }
  
  ssf_path_list <- list(length = length(path_list))
  for(animal in seq(path_list)) {      # for each individual 
    animal_path <- path_list[[animal]]
    if(is.null(animal_path$burstId)) {
      animal_path$burstId <- factor(0)
    }
    burstIds <- animal_path$burstId
    burst_path_list <- lapply(unique(as.character(burstIds)), 
                              function(burst) {
      if (verbose) cat("Processing animal", animal, "burst", burst, "\n")
      path <- animal_path[animal_path$burstId == burst, ]
      path <- move2spdf(path)
      path <- spTransform(path, crs)    
      ssf_path <- fast_sample_path(path, angle_dist, K = K, method = method, 
                                   theta = theta, kappa = kappa)
      ssf_path      
    })
    ssf_path_list[[animal]] <- do.call("rbind", burst_path_list)
  }
  ssf_data <- do.call("rbind", ssf_path_list)  
  ssf_data$stratum <- paste(ssf_data$id, ssf_data$burst, ssf_data$stratum, 
                            sep="_")
  ssf_data <- spTransform(ssf_data, CRS(orig_proj4string))
  ssf_data
}

## ---- clean_move ----
# Specify data files
data <- c("Data/MigratoryBurchellszebra_northernBotswana.csv", 
          "Data/Zebras of Laikipia-Samburu, Kenya.csv")

# Loop for reading each dataset into R
for(i in 1:length(data)){
  if(i == 1){
    # Load Zebra data manually downloaded from Movebank
    zebra_data <-read.csv(data[i], header = TRUE, as.is = TRUE)
    # Clean zebra data
    zebra_data <- subset(zebra_data, 
                         select=-c(utm.northing, utm.easting, utm.zone))
    zebra_data <- subset(zebra_data, !is.na(zebra_data$location.long) & 
                           !is.na(zebra_data$location.lat) & 
                           manually.marked.outlier != "true")
    
    # Convert zebra data to move object
    zebra_data <- move(
      # set columns for coordinates
      x = zebra_data$location.long, y = zebra_data$location.lat,
      # conversion of Date and Time from character to POSIX object
      time = as.POSIXct(zebra_data$timestamp, 
                        format="%Y-%m-%d %H:%M:%S", tz="UTC"),
      # set projection
      proj = "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0",
      # set id
      animal = zebra_data$individual.local.identifier,
      # define data to use and attach to the spatial object
      data = zebra_data)
  }
  else if(i == 2){
    # Load zebra data
    zebra_data <- read.csv(data[i])
    zebra_data <- zebra_data[-5212,] # remove one identical duplicate
    
    # Convert zebra data to move object
    zebra_data <- move(
      # set columns for coordinates
      x = zebra_data$location.long, y = zebra_data$location.lat,
      # conversion of Date and Time from character to POSIX object
      time = as.POSIXct(zebra_data$timestamp, format="%Y-%m-%d %H:%M:%S", 
                        tz="EAT"),
      # set projection
      proj = "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0",
      # set id
      animal = zebra_data$individual.local.identifier,
      # define data to use and attach to the spatial object
      data = zebra_data)
    
    # Convert zebra data to spatial object
    zebra3 <- read.csv("Data/OPC_collar_data.csv")
    colnames(zebra3) <- c("individual.local.identifier", "year", 
                          "month", "day", 
                          "hour", "minute", "utm.easting", "utm.northing")
    levels(zebra3$individual.local.identifier) <- c("PZ10", "PZ14", "PZ6", 
                                                    "PZ8")
    zebra3$utm.northing <- zebra3$utm.northing - 10000000 
    # Correct Northing to Arc projection
    
    # Add latitude and longitude to data
    longlat <- zebra3
    coordinates(longlat) <- ~utm.easting+utm.northing
    projection(longlat) <- CRS("+init=epsg:21097")
    longlat <- spTransform(longlat, crs.wgs84)
    zebra3$location.long <- data.frame(coordinates(longlat))[,1]
    zebra3$location.lat <- data.frame(coordinates(longlat))[,2]; rm(longlat)
    
    # Create date and time vector and combine it to timestamp
    zebra3$date <- 
      as.Date(with(zebra3, paste(year, month, day, sep="-")), 
              "%Y-%m-%d")
    zebra3$time <- with(zebra3, paste(hour, minute, "00", sep=":"))
    zebra3$timestamp <- 
      as.POSIXct(paste(zebra3$date, zebra3$time), 
                 format="%Y-%m-%d %H:%M:%S", tz="EAT")
    
    # Create move object of our data
    zebra3 <- move(
      # set columns for coordinates
      x = zebra3$location.long, y = zebra3$location.lat,
      # conversion of Date and Time from character to POSIX object
      time = as.POSIXct(zebra3$timestamp, 
                        format="%Y-%m-%d %H:%M:%S", tz="EAT"),
      # set projection
      proj = "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0",
      # set id
      animal = zebra3$individual.local.identifier,
      # define data to use and attach to the spatial object
      data = zebra3)
    
    # Merge zebra data together
    zebra_data <- moveStack(c(split(zebra_data), split(zebra3))); rm(zebra3)
  }
  # Create spatial extent of our study area
  studyarea <- 
    spatial.tools::bbox_to_SpatialPolygons(
      extent(zebra_data)*1.05, crs.wgs84);
  studyarea <- as(studyarea, "SpatialPolygonsDataFrame")
  writeOGR(studyarea, dsn = "Results/", 
           layer = paste0("studyarea_", isos[i]), 
           driver="ESRI Shapefile", overwrite=TRUE)
  studyarea_laea <- spTransform(studyarea, crs.laea)
  writeOGR(studyarea_laea, dsn = "Results/", 
           layer = paste0("studyarea_laea_", isos[i]), 
           driver="ESRI Shapefile", overwrite=TRUE); rm(studyarea)
  
  # Spatial extent of study area
  gArea(studyarea_laea)
  
  # Save movement data to file
  saveRDS(zebra_data, paste0(workdir, "Data/zebra_data_", isos[i], ".rds"))
  
  # Transform move object into laea projection
  zebra_data_laea <- spTransform(zebra_data, crs.laea); rm(zebra_data)
  saveRDS(zebra_data_laea, 
          paste0(workdir, "Data/zebra_data_laea_", isos[i], ".rds"))
  rm(zebra_data_laea)
}

## ---- gbif ----

# Get GBIF Zebra data
species_data <- terrestrialGBIF(species=c("Equus quagga", "Equus burchelli", 
                                          "Equus burchellii", "Equus grevyi"))

species_data <- species_data[species_data$decimalLatitude < 10,]
species_data <- species_data[species_data$decimalLongitude > 0,]
species_data <- spTransform(species_data, crs.laea)

zebra_1 <- readRDS(paste0(workdir, "Data/zebra_data_", isos[1], ".rds"))
zebra_1 <- spTransform(zebra_1, crs.laea)
zebra_2 <- readRDS(paste0(workdir, "Data/zebra_data_", isos[2], ".rds"))
zebra_2 <- spTransform(zebra_2, crs.laea)

# Load countries data of Africa
data(countriesHigh, package="rworldxtra")
Africa <- subset(countriesHigh, continent=="Africa"); rm(countriesHigh)
Africa_laea <- spTransform(Africa, crs.laea)

# Add study area shapefiles to plot
studyarea_1 <- readOGR(dsn=paste0(workdir, "Results/"), 
                       layer="studyarea_laea_BWA")
studyarea_2 <- readOGR(dsn=paste0(workdir, "Results/"), 
                       layer="studyarea_laea_KEN")

# Add protected areas to plot
if(!dir.exists(paste0(workdir, "Results/wdpa_laea_Africa"))){
  wdpa <- readOGR(dsn=paste0(filedir, "WDPA/WDPA_Oct2016_Africa-shapefile"), 
                  layer="WDPA_Oct2016_Africa-shapefile-polygons", 
                  verbose=FALSE)
  wdpa <- wdpa[wdpa$MARINE == 0,]
  # Add only IUCN categorised ones!
  wdpa <- crop(wdpa, Africa); rm(Africa)
  wdpa_laea <- spTransform(wdpa, crs.laea)
  writeOGR(wdpa_laea, dsn = paste0(workdir, "Results/wdpa_laea_Africa"), 
           layer="wdpa", driver="ESRI Shapefile", check_exists=TRUE, 
           overwrite_layer=TRUE)
}
wdpa_laea <- readOGR(dsn=paste0(workdir, "Results/wdpa_laea_Africa"), 
                     layer="wdpa_laea_Africa", verbose=FALSE)

# Create zebra plot
zebra <- readPNG("Pictures/Zebra-PNG-Image.png")
j <- rasterGrob(zebra, interpolate=TRUE); rm(zebra)
k <- ggplot(data=data.frame(x = 1:5, y = 1:10), aes(x=x, y=y)) + 
  annotation_custom(j, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + 
  theme(panel.background = element_rect(fill = "transparent",colour = NA), 
        panel.grid.minor = element_line(colour = NA), 
        panel.grid.major = element_line(colour = NA), 
        axis.ticks = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.x = element_blank(), axis.text.y = element_blank(), 
        panel.border = element_rect(fill=NA, colour = NA), 
        plot.background = element_rect(fill = "transparent",colour = NA), 
        legend.spacing=unit(-0.5, "lines")) + labs(x=NULL, y=NULL); rm(j)

# Plot species data
a <- ggplot() + geom_polygon(data=Africa_laea, 
                             aes(x=long, y=lat, group=group), 
                             fill="gray", colour="black") + 
  geom_polygon(data=wdpa_laea, aes(x=long, y=lat, group=group), 
               fill="darkseagreen2") + 
  geom_polygon(data=studyarea_1, aes(x=long, y=lat, group=group),
               fill="transparent", colour="red") + 
  geom_polygon(data=studyarea_2, aes(x=long, y=lat, group=group), 
               fill="transparent", colour="red") + 
  geom_point(data = data.frame(species_data), aes(x = decimalLongitude, 
                                                  y = decimalLatitude, 
                                                  colour = "gbif"), 
             size=0.7) + geom_path(aes(x = coords.x1, y = coords.x2, 
                                       colour="move"), size=0.7, 
                                            data = data.frame(zebra_1)) + 
  geom_path(aes(x = coords.x1, y = coords.x2, colour="move"), size=0.7, 
            data = data.frame(zebra_2)) + 
  scale_x_continuous(name="Eastings (m)", limits=c(-4750000, 3750000), 
                     expand=c(0,0)) + 
  scale_y_continuous(name="Northings (m)", limits=c(-4750000, 3750000), 
                     expand=c(0,0)) + 
  scale_colour_manual(name="", values=c(gbif="#F8766D", move="#00BFC4"), 
                      labels=c(gbif="GBIF", move="Movebank")) + 
  guides(color=guide_legend(override.aes=list(shape=c(16,NA), 
                                              linetype=c(0,1)))) + 
  theme(legend.justification=c(1,1), legend.position=c(1,1), 
        legend.background = element_rect(fill="transparent", colour=NA), 
        axis.text.y = rotatedAxisElementText(90,"y")) + 
  coord_equal(expand=FALSE)

b <- ggScaleBar(a, extent(Africa_laea)/ 1.5, 
                addParams = list(noBins = 3, placement="bottomleft", 
                                 unit="km", sbLengthPct = 0.3, 
                                 sbHeightvsWidth = 1/30, size=4))

cairo_pdf(paste0("Figures/Dist_Zebra.pdf"), width=6, height=6, 
          bg="transparent")
grid.newpage()
print(b, vp = viewport(width = 1, height = 1, x = 0.5, y = 0.5))
print(h, vp = viewport(width = 0.1, height = 0.2, x = 0.15, y = 0.15))
print(k, vp = viewport(width = 0.1, height = 0.2, x = 0.65, y = 0.4))
dev.off()

## ---- summary ----
for(i in 1:length(isos)){
  zebra_data_laea <- readRDS(
    paste0(workdir, "Data/zebra_data_laea_", isos[i], ".rds"))
  
  # Create Latex table from Info file
  (summary <- 
      data.frame(ID = unique(zebra_data_laea$individual.local.identifier), 
                 do.call("rbind", 
                         by(as.Date(timestamps(zebra_data_laea)), 
                            zebra_data_laea$individual.local.identifier, 
                            FUN= function(x) 
                              range(format(x, 
                                           format = "%Y-%m-%d", tz="GMT")))), 
                         unlist(lapply(timeLag(zebra_data_laea, 
                                               units='days'), sum)), 
                         unlist(lapply(distance(zebra_data_laea), sum))/1000, 
                         unlist(lapply(timeLag(zebra_data_laea, units='mins'), 
                                       function(x) round(median(x), 0)))))
  colnames(summary) <- c("ID", "Start date", "End date", 
                         "Duration (days)", "Distance (km)", 
                         "Time interval (min)")
  print(xtable(summary, 
caption=c(paste0("Start and end date of the tracking period, 
tracking duration (days), total distance (km) and 
median time interval (min) of each individual zebra."), 
          paste0("Tracking summary of each individual, ", region[i])), 
label=paste0("table:summary_zebra_move_", isos[i])), type="latex", 
        file = paste0("Tables/summary_zebra_move_", isos[i], ".tex"), 
caption.placement="top", 
booktabs=TRUE, include.rownames=FALSE, table.placement="H")
  rm(zebra_data_laea, summary)
}

## ---- studyarea ----

# Load countries data of Africa
data(countriesHigh, package="rworldxtra")
Africa <- subset(countriesHigh, continent=="Africa"); rm(countriesHigh)

# Load capital cities of Africa
data(world.cities, package="maps")
cities_afr <- world.cities[world.cities$country.etc %in% 
                             c(as.character(Africa$NAME)),]
cities_afr <- na.omit(cities_afr); rm(world.cities)
cap_cities_afr <- cities_afr[cities_afr$capital == 1,]

for(i in 1:length(isos)){
  # Load country boundaries
  gadm <- getData("GADM", country = getData("ISO3")[
    getData("ISO3")$ISO3 == isos[i],2], 
                  level=0, path=paste0(filedir, "GADM"))
  gadm_l1 <- getData("GADM", country = getData("ISO3")[
    getData("ISO3")$ISO3 == isos[i],2], 
                     level=1, path=paste0(filedir, "GADM"))
  
  # Extract capital city of country
  cap_city <- cap_cities_afr[cap_cities_afr$country.etc == 
                               getData("ISO3")[
                                 getData("ISO3")$ISO3 == isos[i],2],]
  cities <- cities_afr[cities_afr$country.etc == getData("ISO3")[
    getData("ISO3")$ISO3 == isos[i],2],]
  
  # Add symbol according to population size
  cities$symbol <- NA
  cities[cities$pop <= 500000,]$symbol <- 21
  if(i == 2){
    cities[cities$pop > 500000 & cities$pop <= 1000000,]$symbol <- 22
    cities[cities$pop > 1000000,]$symbol <- 23
  }
  # Get mean position of districts
  districts <- aggregate(cbind(long, lat) ~ id, data= fortify(gadm_l1), 
                         FUN=function(x)mean(range(x)))
  gadm_l1@data$id <- rownames(gadm_l1@data)
  districts <- plyr::join(districts, gadm_l1@data, by="id")
  districts$capital <- FALSE
  districts$capital[districts$NAME_1 == cap_city$name] <- TRUE
  
  # Load studyarea
  studyarea <- readOGR(dsn = paste0(workdir, "Results/"), 
                       layer = paste0("studyarea_", isos[i]), verbose=FALSE)
  
  # Plot of Africa
  a <- ggplot() + geom_polygon(data=Africa, 
                               aes(x=long, y=lat, group=group), 
                               fill="gray", colour=NA) + 
    geom_polygon(data=gadm, aes(x=long, y=lat, group=group), 
                 fill="grey40", colour="transparent") + 
    geom_polygon(data=studyarea, aes(x=long, y=lat, group=group), 
                 fill=NA, color="red") + 
    lims(x=c(-22.5,55), y=c(-40,40)) + 
    theme(panel.background = element_blank(), axis.ticks = element_blank(), 
          axis.title.x = element_blank(), axis.title.y = element_blank(), 
          axis.text.x = element_blank(), axis.text.y = element_blank(), 
          panel.border = element_rect(fill="transparent", colour = "black"), 
          legend.spacing=unit(-0.5, "lines")) + labs(x=NULL, y=NULL) + 
    coord_map()
  
  # Plot of Country
  b <- ggplot() + 
    scale_x_continuous(name=expression(paste("Longitude (",degree,")")), 
                       limits=c(floor(extent(gadm))[1], 
                                floor(extent(gadm))[2]), expand=c(0,0)) + 
    scale_y_continuous(name=expression(paste("Latitude (",degree,")")), 
                       limits=c(floor(extent(gadm))[3], 
                                floor(extent(gadm))[4]), expand=c(0,0)) + 
    geom_polygon(data=gadm_l1, aes(x=long, y=lat, group=group), 
                 fill="gray", colour="white", linetype="dashed") + 
    geom_polygon(data=gadm, aes(x=long, y=lat, group=group), 
                 fill=NA, colour="black") + 
    geom_polygon(data=studyarea, aes(x=long, y=lat, group=group), 
                 fill=NA, color="red") + 
    geom_text_repel(data=districts, aes(x=long, y=lat, 
                                        label = NAME_1, color= capital)) + 
    scale_color_manual(values =c('black', 'coral1'), guide="none") + 
    geom_point(data=cities, aes(x=long, y=lat), pch=cities$symbol, 
               color="black") + 
    geom_point(data = cap_city, aes(x = long, y = lat), 
               color = "coral1", pch=19) + 
    theme_bw() + theme(legend.position = "none") + coord_map()
  if(i == 1){
    # Map of Country with study area and location of country 
    # within Africa as Inlay
    cairo_pdf(paste0(workdir, "Figures/Studyarea_", isos[i], ".pdf"), 
              width=9, height=8, bg="transparent")
    grid.newpage()
    print(b, vp = viewport(width = 1, height = 1, x = 0.5, y = 0.5))
    print(a, vp = viewport(width = 0.2, height = 0.5, x = 0.85, y = 0.21))
    dev.off()
  }else if(i == 2){
    cairo_pdf(paste0(workdir, "Figures/Studyarea_", isos[i], ".pdf"), 
              width=7, height=8, bg="transparent")
    grid.newpage()
    print(b, vp = viewport(width = 1, height = 1, x = 0.5, y = 0.5))
    print(a, vp = viewport(width = 0.2, height = 0.5, x = 0.25, y = 0.21))
    dev.off()
  }
}

## ---- wdpa ----
for(i in 1:length(isos)){
  if(!file.exists(paste0(filedir, "WDPA/WDPA_July2016_", isos[i], 
                         "-shapefile"))){
   #Download shapefile of PAs
   download.file(paste0("http://d1gam3xoknrgr2.cloudfront.net/current/WDPA_", 
                        format(Sys.time(), "%B%Y"), "_" , isos[i], 
                        "-shapefile.zip"), 
                 destfile=paste0("WDPA_", format(Sys.time(), "%B%Y"), 
                                 "_", isos[i], ".zip"))
  
  # Unzip shapefile
  unzip(paste0("WDPA_", format(Sys.time(), "%B%Y"), "_", isos[i], ".zip"), 
        exdir=paste0(filedir, "WDPA/WDPA_", 
                     format(Sys.time(), "%B%Y"), "_", isos[i], "-shapefile"))
  }
  if(!file.exists(paste0(workdir, "Results/wdpa_sa_", isos[i], ".tif"))){
    # Read WDPA data
    wdpa <- readOGR(dsn=paste0(filedir, "WDPA/WDPA_July2016_", 
                               isos[i], "-shapefile"), layer = 
                      paste0("WDPA_July2016_", isos[i], 
                             "-shapefile-polygons"), 
                    verbose=FALSE)
    studyarea <- readOGR(dsn = "Results/", 
                         layer = paste0("studyarea_", isos[i]), 
                         verbose=FALSE)
    wdpa <- crop(wdpa, studyarea)
    writeOGR(wdpa, dsn= paste0(workdir, "Results/wdpa_sa_", isos[i]), 
             layer = "wdpa", driver="ESRI Shapefile", 
             check_exists=TRUE, overwrite_layer=TRUE)
  }
  
  # Load study area shapefile  
  studyarea_laea <- readOGR(dsn=paste0(workdir, "Results/"), 
                            layer=paste0("studyarea_laea_", isos[i]))
  
  # Load cropped WDPA
  if(!file.exists(paste0(workdir, "Results/wdpa_laea_", isos[i], ".tif"))){
    wdpa <- readOGR(dsn=paste0(workdir, "Results/wdpa_sa_", isos[i]), 
                    layer = "wdpa", verbose=FALSE)
    wdpa_laea <- spTransform(wdpa, crs(studyarea_laea))
    writeOGR(wdpa_laea, dsn = paste0(workdir, "Results/wdpa_laea_", isos[i]), 
             layer="wdpa", driver="ESRI Shapefile", 
             check_exists=TRUE, overwrite_layer=TRUE)
  } else{
    wdpa_laea <- readOGR(wdsn = paste0(workdir, 
                                       "Results/wdpa_laea_", isos[i]), 
                         layer="wdpa", verbose=FALSE)
  }
  
  # Create area with 5 km buffer around PAs
  wdpa_buf_5km <- gBuffer(wdpa_laea, width=5000, byid=TRUE)
  wdpa_buf_5km <- gIntersection(wdpa_buf_5km, studyarea_laea)
  
  # Create shapefile of unprotected areas
  unp <- gDifference(studyarea_laea, wdpa_buf_5km)
  
  # Just select buffer
  wdpa_buf_5km <- gDifference(wdpa_buf_5km, wdpa_laea)
  wdpa_buf_5km <- as(wdpa_buf_5km, "SpatialPolygonsDataFrame")
  
  # Save 5km buffer to file
  writeOGR(wdpa_buf_5km, dsn = paste0(workdir, 
                                      "Results/buf_5km_laea_", isos[i]), 
           layer="buf_5km", driver="ESRI Shapefile", 
           check_exists=TRUE, overwrite_layer=TRUE)
  
  # Calculate percentage of area, which is protected and 
  # which lies in buffer zone
  gArea(wdpa_laea)/gArea(studyarea_laea)*100
  gArea(wdpa_buf_5km)/gArea(studyarea_laea)*100
  
  if(file.exists(paste0(workdir, "Results/wdpa_30m_", isos[i], ".tif"))){
    wdpa_30m <- raster(paste0(workdir, "Results/wdpa_30m_", 
                              isos[i], ".tif"))
  }else{
    # Rasterize protected area shapefile
    writeRaster(raster(resolution = c(30.7, 29.3), 
                       extent(studyarea_laea), vals=0, crs=crs.laea), 
                paste0(workdir, "/Results/wdpa_30m_", isos[i], ".tif"), 
                overwrite=TRUE)
    wdpa_30m <- gdal_rasterize(paste0(workdir, "Results/wdpa_laea_", 
                                      isos[i]), 
                               paste0(workdir, "Results/wdpa_30m_", 
                                      isos[i], ".tif"), 
                               burn=1, l="wdpa", output_Raster=TRUE)
  }
  
  if(file.exists(paste0(workdir, "Results/buf_5km_30m_", isos[i], ".tif"))){
    buf_5km_30m <- raster(paste0(workdir, "Results/buf_5km_30m_", 
                                 isos[i], ".tif"))
  }else{
    # Rasterize 5km buffer shapefile
    writeRaster(raster(resolution = c(30.7, 29.3), 
                       extent(studyarea_laea), vals=0, crs=crs.laea), 
                paste0(workdir, "/Results/buf_5km_30m_", isos[i], ".tif"), 
                overwrite=TRUE)
    buf_5km_30m <- gdal_rasterize(paste0(workdir, "Results/buf_5km_laea_", 
                                         isos[i]), 
                                  paste0(workdir, "Results/buf_5km_30m_", 
                                         isos[i], ".tif"), 
                                  burn=2, l="buf_5km", output_Raster=TRUE)
  }
  
  # Read zebra movement data
  zebra_data_laea <- readRDS(paste0(workdir, "Data/zebra_data_laea_", 
                                    isos[i], ".rds"))
  if(class(zebra_data_laea) != "MoveStack"){
    zebra_data_laea <- moveStack(zebra_data_laea)
  }
  
  # Combine wpda and buf
  wdpa_buf <- mosaic(wdpa_30m, buf_5km_30m, fun=sum)
  
  # Extract wdpa data for zebra locations
  zebra_data_laea$status <- raster::extract(wdpa_buf, zebra_data_laea)
  
  # Save zebra data
  saveRDS(zebra_data_laea, paste0(workdir, 
                                  "Data/zebra_data_laea_", isos[i],
                                  ".rds"))
  rm(wdpa_buf)

  if(length(wdpa[wdpa$DESIG == "Ramsar Site, 
                 Wetland of International Importance",]) > 0){
    wdpa$DESIG <- as.character(wdpa$DESIG)
    wdpa$DESIG[wdpa$DESIG == "Ramsar Site, 
               Wetland of International Importance"] <- "Ramsar Site"
  }
  
  wdpa$STATUS <- as.character(wdpa$STATUS)
  wdpa$STATUS[wdpa$STATUS == "Not Applicable"] <- "NA"
  wdpa$STATUS[wdpa$STATUS == "Not Reported"] <- "-"
  
  # Create a summary table of protected areas and save it to file
  summary_wdpa <- wdpa@data[c("NAME", "DESIG", "IUCN_CAT", 
                              "REP_AREA", "STATUS_YR")]
  colnames(summary_wdpa) <- c("Name", "Designation", "IUCN Category", 
                              "Area (km2)", "Year")
  print(xtable(summary_wdpa, caption=c("Summary of the Protected 
                                       Areas of the study area.", 
                                       paste0("Summary of the 
                                              Protected Areas, ", 
                                              region[i])), 
               label=paste0("table:summary_wdpa_", isos[i]), digits=0), 
        type="latex", file = paste0("Tables/summary_wdpa_", isos[i], ".tex"), 
        caption.placement="top", booktabs=TRUE, table.placement="H", 
        include.rownames=FALSE) 

  # Transform buf to crs.wgs84
  buf_5km_wgs84 <- spTransform(wdpa_buf_5km, crs.wgs84)
  
  # Convert protect areas shapefile into a dataframe
  wdpa@data$id <- rownames(wdpa@data)
  df_wdpa <- fortify(wdpa, region="id")
  df_wdpa <- plyr::join(df_wdpa, wdpa@data, by="id")
  if(i == 1){
    # Map of protected areas, colour coded by year of implementation
    ggplot() + geom_polygon(data=studyarea, aes(x=long, y=lat, group=group), 
                            fill="snow2", color=NA) + 
      geom_polygon(data=buf_5km_wgs84, aes(x=long, y=lat, group=group), 
                   fill="darkseagreen2", color="darkseagreen2") + 
      geom_polygon(data=df_wdpa, aes(x=long, y=lat, group=group, 
                                     fill=factor(STATUS_YR)), 
                   colour="forestgreen") + 
      geom_polygon(data=subset(df_wdpa, NAME=="Maun"), 
                   aes(x=long, y=lat, group=group, fill=factor(STATUS_YR)), 
                   colour="forestgreen") + 
      geom_polygon(data=subset(df_wdpa, NAME=="Okavango Delta"), 
                   aes(x=long, y=lat, group=group, fill=factor(STATUS_YR))) + 
      geom_polygon(data=subset(df_wdpa, NAME=="Moremi"), 
                   aes(x=long, y=lat, group=group, fill=factor(STATUS_YR)), 
                   colour="forestgreen") + 
      geom_polygon(data=studyarea, aes(x=long, y=lat, group=group), 
                   fill="transparent", color="red") + 
      labs(x=expression(paste("Longitude (",degree,")")), 
          y=expression(paste("Latitude (",degree,")"))) + 
      scale_fill_manual(values=viridis::viridis(nlevels(
        factor(df_wdpa$STATUS_YR))), 
                        name="Status year") + coord_map() + 
      theme(legend.position = "bottom")
    ggsave(paste0("Figures/Temp_WDPA_", isos[i], ".pdf"), 
           width=9, height=8)
  } else if(i == 2){
    # Map of protected areas, colour coded by year of implementation
    p1<- ggplot() + geom_polygon(data=studyarea, 
                                 aes(x=long, y=lat, group=group), 
                            fill="snow2", color=NA) + 
      geom_polygon(data=buf_5km_wgs84, aes(x=long, y=lat, group=group), 
                   fill="darkseagreen2", color="darkseagreen2") + 
      geom_polygon(data=df_wdpa, 
                   aes(x=long, y=lat, group=group, fill=factor(STATUS_YR)), 
                   colour="forestgreen") + 
      geom_polygon(data=studyarea, aes(x=long, y=lat, group=group), 
                   fill="transparent", color="red") + 
      xlab(expression(paste("Longitude (",degree,")"))) + 
      ylab(expression(paste("Latitude (",degree,")"))) + 
      scale_fill_manual(
        values=viridis::viridis(nlevels(factor(df_wdpa$STATUS_YR))), 
                        name="Status year") + coord_map() + 
      theme(legend.position = "bottom") + 
      guides(fill=guide_legend(nrow=2,byrow=TRUE))
  }
  
  # Only take the first entry of each list and add it to labels!!!
  labels_wdpa <- strsplit(levels(factor(df_wdpa$DESIG)), split=",")
  labels_wdpa <- c(unlist(lapply(labels_wdpa, function(x) x[[1]][1])))

  # Map of country with protected areas according to protection status
  if(i == 1){
    ggplot() + geom_polygon(data=studyarea, aes(x=long, y=lat, group=group), 
                            fill="snow2", color=NA) + 
      geom_polygon(data=buf_5km_wgs84, aes(x=long, y=lat, group=group), 
                            fill="darkseagreen2", color="darkseagreen2") + 
      geom_polygon(data=df_wdpa, aes(x=long, y=lat, 
                                     group=group, fill=DESIG), 
                            colour="forestgreen") + 
      geom_polygon(data=subset(df_wdpa, NAME=="Maun"), 
                   aes(x=long, y=lat, group=group, fill=DESIG), 
                   colour="forestgreen") + 
      geom_polygon(data=subset(df_wdpa, NAME=="Okavango Delta"), 
                   aes(x=long, y=lat, group=group, fill=DESIG), 
                   colour="forestgreen") + 
      geom_polygon(data=subset(df_wdpa, NAME=="Moremi"), 
                   aes(x=long, y=lat, group=group, fill=DESIG), 
                   colour="forestgreen") + 
      geom_polygon(data=studyarea, aes(x=long,y=lat,group=group), 
                   fill="transparent", colour="red") + 
      xlab(expression(paste("Longitude (",degree,")"))) + 
      ylab(expression(paste("Latitude (",degree,")"))) + 
      scale_fill_manual(
        values=viridis::viridis(nlevels(factor(df_wdpa$DESIG))), 
                        name="Designation", labels=
          levels(factor(labels_wdpa))) + 
      coord_map() + theme(legend.position = "bottom")
    ggsave(paste0("Figures/Status_WDPA_", isos[i], ".pdf"), 
           width=9, height=8)
  } else if(i == 2){
    p2 <- ggplot() + geom_polygon(data=studyarea, aes(x=long, y=lat, 
                                                      group=group), 
                          fill="snow2", color=NA) + 
      geom_polygon(data=buf_5km_wgs84, aes(x=long, y=lat, group=group), 
                            fill="darkseagreen2", color="darkseagreen2") + 
      geom_polygon(data=df_wdpa, aes(x=long, y=lat, group=group, 
                                     fill=DESIG), 
                            colour="forestgreen") + 
      geom_polygon(data=studyarea, aes(x=long,y=lat,group=group), 
                   fill="transparent", colour="red") + 
      xlab(expression(paste("Longitude (",degree,")"))) + 
      ylab(expression(paste("Latitude (",degree,")"))) + 
      scale_fill_manual(values=
                          viridis::viridis(nlevels(factor(df_wdpa$DESIG))), 
                        name="Designation", labels=
                          levels(factor(labels_wdpa))) + 
      coord_map() + theme(legend.position = "bottom") + 
      guides(fill=guide_legend(nrow=2,byrow=TRUE))
    cowplot::plot_grid(p1, p2, ncol=2, labels="auto")
    ggsave(paste0("Figures/Status_Temp_WDPA_", isos[i], ".pdf"), 
           width=9, height=8)
  }
}

## ---- worldclim ----
for(i in 1:length(isos)){
  # Load study area shapefile
  studyarea <- readOGR(dsn = "Results/", layer = 
                         paste0("studyarea_", isos[i]))
  
  # Read cropped WDPA shapefile
  wdpa <- readOGR(dsn=paste0(workdir, "Results/wdpa_sa_", isos[i]), 
                  layer = "wdpa", verbose=FALSE)
  
  # Get climate & elevation data
  if(i == 1){
    prec <- raster::getData("worldclim", var="prec", res=.5,
                            lon=coordinates(studyarea)[,1],
                            lat=coordinates(studyarea)[,2], 
                            path=paste0(filedir, "worldclim"), download=T)
    tmean <- raster::getData("worldclim", var="tmean", res=.5,
                             lon=coordinates(studyarea)[,1],
                             lat=coordinates(studyarea)[,2], 
                             path=paste0(filedir, "worldclim"), download=T)
  } else{
    prec <- merge(raster::getData("worldclim", var="prec", res=.5,
                                  lon=coordinates(studyarea)[,1],
                                  lat=coordinates(studyarea)[,2], 
                                  path=paste0(filedir, "worldclim"), 
                                  download=T), 
                  raster::getData("worldclim", var="prec", res=.5,
                                  lon=coordinates(studyarea)[,1],
                                  lat=coordinates(studyarea)[,2]-5, 
                                  path=paste0(filedir, "worldclim"),
                                  download=T))
    tmean <- merge(raster::getData("worldclim", var="tmean", res=.5,
                                   lon=coordinates(studyarea)[,1],
                                   lat=coordinates(studyarea)[,2], 
                                   path=paste0(filedir, "worldclim"), 
                                   download=T), 
                   raster::getData("worldclim", var="tmean", res=.5,
                                   lon=coordinates(studyarea)[,1],
                                   lat=coordinates(studyarea)[,2]-5, 
                                   path=paste0(filedir, "worldclim"), 
                                   download=T))
  }
  
  # Mask and crop precipitation data by country
  prec <- mask(crop(prec, studyarea), studyarea)
  tmean <- mask(crop(tmean,studyarea), studyarea)
  
  # Standardise temperature to degree Celsius
  tmean <- calc(tmean, fun=function(x){x/10})
  
  # Monthly precipitation
  prec_df <- as.data.frame(rasterToPoints(prec))
  prec_df <- tidyr::gather(prec_df, month, prec, -c(x,y))
  prec_df$month <- factor(prec_df$month, 
                          labels = c("January", "February", "March", 
                                     "April", "May", "June", "July", 
                                     "August", "September","October", 
                                     "November","December"))
  
  # Spatial average mean precipitation & temperature data
  prec_sp_avg <- cellStats(prec, stat="mean")
  tmean_sp_avg <- cellStats(tmean, stat="mean")
  clim_sp_avg <- as.data.frame(cbind(prec_sp_avg, tmean_sp_avg))
  
  # Mean prec and temp
  mprec <- fortify(overlay(prec, fun=mean), maxpixels=maxpixels)
  mtemp <- fortify(overlay(tmean, fun=mean), maxpixels=maxpixels)
  
  # Map of mean annual temperature and annual precipitation (mm)
  p1 <- ggplot() + geom_raster(data=mprec, aes(x,y, fill = layer)) + 
    scale_fill_gradient(name="Prec (mm)", low = "white", high = "blue", 
                        na.value = "transparent", guide = "colourbar") + 
    labs(x=expression(paste("Longitude (",degree,")")), 
         y=expression(paste("Latitude (",degree,")"))) + 
    geom_polygon(data=wdpa, aes(x=long, y=lat, group=group), 
                 fill="transparent", color="green3") + 
    geom_polygon(data=studyarea, aes(x=long, y=lat, group=group), 
                 fill="transparent", color="red")
  p2 <- ggplot() + geom_raster(data=mtemp, aes(x,y, fill = layer)) + 
    scale_fill_gradient2(name=expression(paste("Temp (",degree,"C)")), 
                         low="blue", mid="white", high="red", 
                         midpoint=mean(mtemp$layer), na.value="transparent") + 
    labs(x=expression(paste("Longitude (",degree,")")), 
         y=expression(paste("Latitude (",degree,")"))) + 
    geom_polygon(data=wdpa, aes(x=long, y=lat, group=group), 
                 fill="transparent", color="green3") + 
    geom_polygon(data=studyarea, aes(x=long, y=lat, group=group), 
                 fill="transparent", color="red")
  cowplot::plot_grid(p1, p2, ncol=2, labels = "auto")
  ggsave(paste0("Figures/Env_", isos[i], ".pdf"), width=10, height=4)
  
  # Monthly precipitation
  if(i == 1){
    ggplot() + geom_raster(data=prec_df, aes(x,y, fill=prec)) + 
      facet_wrap(~ month) + 
      scale_fill_gradient(name="Precipitation (mm)", low = "white", 
                          high = "blue", 
                          na.value = "transparent", guide = "colourbar") + 
      labs(x=expression(paste("Longitude (",degree,")")), 
           y=expression(paste("Latitude (",degree,")"))) + 
      geom_polygon(data=wdpa, aes(x=long, y=lat, group=group), 
                   fill="transparent", color="green3") + 
      geom_polygon(data=studyarea, aes(x=long, y=lat, group=group), 
                   fill="transparent", color="red") + 
      theme(legend.position = "bottom")
    ggsave(paste0("Figures/Monthly_Prec_", isos[i], ".pdf"), 
           width=8, height=6)
  } else if(i == 2){
    ggplot() + geom_raster(data=prec_df, aes(x,y, fill=prec)) + 
      facet_wrap(~ month) + 
      scale_fill_gradient(name="Precipitation (mm)", low = "white", 
                          high = "blue", 
                          na.value = "transparent", guide = "colourbar") + 
      labs(x=expression(paste("Longitude (",degree,")")), 
           y=expression(paste("Latitude (",degree,")"))) + 
      geom_polygon(data=wdpa, aes(x=long, y=lat, group=group),
                   fill="transparent", color="green3") + 
      geom_polygon(data=studyarea, aes(x=long, y=lat, group=group), 
                   fill="transparent", color="red")
    ggsave(paste0("Figures/Monthly_Prec_", isos[i], ".pdf"), 
           width=8, height=8)
  }

  # Plot average precipitation and mean temperature data
  cairo_pdf(paste0(workdir, "Figures/Climograph_", isos[i], ".pdf"), 
            width=10, height=6)
  par(mar = c(4.5,4.5,1.5,4.5)+.1, bg="white", cex.lab=1.5, cex.axis=1.25)
  barplot(clim_sp_avg$prec_sp_avg, space=0.5, 
          ylim=c(0,max(ceiling(clim_sp_avg$prec_sp_avg/20)*20)), 
          pch=19, axes=FALSE, 
          col="blue", xlab="Month of the Year", ylab="")
  axis(side=2, at=seq(0, max(ceiling(clim_sp_avg$prec_sp_avg/20)*20), by=20), 
       pos= -0.2)
  axis(side=1, at=c(1,2.5,4,5.5,7,8.5,10,11.5,13,14.5,16,17.5), 
       labels=c("J","F","M","A","M","J","J","A","S","O","N","D"))
  lines((clim_sp_avg$tmean_sp_avg*2)~c(1,2.5,4,5.5,7,8.5,
                                       10,11.5,13,14.5,16,17.5), 
        col="red", lwd=2.5)
  axis(side=4, at=seq(0, max(ceiling(clim_sp_avg$tmean_sp_avg/10)*20), 
                      by=20), 
       labels=seq(0, max(ceiling(clim_sp_avg$tmean_sp_avg/10)*10), by=10))
  mtext(side=2, line=3, "Precipitation (mm)", col="blue", cex=1.5)
  mtext(side =4, line=3, expression(paste("Temperature (",degree,"C)")), 
        cex=1.5, col="red")
  box()
  dev.off()
}

## ---- fence ----

# Read fence shapefiles
fence <- readOGR(paste0(workdir, "Data/Fence.kml"), 
                 "Fence", verbose=FALSE, pointDropZ=TRUE)
cattle_ranches <- readOGR(paste0(workdir, "Data/CattleRanches.kml"), 
                          "CattleRanches", verbose=FALSE, pointDropZ=TRUE)

# Merge shapefiles
fences <- gUnion(fence, cattle_ranches); rm(fence, cattle_ranches)
fences_laea <- spTransform(fences, crs.laea); rm(fences)

# Read study area shapefile 
studyarea_laea <- readOGR(dsn = paste0(workdir, "Results/"), 
                          layer = paste0("studyarea_laea_", isos[1]), 
                          verbose=FALSE)

# Intersect fence with studyarea
fences <- as(gIntersection(fences_laea, studyarea_laea), 
             "SpatialPolygonsDataFrame")

# Save shapefile
writeOGR(fences, dsn= paste0(workdir, "Results/"), 
         layer = "fences", driver="ESRI Shapefile", 
         check_exists=TRUE, overwrite_layer=TRUE)
fences <- readOGR(paste0(workdir, "Results/"), layer="fences")

# Create to empty raster with 30 m resolution
writeRaster(raster(resolution = c(30.7, 29.3), 
                   extent(studyarea_laea), vals=0, 
                   crs=crs.laea), 
            paste0(workdir, "/Results/fences_30m.tif"), 
            overwrite=TRUE)

# Read location of shp and tif file
fences_30m <- gdal_rasterize(paste0(workdir ,"Results/fences.shp"), 
                             paste0(workdir, "Results/fences_30m.tif"), 
                             burn=1, l="fences", output_Raster=TRUE)

# Read zebra movement data
zebra_data_laea <- readRDS(paste0(workdir, 
                                  "Data/zebra_data_laea_", 
                                  isos[1], ".rds"))

# Extract fence data for zebra locations
zebra_data_laea$fences <- raster::extract(fences_30m, zebra_data_laea)

# Save zebra data
saveRDS(zebra_data_laea, 
        paste0(workdir, 
               "Data/zebra_data_laea_", 
               isos[1], ".rds"))

# Plot Fence
fence_df <- fortify(fences_30m, maxpixels=maxpixels); rm(fences_30m)
names(fence_df) <- c("x","y","fence")
ggplot() + geom_raster(data=fence_df, aes(x, y, fill = factor(fence))) + 
  scale_fill_manual(name="Fenced", na.value="transparent", 
                    values=viridis::viridis(2), breaks=c(0, 1), 
                    labels=c("No", "Yes")) + 
  geom_polygon(data=studyarea_laea, aes(x=long, y=lat, group=group), 
               fill="transparent", color="red") + 
  theme(legend.background = element_rect(fill="transparent", colour=NA), 
        axis.text.y = rotatedAxisElementText(90,"y")) + 
  labs(x="Eastings (m)", y="Northings (m)") + coord_equal(expand=FALSE)
ggsave(paste0("Figures/fence_", isos[1], ".pdf"), 
       width=6, height=4); rm(fence_df)

## ---- move_plot ----
for(i in 1:2){
  # Read study area shapefile 
  studyarea_laea <- readOGR(dsn = "Results/", 
                            layer = paste0("studyarea_laea_", 
                                           isos[i]))
  
  # Read cropped WDPA shapefile
  wdpa_laea <- readOGR(dsn=paste0(workdir, "Results/wdpa_laea_", 
                                  isos[i]), 
                       layer = "wdpa", verbose=FALSE)
  
  # Read zebra movement data
  zebra_data_laea <- readRDS(paste0(workdir, 
                                    "Data/zebra_data_laea_", 
                                    isos[i], ".rds"))
  if(class(zebra_data_laea) != "MoveStack"){
    zebra_data_laea <- moveStack(zebra_data_laea)
  }
  
  # Get mean position of pas
  pa_names <- aggregate(cbind(long, lat) ~ id, 
                        data=fortify(wdpa_laea), FUN=function(x) 
                          mean(range(x)))
  wdpa_laea@data$id <- rownames(wdpa_laea@data)
  pa_names <- plyr::join(pa_names, wdpa_laea@data, by="id")
  
  df_wdpa <- fortify(wdpa_laea, region="id")
  df_wdpa <- plyr::join(df_wdpa, wdpa_laea@data, by="id")
  
  # Extract start point & end point for each individual
  zebra_start <- do.call("rbind", 
                         lapply(move::split(zebra_data_laea), function(x) 
                           data.frame(x[1])))
  zebra_end <- do.call("rbind", 
                       lapply(move::split(zebra_data_laea), function(x) 
                         data.frame(x[length(x)])))
  
  # Create plot
  if(i == 1){
    fences_laea <- readOGR(paste0(workdir, "Results/"), layer="fences")
    a <- ggplot(data=studyarea_laea, aes(x=long, y=lat)) + 
      geom_polygon(aes(group=group), fill="snow2", colour=NA) + 
      geom_polygon(data=wdpa_laea, aes(x=long, y=lat, group=group), 
                   fill="darkseagreen2", colour="forestgreen") + 
      geom_polygon(data=subset(df_wdpa, NAME=="Maun"), 
                   aes(x=long, y=lat, group=group), 
                   fill="darkseagreen2", colour="forestgreen") + 
      geom_polygon(data=subset(df_wdpa, NAME=="Okavango Delta"), 
                   aes(x=long, y=lat, group=group), 
                   fill="darkseagreen2", colour="forestgreen") + 
      geom_polygon(data=subset(df_wdpa, NAME=="Moremi"), 
                   aes(x=long, y=lat, group=group), 
                   fill="darkseagreen2", colour="forestgreen") + 
      geom_polygon(data=fences_laea, aes(x=long, y=lat, group=group), 
                     fill=NA, colour="brown", linetype="dashed") + 
      geom_polygon(data=studyarea_laea, aes(x=long, y=lat, group=group), 
                   fill=NA, colour="red")
  } else if(i == 2){
    a <- ggplot(data=studyarea_laea, aes(x=long, y=lat)) + 
      geom_polygon(aes(group=group), fill="snow2", colour=NA) + 
      geom_polygon(data=wdpa_laea, aes(x=long, y=lat, group=group), 
                   fill="darkseagreen2", colour="forestgreen") + 
      geom_polygon(data=studyarea_laea, aes(x=long, y=lat, group=group), 
                   fill=NA, colour="red")
  }
  a <- a + geom_point(data = data.frame(zebra_data_laea), 
               aes(x = coords.x1, y = coords.x2, 
                   colour=individual.local.identifier), 
               shape=1, size=0.7) + 
    geom_path(aes(x = coords.x1, y = coords.x2, 
                  colour=individual.local.identifier), 
              data = data.frame(zebra_data_laea)) + 
    geom_point(aes(x = coords.x1, y = coords.x2, shape="start"), size=1, 
               data = zebra_start, colour="black") + 
    geom_point(aes(x = coords.x1, y = coords.x2, shape="end"), size=1, 
               data = zebra_end, colour="black") + 
    scale_x_continuous(name="Eastings (m)", expand=c(0,0)) + 
    scale_y_continuous(name="Northings (m)", expand=c(0,0)) + 
    scale_colour_discrete(name="ID") + 
    scale_shape_manual(name="", values=c(start=2, end=0), 
                       labels=c(start="Start point", end="End point")) + 
    guides(color = guide_legend(order=1), shape = guide_legend(order=2)) + 
    coord_equal(expand=FALSE)
  
  # Need to specify coordinate names, placement location and unit
  if(i == 1){
    a <- a + 
      geom_text_repel(data=pa_names, aes(x=long, y=lat, label = NAME), 
                      color = "darkgreen", size=4) + 
      theme(legend.justification=c(1,1), legend.position=c(1,1), 
            legend.background = element_rect(fill="transparent", colour=NA), 
            axis.text.y = rotatedAxisElementText(90,"y"))
    a <- ggScaleBar(a, extent(studyarea_laea), 
                    addParams = list(noBins = 3, placement="bottomleft", 
                                     unit="km", 
                                     sbLengthPct = 0.2, 
                                     sbHeightvsWidth = 1/40, 
                                     size=4, margin=7000))
    cairo_pdf(paste0("Figures/Zebra_", isos[i], "_Overview.pdf"), 
              width=6, height=5, bg="transparent")
    grid.newpage()
    print(a, vp = viewport(width = 1, height = 1, x = 0.5, y = 0.5))
    print(h, vp = viewport(width = 0.1, height = 0.2, x = 0.19, y = 0.25))
    dev.off()
  } else if(i == 2){
    a <- a + 
      geom_text_repel(data=pa_names, aes(x=long, y=lat, label = NAME), 
                      color = "darkgreen", size=4) + 
      theme(legend.background = element_rect(fill="transparent", colour=NA), 
            axis.text.y = rotatedAxisElementText(90,"y"))
    a <- ggScaleBar(a, extent(studyarea_laea), 
                    addParams = list(noBins = 3, placement="bottomleft", 
                                     unit="km", 
                                     sbLengthPct = 0.2, 
                                     sbHeightvsWidth = 1/20, 
                                     size=4, margin=3500))
    cairo_pdf(paste0("Figures/Zebra_", isos[i], "_Overview.pdf"), 
              width=5, height=6, bg="transparent")
    grid.newpage()
    print(a, vp = viewport(width = 1, height = 1, x = 0.5, y = 0.5))
    print(h, vp = viewport(width = 0.1, height = 0.2, x = 0.23, y = 0.19))
    dev.off()
  }
  
  # Zebra map facetted by individual
  a <- ggplot(data=studyarea_laea, aes(x=long, y=lat)) + 
    geom_polygon(aes(group=group), fill="snow2", colour=NA) + 
    geom_polygon(data=wdpa_laea, aes(x=long, y=lat, group=group), 
                 fill="darkseagreen2", colour="forestgreen") + 
    geom_polygon(data=studyarea_laea, aes(x=long, y=lat, group=group), 
                 fill=NA, colour="red") + 
    geom_point(aes(x = coords.x1, y = coords.x2, 
                   colour=individual.local.identifier), 
               shape=1, size=0.7, data = data.frame(zebra_data_laea)) + 
    geom_path(aes(x = coords.x1, y = coords.x2, 
                  colour=individual.local.identifier), 
              data = data.frame(zebra_data_laea)) + 
    geom_point(aes(x = coords.x1, y = coords.x2, shape="start"), size=1, 
               data = zebra_start, colour="black") + 
    geom_point(aes(x = coords.x1, y = coords.x2, shape="end"), size=1, 
               data = zebra_end, colour="black") + 
    scale_x_continuous(name="Eastings (m)", expand=c(0,0)) + 
    scale_y_continuous(name="Northings (m)", expand=c(0,0)) + 
    coord_equal(expand=FALSE)
  if(i == 1){
    a + facet_wrap(~ individual.local.identifier, ncol=3) + 
      theme(legend.position = "none")
    ggsave(paste0("Figures/Zebra_", isos[i], "_Individual.pdf"), 
           width=6, height=5, bg="transparent")
  } else if(i == 2){
    a + facet_wrap(~ individual.local.identifier, ncol=5) + 
      theme(legend.position = "none")
    ggsave(paste0("Figures/Zebra_", isos[i], "_Individual.pdf"), 
           width=6, height=6, bg="transparent")
  }

  # Create monthly vector and plot
  zebra_data_laea$month <- factor(month(timestamps(zebra_data_laea)), 
                                  labels=c("January", "February", "March", 
                                           "April", "May", "June", "July", 
                                           "August", "September", "October", 
                                           "November", "December"))
  a <- ggplot() + geom_polygon(data=studyarea_laea, aes(x=long, y=lat, 
                                                        group=group), 
                          fill="snow2", colour=NA) + 
    geom_polygon(data=wdpa_laea, aes(x=long, y=lat, group=group), 
                 fill="darkseagreen2", colour="forestgreen") + 
    geom_polygon(data=studyarea_laea, aes(x=long, y=lat, group=group), 
                 fill="transparent", colour="red", lwd=1.5) + 
    geom_point(aes(x = coords.x1, y = coords.x2, 
                   colour=individual.local.identifier), 
               shape=1, size=0.7, data = data.frame(zebra_data_laea)) + 
    geom_path(aes(x = coords.x1, y = coords.x2, 
                  colour=individual.local.identifier), 
              data = data.frame(zebra_data_laea))+ 
    facet_wrap(~ month, ncol=4) + 
    scale_x_continuous(name="Eastings (m)", expand=c(0,0)) + 
    scale_y_continuous(name="Northings (m)", expand=c(0,0)) + 
    scale_colour_discrete(name="ID") + coord_equal(expand=FALSE)
  if(i == 1){
    a
    ggsave(paste0("Figures/Zebra_Month_", isos[i], ".pdf"), 
           width=6, height=4, bg="transparent")
  } else if (i == 2){
    a
    ggsave(paste0("Figures/Zebra_Month_", isos[i], ".pdf"), 
           width=6, height=6, bg="transparent")
  }
}

## ---- segmentation ----

# Read zebra movement data
zebra_data_laea <- readRDS(paste0(workdir, "Data/zebra_data_laea_", 
                                  isos[1], ".rds"))

# Mean starting location is mean location between oct 20-31, 
# here we look at population
(start2007 <- with(
  data.frame(subset(zebra_data_laea, 
                    timestamps(zebra_data_laea) >= "2007-10-20" & 
                      timestamps(zebra_data_laea) < "2007-11-01")), 
                   mean(complex(re = coords.x1, im =  coords.x2))))
(start2008 <- with(
  data.frame(subset(zebra_data_laea, 
                    timestamps(zebra_data_laea) >= "2008-10-20" & 
                      timestamps(zebra_data_laea) < "2008-11-01")), 
                   mean(complex(re =  coords.x1, im =  coords.x2))))
startsouth <- c(start2007, start2008)

# Mean ending location is mean location between nov 25 - nov 30
(end2007 <- with(data.frame(
  subset(zebra_data_laea, timestamps(zebra_data_laea) >= "2007-11-25" & 
           timestamps(zebra_data_laea) < "2007-12-01")), 
  mean(complex(re =  coords.x1, im = coords.x2))))
(end2008 <- with(data.frame(
  subset(zebra_data_laea, timestamps(zebra_data_laea) >= "2008-11-25" & 
           timestamps(zebra_data_laea) < "2008-12-01")), 
  mean(complex(re =  coords.x1, im =  coords.x2))))
endsouth <- c(end2007, end2008)

getZebSouthMigList = function(year = c("2007", "2008"))
{year = match.arg(year)
if (year == "2007"){
  ids = 
    unique(zebra_data_laea$individual.local.identifier[
      timestamps(zebra_data_laea) < "2007-11-01"])
  zebmig = subset(zebra_data_laea, 
                  study.local.timestamp >= "2007-10-20" & 
                    study.local.timestamp < "2007-12-01" & 
                    individual.local.identifier %in% ids)
  start <- start2007
  end <- end2007
} else if (year == "2008"){
  ids = unique(zebra_data_laea$individual.local.identifier[
    timestamps(zebra_data_laea) >= "2008-10-24"])
  zebmig = subset(zebra_data_laea, 
                  study.local.timestamp >= "2008-10-20" & 
                    study.local.timestamp < "2008-12-01" & 
                    individual.local.identifier %in% ids)
  start <- start2008
  end <- end2008
}

zebmig$z = complex(re = zebmig$coords.x1, im = zebmig$coords.x2)
zebmig$distsStart = Mod(start - zebmig$z)
zebmig$distsEnd = Mod(end - zebmig$z)
zebmig$angle = Arg(zebmig$z - start) # angle from start to this point

zebmig = split(data.frame(zebmig), zebmig$individual.local.identifier)
return(zebmig)
}

getSouthMigDates = function()
{
  zeb2008mig = getZebSouthMigList(year = "2008")
  zeb2007mig = getZebSouthMigList(year = "2007")
  
  # criteria: distance from start > 25 km
  (migrStartIdxs2008 = sapply(zeb2008mig, with, 
                              which(distsStart > 25000)[1]))
  (migrStartIdxs2007 = sapply(zeb2007mig, with, 
                              which(distsStart > 25000)[1]))
  
  # could then trace back to point closest to the start pre-departure 
  migrStartDates2008 = lapply(1:length(zeb2008mig), 
                              function(x) 
                                zeb2008mig[[x]]$study.local.timestamp[
                                  migrStartIdxs2008[x]])
  migrStartDates2007 = lapply(1:length(zeb2007mig), 
                              function(x) 
                                zeb2007mig[[x]]$study.local.timestamp[
                                  migrStartIdxs2007[x]])
  
  # criteria: distance from end < 5 km or min distance 
  # (second zebra gets within 6.3 km)
  (migrEndIdxs2008 = sapply(zeb2008mig, with, 
                            min(which(distsEnd < 5000)[1], 
                                which(distsEnd == min(distsEnd)), 
                                na.rm = TRUE)))
  (migrEndIdxs2007 = sapply(zeb2007mig, with, 
                            min(which(distsEnd < 5000)[1], 
                                which(distsEnd == min(distsEnd)), 
                                na.rm = TRUE)))
  
  migrEndDates2008 = lapply(1:length(zeb2008mig), 
                            function(x) 
                              zeb2008mig[[x]]$study.local.timestamp[
                                migrEndIdxs2008[x]])
  migrEndDates2007 = lapply(1:length(zeb2007mig), 
                            function(x) 
                              zeb2007mig[[x]]$study.local.timestamp[
                                migrEndIdxs2007[x]])
  
  # need to specify timezone or else becomes current computer timezone, 
  # using UTC as defaut, though really CAS
  zebDates <- 
    data.frame(id = c(sapply(zeb2007mig, with, 
                             individual.local.identifier[1]), 
                      sapply(zeb2008mig, with, 
                             individual.local.identifier[1])), 
                        startIdx = c(migrStartIdxs2007, migrStartIdxs2008), 
                        startDate = 
                 as.POSIXct(do.call(c, c(migrStartDates2007, 
                                         migrStartDates2008)), tz = "UTC"), 
                        endIdx = c(migrEndIdxs2007, migrEndIdxs2008), 
                        endDate = as.POSIXct(do.call(c, 
                                                     c(migrEndDates2007, 
                                                       migrEndDates2008)), 
                                             tz = "UTC"), 
                        year = c(rep(2007, length(migrStartIdxs2007)), 
                                 rep(2008, length(migrStartIdxs2008))))
  zebDates$duration = zebDates$endDate - zebDates$startDate
  return(zebDates)
}

(zebSouthMigration <- getSouthMigDates())

# Calculate start and end of northern migration
# mean starting location is mean location between march 1st - april 1st, 
# here we look at population
(start2008 <- with(
  data.frame(subset(zebra_data_laea, 
                    timestamps(zebra_data_laea) >= "2008-03-01" & 
                      timestamps(zebra_data_laea) < "2008-04-01")), 
                   mean(complex(re = coords.x1, im =  coords.x2))))
(start2009 <- with(
  data.frame(subset(zebra_data_laea,
                    timestamps(zebra_data_laea) >= "2009-03-01" & 
                      timestamps(zebra_data_laea) < "2009-04-01")), 
                   mean(complex(re =  coords.x1, im =  coords.x2))))
startnorth <- c(start2008, start2009)

# mean ending location is mean location between june 1st - july 1
(end2008 <- with(
  data.frame(subset(zebra_data_laea, 
                    timestamps(zebra_data_laea) >= "2008-06-01" & 
                      timestamps(zebra_data_laea) < "2008-07-01")), 
                 mean(complex(re =  coords.x1, im = coords.x2))))
(end2009 <- with(
  data.frame(subset(zebra_data_laea, 
                    timestamps(zebra_data_laea) >= "2009-06-01" & 
                      timestamps(zebra_data_laea) < "2009-07-01")), 
                 mean(complex(re =  coords.x1, im =  coords.x2))))
endnorth <- c(end2008, end2009)

getZebNorthMigList = function(year = c("2008", "2009"))
{
  year = match.arg(year)
  if (year == "2008")
  {
    ids = unique(zebra_data_laea$individual.local.identifier[
      timestamps(zebra_data_laea) <= "2008-07-01"])
    zebmig = subset(zebra_data_laea, study.local.timestamp >= "2008-03-01" & 
                      study.local.timestamp < "2008-07-01" & 
                      individual.local.identifier %in% ids)
    start = start2008
    end = end2008
  } else if (year == "2009")
  {
    ids = unique(zebra_data_laea$individual.local.identifier[
      timestamps(zebra_data_laea) >= "2009-03-1"])
    zebmig = subset(zebra_data_laea, study.local.timestamp >= "2009-03-01" & 
                      study.local.timestamp < "2009-07-01" & 
                      individual.local.identifier %in% ids)
    start = start2009
    end = end2009
  }
  
  zebmig$z = complex(re = zebmig$coords.x1, im = zebmig$coords.x2)
  zebmig$distsStart = Mod(start - zebmig$z)
  zebmig$distsEnd = Mod(end - zebmig$z)
  zebmig$angle = Arg(zebmig$z - start) # angle from start to this point
  
  zebmig <- split(data.frame(zebmig), zebmig$individual.local.identifier)
  return(zebmig)
}

getNorthMigDates = function()
{
  zeb2009mig = getZebNorthMigList("2009")
  zeb2008mig = getZebNorthMigList("2008")
  
  # criteria: distance from start > 25 km
  (migrStartIdxs2009 = sapply(zeb2009mig, with, 
                              which(distsStart > 25000)[1]))
  (migrStartIdxs2008 = sapply(zeb2008mig, with, 
                              which(distsStart > 25000)[1]))
  
  (migrStartDates2009 = 
      lapply(1:length(zeb2009mig), 
             function(x) zeb2009mig[[x]]$study.local.timestamp[
               migrStartIdxs2009[x]]))
  (migrStartDates2008 = 
      lapply(1:length(zeb2008mig), 
             function(x) zeb2008mig[[x]]$study.local.timestamp[
               migrStartIdxs2008[x]]))
  
  # criteria: distance from end < 5 km or min distance 
  #(second zebra gets within 6.3 km)
  (migrEndIdxs2009 = 
      sapply(zeb2009mig, with, 
             min(which(distsEnd < 5000)[1], 
                 which(distsEnd == min(distsEnd)), na.rm = TRUE)))
  (migrEndIdxs2008 = 
      sapply(zeb2008mig, with, 
             min(which(distsEnd < 5000)[1], 
                 which(distsEnd == min(distsEnd)), na.rm = TRUE)))
  
  (migrEndDates2009 = 
      lapply(1:length(zeb2009mig), 
             function(x) 
               zeb2009mig[[x]]$study.local.timestamp[migrEndIdxs2009[x]]))
  (migrEndDates2008 = 
      lapply(1:length(zeb2008mig), 
             function(x) 
               zeb2008mig[[x]]$study.local.timestamp[migrEndIdxs2008[x]]))
  
  # need to specify timezone or else becomes current computer timezone, 
  # using UTC as defaut, though really CAS
  zebDates = data.frame(id = 
                          c(sapply(zeb2008mig, with, 
                                   individual.local.identifier[1]), 
                               sapply(zeb2009mig, with, 
                                      individual.local.identifier[1])), 
                        startIdx = c(migrStartIdxs2008, migrStartIdxs2009), 
                        startDate = as.POSIXct(
                          do.call(c, c(migrStartDates2008, 
                                       migrStartDates2009)), tz = "UTC"), 
                        endIdx = c(migrEndIdxs2008, migrEndIdxs2009), 
                        endDate = as.POSIXct(
                          do.call(c, c(migrEndDates2008, migrEndDates2009)), 
                          tz = "UTC"), 
                        year = c(rep(2008, length(migrStartIdxs2008)), 
                                 rep(2009, length(migrStartIdxs2009))))
  zebDates$duration = zebDates$endDate - zebDates$startDate
  return(zebDates)
}

(zebNorthMigration <- getNorthMigDates())

# Only take complete cases into account
zebSouthMigration <- zebSouthMigration[c(1:3,11:14),]
zebNorthMigration <- zebNorthMigration[c(1:3,11:14),]

# Define states
ids <- unique(zebra_data_laea$individual.local.identifier)
states <- list()
for(i in 1:length(ids)){
  zebra_sub <- zebra_data_laea[
    zebra_data_laea$individual.local.identifier == ids[i],]
  zebra_sub$states <- NA
  zebra_sub$states[timestamps(zebra_sub) <= 
                     zebSouthMigration$startDate[i]] <- 1
  zebra_sub$states[timestamps(zebra_sub) > zebSouthMigration$startDate[i] & 
                     timestamps(zebra_sub) <= 
                     zebSouthMigration$endDate[i]] <- 2
  if(!is.na(zebNorthMigration$startDate[i])){
    zebra_sub$states[timestamps(zebra_sub) > 
                       zebSouthMigration$endDate[i] & 
                       timestamps(zebra_sub) <= 
                       zebNorthMigration$startDate[i]] <- 3
  }else{  
    zebra_sub$states[timestamps(zebra_sub) > 
                       zebSouthMigration$endDate[i]] <- 3
  }
  zebra_sub$states[timestamps(zebra_sub) > 
                     zebNorthMigration$startDate[i] & 
                     timestamps(zebra_sub) <= 
                     zebNorthMigration$endDate[i]] <- 4
  zebra_sub$states[timestamps(zebra_sub) > 
                     zebNorthMigration$endDate[i]] <- 1
  states_ind <- list(zebra_sub$states)
  states <- append(states, states_ind); rm(states_ind, zebra_sub)
}

zebra_data_laea$state <- unlist(states)
# Save move object to file
saveRDS(zebra_data_laea, paste0(workdir, 
                                "Data/zebra_data_laea_", 
                                isos[i], ".rds"))

# Read study area shapefile 
studyarea_laea <- readOGR(dsn = "Results/", 
                          layer = paste0("studyarea_laea_", 
                                         isos[1]))

# Read cropped WDPA shapefile
wdpa_laea <- readOGR(dsn=paste0(workdir, "Results/wdpa_laea_", 
                                isos[1]), 
                     layer = "wdpa", verbose=FALSE)
wdpa_laea@data$id <- rownames(wdpa_laea@data)
df_wdpa <- fortify(wdpa_laea, region="id")
df_wdpa <- plyr::join(df_wdpa, wdpa_laea@data, by="id")

ggplot(data=studyarea_laea, aes(x=long, y=lat)) + 
  geom_polygon(aes(group=group), fill="snow2", colour=NA) + 
  geom_polygon(data=df_wdpa, aes(x=long, y=lat, group=group), 
               fill="forestgreen") + 
  geom_polygon(data=subset(df_wdpa, NAME=="Maun"), 
               aes(x=long, y=lat, group=group), 
               fill="forestgreen") + 
  geom_polygon(data=subset(df_wdpa, NAME=="Okavango Delta"), 
               aes(x=long, y=lat, group=group), fill="forestgreen") + 
  geom_polygon(data=subset(df_wdpa, NAME=="Moremi"), 
               aes(x=long, y=lat, group=group), fill="forestgreen") + 
  geom_polygon(data=studyarea_laea, aes(x=long,y=lat,group=group), 
               fill="transparent", colour="red") + 
  geom_point(data=data.frame(zebra_data_laea), 
                      aes(x=coords.x1, y=coords.x2, color=factor(state))) + 
  scale_colour_discrete(name="Segments", 
                      labels=c("Northern Range", "Southern Migration", 
                               "Southern Range", 
                               "Northern Migration")) + 
  labs(x="Eastings (m)", y="Northings (m)") + theme_bw() + 
  theme(legend.justification=c(1,1), legend.position=c(1,1),
        legend.background = element_rect(fill="transparent", colour=NA), 
        axis.text.y = rotatedAxisElementText(90,"y")) + 
  coord_equal(expand=FALSE)
ggsave(paste0("Figures/Zebra_Migration_", isos[1], "_Overview.pdf"), 
       width=6, height=5)

ggplot(data=studyarea_laea, aes(x=long, y=lat)) + 
  geom_polygon(aes(group=group), fill="snow2", colour=NA) + 
  geom_polygon(data=df_wdpa, aes(x=long, y=lat, group=group), 
               fill="forestgreen") + 
  geom_polygon(data=subset(df_wdpa, NAME=="Maun"), 
               aes(x=long, y=lat, group=group), 
               fill="forestgreen") + 
  geom_polygon(data=subset(df_wdpa, NAME=="Okavango Delta"), 
               aes(x=long, y=lat, group=group), fill="forestgreen") + 
  geom_polygon(data=subset(df_wdpa, NAME=="Moremi"), 
               aes(x=long, y=lat, group=group), fill="forestgreen") + 
  geom_polygon(data=studyarea_laea, aes(x=long,y=lat,group=group), 
               fill="transparent", colour="red") + 
  geom_point(data=data.frame(zebra_data_laea), 
                     aes(x=coords.x1, y=coords.x2, color=factor(state))) + 
  scale_colour_discrete(name="Segments", 
                      labels=c("Northern Range", "Southern Migration", 
                               "Southern Range", 
                               "Northern Migration")) + 
  facet_wrap(~ individual.local.identifier) + 
  labs(x="Eastings (m)", y="Northings (m)") + theme_bw() + 
  theme(axis.text.y = rotatedAxisElementText(90,"y")) + 
  coord_equal(expand=FALSE)
ggsave(paste0("Figures/Zebra_Migration_", isos[1], "_Individual.pdf"), 
       width=10, height=8)

# Calculate migration distances
dist_southmig <- sapply(
  distance(zebra_data_laea[zebra_data_laea$state == 2]), FUN=sum) 
# Southern Migration
dist_northmig <- sapply(
  distance(zebra_data_laea[zebra_data_laea$state == 4]), FUN=sum) 
# Northern Migration
dist_mig <- c(dist_southmig, NA, dist_northmig)

# Create table of migration dates for each individual
zebSouthMigration$dir <- "South"
zebNorthMigration$dir <- "North"
migration <- rbind(zebSouthMigration, zebNorthMigration)
migration$id <- levels(trackId(zebra_data_laea))[migration$id] 
migration$duration <- as.numeric(round(migration$duration))
migration <- migration[,c(1,8,3,5,7)]
migration[,3] <- as.Date(migration[,3])
migration[,4] <- as.Date(migration[,4])
migration$distance <- round(dist_mig)/1000
colnames(migration) <- c("ID", "Direction", "Start Date", 
                         "End Date", "Duration (days)", "Distance (km)")
migration[,3] <- as.factor(migration[,3])
migration[,4] <- as.factor(migration[,4])
migration[8,] <- c("Z3684", "North", NA, NA, NA, NA)
print(xtable(migration, 
caption=c("Summary of migration (direction, start and end date and duration) 
             for each individual.", 
          paste0("Migration summary for each individual, ", region[1])), 
label=paste0("table:summary_migrate_", isos[1])), type="latex", 
      file=paste0("Tables/Summary_migration_BWA.tex"), 
caption.placement="top", 
      include.rownames=FALSE, bookstab=TRUE, table.placement="H"); 
rm(migration)

## ---- goal ----
studyarea_laea <- readOGR(dsn=paste0(workdir, "Results/"), 
                          layer=paste0("studyarea_laea_", isos[1]))

# Get start and end point coordinates
start <- mean(c(startsouth, endnorth))
end <- mean(c(endsouth, startnorth))

# Create a start and a goal vector
goal <- SpatialPoints(
  coordinates(
    data.frame(x=c(Re(start),Re(end)), y=c(Im(start), Im(end)))), 
  proj4string = crs.laea)

# Create a distance to goal raster layer
if(!file.exists(paste0(workdir, "Results/goal_30m_", isos[1], ".tif"))){
  writeRaster(raster(resolution = c(30.7, 29.3), 
                     extent(studyarea_laea), vals=NA, 
                     crs=crs.laea), 
              paste0(workdir, "/Results/goal_30m_", isos[1], ".tif"), 
              overwrite=TRUE)
  goal <- distanceFromPoints(raster(paste0(workdir, "Results/goal_30m_", 
                                           isos[1], ".tif")), 
                             xy=goal, 
                             filename=paste0(workdir, "Results/goal_30m_", 
                                             isos[1], ".tif"), 
                             overwrite=TRUE)

  # Add goal data to zebra data
  zebra_data_laea <- readRDS(paste0(workdir, "Data/zebra_data_laea_", 
                                    isos[1], ".rds"))
  zebra_data_laea$goal <- raster::extract(goal, zebra_data_laea)
  saveRDS(zebra_data_laea, paste0(workdir, "Data/zebra_data_laea_", 
                                  isos[1], ".rds"))
}

# Plot goal layer
goal <- raster(paste0(workdir, "Results/goal_30m_", isos[1], ".tif"))
goal_df <- fortify(goal, maxpixels=maxpixels)
names(goal_df) <- c("x","y","goal")
ggplot() + geom_raster(data=goal_df, aes(x, y, fill = goal)) + 
  scale_fill_gradientn(name="Distance\nto Goal (m)", na.value="transparent", 
                       colours=colourtheme) + 
  geom_polygon(data=studyarea_laea, aes(x=long, y=lat, group=group), 
               fill="transparent", color="red") + 
  theme(legend.background = element_rect(fill="transparent", colour=NA), 
        axis.text.y = rotatedAxisElementText(90,"y")) + 
  labs(x="Eastings (m)", y="Northings (m)") + coord_equal(expand=FALSE)
ggsave(paste0("Figures/goal_30m_", isos[1], ".pdf"), height=4, width=6, 
       bg="transparent")

## ---- move_density ----

# Read study area shapefile 
studyarea_laea <- readOGR(dsn = "Results/", 
                          layer = paste0("studyarea_laea_", 
                                                           isos[2]))
studyarea <- readOGR(dsn = "Results/", 
                     layer = paste0("studyarea_", isos[2]))

# Read cropped WDPA shapefile
wdpa <- readOGR(dsn=paste0("Results/wdpa_laea_", isos[2]), layer = "wdpa", 
                verbose=FALSE)

# Read zebra movement data
zebra_data_laea <- readRDS(paste0("Data/zebra_data_laea_", isos[2], ".rds"))
zebra_data <- readRDS(paste0("Data/zebra_data_", isos[2], ".rds"))
if(class(zebra_data_laea) != "MoveStack"){
  zebra_data_laea <- moveStack(zebra_data_laea)
}
# Convert zebra data to spatial points
xy <- zebra_data_laea@coords
xy_wgs84 <- zebra_data@coords 

# Create empty raster
res_raster <- list(raster(resolution = c(30.7, 29.3), 
                          extent(studyarea_laea), crs=crs.laea), 
                   raster(resolution = c(475, 451), extent(studyarea_laea), 
                          crs=crs.laea), 
                   raster(resolution = c(237, 226), extent(studyarea_laea), 
                          crs=crs.laea))
avhrr_raster <- raster(resolution = c(0.08333333, 0.08333333), 
                       extent(studyarea), crs=crs.wgs84)

# Create density layer from empty raster
zebra_raster <- lapply(res_raster, FUN=function(x) 
  rasterize(xy, x, fun=function(x,...)length(x), na.rm=TRUE))
zebra_raster2 <- rasterize(xy_wgs84, 
                           avhrr_raster, fun=function(x,...)length(x), 
                           na.rm=TRUE)

# Convert to data frame
zebra_raster <- lapply(zebra_raster, FUN=function(x) 
  as.data.frame(rasterToPoints(x)))

# Get maximum density point
zebra_sp_max <- lapply(zebra_raster, FUN=function(x) x[
  x$count == max(x$count),])
zebra_avhrr_max <- zebra_raster2[
  zebra_raster2$count == max(zebra_raster2$count),]
coordinates(zebra_avhrr_max) <- ~x+y
projection(zebra_avhrr_max) <- crs.wgs84
zebra_avhrr_max <- spTransform(zebra_avhrr_max, crs.laea)
zebra_sp_max[[4]] <- c(x=zebra_avhrr_max@coords[1], 
                       y=zebra_avhrr_max@coords[2], 
                       zebra_avhrr_max@data)
zebra_sp_max <- do.call("rbind", zebra_sp_max)

## ---- corridor ----
for(i in 1:length(isos)){
  zebra_data <- readRDS(paste0(workdir, "Data/zebra_data_", isos[i], ".rds"))
  zebra_data <- split(zebra_data)
  stackcor <- lapply(zebra_data, corridor)
  corridors <- lapply(stackcor, FUN=function(x) c(NA, x@burstId))
  id <- lapply(stackcor, 
               FUN=function(x) c(as.character(
                 x@data$individual.local.identifier)))
  corridors_sp <- data.frame(
    do.call("rbind", mapply(x=stackcor, y=corridors, 
                            FUN=function(x,y) if(length(unique(y))==2){
                              x@coords[y == 2,]} else{x@coords[y == 1,]})))
  corridors_id <- data.frame(
    do.call("c", mapply(x=id, y=corridors, 
                        FUN=function(x,y) if(length(unique(y))==2){x[y == 2]} 
                        else{x[y == 1]})))
  corridors_sp <- cbind(corridors_sp, corridors_id)
  names(corridors_sp) <- c("coords.x1", "coords.x2", "ID")
  corridors_sp <- droplevels(corridors_sp)
  corridors_sp <- corridors_sp[!is.na(corridors_sp$ID),]
  coordinates(corridors_sp) <- ~coords.x1+coords.x2
  projection(corridors_sp) <- crs.wgs84
  # Load studyarea
  studyarea <- readOGR(dsn = paste0(workdir, "Results/"), 
                       layer = paste0("studyarea_", isos[i]), verbose=FALSE)
  # Read cropped WDPA shapefile
  wdpa <- readOGR(dsn=paste0(workdir, "Results/wdpa_sa_", isos[i]), 
                  layer = "wdpa", verbose=FALSE)
  
  # Number of corridors and number of corridors within PA
  nrow(corridors_sp)
  corridor_pa <- na.omit(over(corridors_sp, wdpa))
  nrow(corridor_pa)/nrow(corridors_sp)*100
  
  # Plot corridors
  ggplot() + geom_polygon(data=studyarea, aes(x=long, y=lat, group=group), 
                          fill="snow2", colour=NA) + 
    geom_polygon(data=wdpa, aes(x=long, y=lat, group=group), 
                 fill="darkseagreen2", colour="forestgreen") + 
    geom_polygon(data=studyarea, aes(x=long, y=lat, group=group), 
                 fill="transparent", colour="red") + 
    geom_point(data=corridors_sp, aes(x=coords.x1, y=coords.x2, colour=ID)) + 
    labs(x=expression(paste("Longitude (",degree,")")), 
         y=expression(paste("Latitude (",degree,")"))) + 
    coord_map()
  ggsave(paste0("Figures/Corridor_Points_", isos[i], ".pdf"), 
         width=8, height=6)
}

## ---- gimms ----

# Get GIMMS NDVI time series data
for(i in 1:length(isos)){
  
  if(file.exists(paste0(workdir, "Results/gimms_sa_", isos[i], ".grd"))){
    # Read NDVI sa data
    ndvi_sa <- stack(paste0(workdir, "Results/gimms_sa_", isos[i], ".grd"))
  }else{
    # Load gimms3g function
    source("Functions/gimms3g_function.R")
    
    # Read study area shapefile 
    studyarea <- readOGR(dsn = "Results/", 
                         layer = paste0("studyarea_", isos[i]), 
                         verbose=FALSE)
    # Get NDVI data for study area and save to file
    ndvi_sa <- gimms3g(extent=studyarea, 
                       path=paste0(filedir, "GIMMS"), 
                       snap="out")
    writeRaster(ndvi_sa, filename = 
                  paste0(workdir, "Results/gimms_sa_", isos[i], ".grd"), 
                bandorder = "BIL", overwrite=TRUE, 
                options = c("COMPRESS=NONE"), 
                dataType = "INT2U")
  }
  
  # Convert to crs.laea
  ndvi_sa <- projectRaster(ndvi_sa, crs=crs.laea)
  
  # Load WDPA laea and study area laea shapefile 
  studyarea_laea <- readOGR(dsn = "Results/", 
                            layer = paste0("studyarea_laea_", isos[i]), 
                            verbose=FALSE)
  wdpa_laea <- readOGR(dsn=paste0(workdir, "Results/wdpa_laea_", isos[i]), 
                       layer = "wdpa", verbose=FALSE)
  
  # Calculate linear regression for each pixel
  ndvi_sa_df_t <- t(data.frame(rasterToPoints(ndvi_sa)))
  date <- seq(from = as.Date(paste(1982, "01-15", sep="-")), 
              to = as.Date(paste(2013, "12-30", sep="-")), 
              length=24*(2013-1982+1))
  lin_reg_ndvi <- data.frame(x = ndvi_sa_df_t[1,], 
                             y = ndvi_sa_df_t[2,])
  for(j in 1:ncol(ndvi_sa_df_t)){
    # Run individual model
    m <- lm(ndvi_sa_df_t[c(-1,-2),j] ~ date)
    
    # Extract R2 value and save to dataframe
    lin_reg_ndvi$slope[j] <- m$coefficients[2]
    lin_reg_ndvi$adjrsq[j] <- summary(m)$adj.r.squared
    lin_reg_ndvi$pvalue[j] <- summary(m)$coefficients[2,4]
  }
  
  # Calculate % cover of significant positive slope
  lin_reg_ndvi_sig <- lin_reg_ndvi[lin_reg_ndvi$pvalue < 0.05,]
  (nrow(lin_reg_ndvi_sig[lin_reg_ndvi_sig$slope > 0,])*res(ndvi_sa)[1]*
      res(ndvi_sa)[2]*1e-06)/(nrow(lin_reg_ndvi)*res(ndvi_sa)[1]*
                                res(ndvi_sa)[2]*1e-06)*100
  
  # Create area with 5 km buffer around PAs
  wdpa_buf_5km <- gBuffer(wdpa_laea, width=5000, byid=TRUE)
  wdpa_buf_5km <- gIntersection(wdpa_buf_5km, studyarea_laea)
  
  # Create shapefile of unprotected areas
  unp <- gDifference(studyarea_laea, wdpa_buf_5km)
  
  # Just select buffer
  wdpa_buf_5km <- gDifference(wdpa_buf_5km, wdpa_laea)
  
  #Plot slope and adj. R2 for each pixel
  p1 <- ggplot(data = lin_reg_ndvi, aes(x = x, y = y), maxpixels=maxpixels) + 
    geom_raster(aes(fill=slope)) + 
    scale_fill_gradientn(name="Slope", colors = colourtheme, 
                         na.value="transparent") + 
    geom_polygon(data=wdpa_buf_5km, aes(x=long, y=lat, group=group), 
                 fill="transparent", colour="darkseagreen2") + 
    geom_polygon(data=wdpa_laea, aes(x=long, y=lat, group=group),
                 fill="transparent", colour="forestgreen") + 
    geom_polygon(data=studyarea_laea, aes(x=long, y=lat, group=group), 
                 fill="transparent", color="red") + 
    labs(x="Eastings (m)", y="Northings (m)") + 
    theme(legend.background = element_rect(fill="transparent", colour=NA), 
          axis.text.y = rotatedAxisElementText(90,"y")) + 
    guides(color = guide_legend(order=1), shape = guide_legend(order=2)) + 
    coord_equal(expand=FALSE)
  p2 <- ggplot(data = lin_reg_ndvi, aes(x = x, y = y), maxpixels=maxpixels) + 
    geom_raster(aes(fill=adjrsq)) + 
    scale_fill_gradientn(name="Adjusted R2", 
                         colors = colourtheme, na.value="transparent") + 
    geom_polygon(data=wdpa_buf_5km, aes(x=long, y=lat, group=group), 
                 fill="transparent", colour="darkseagreen2") + 
    geom_polygon(data=wdpa_laea, aes(x=long, y=lat, group=group), 
                 fill="transparent", colour="forestgreen") + 
    geom_polygon(data=studyarea_laea, aes(x=long, y=lat, group=group), 
                 fill="transparent", color="red") + 
    labs(x="Eastings (m)", y="Northings (m)") + 
    theme(legend.background = element_rect(fill="transparent", colour=NA), 
          axis.text.y = rotatedAxisElementText(90,"y")) + 
    guides(color = guide_legend(order=1), shape = guide_legend(order=2)) + 
    coord_equal(expand=FALSE)
  cowplot::plot_grid(p1, p2, align="h", nrow=1, labels="auto")
  if(i == 1){
    ggsave(paste0("Figures/Gimms_NDVI_", isos[i], ".pdf"), width=12, height=4)
  } else if(i == 2){
    ggsave(paste0("Figures/Gimms_NDVI_", isos[i], ".pdf"), width=8, height=4)
  }
  
  # Convert dataframe to raster
  coordinates(lin_reg_ndvi) <- ~x+y
  lin_reg_ndvi <- raster::rasterize(x=lin_reg_ndvi, y=ndvi_sa)[[2:4]]
  
  # Set pixels with not-significant relationship to NA
  is.na(lin_reg_ndvi$slope) <- lin_reg_ndvi$pvalue > 0.05
  is.na(lin_reg_ndvi$adjrsq) <- lin_reg_ndvi$pvalue > 0.05 
  
  # Mask data by protected, buffer and non-protected areas
  lin_reg_ndvi <- projectRaster(lin_reg_ndvi, crs=crs.laea)
  lin_reg_ndvi <- lin_reg_ndvi[[1:2]]
  lin_regs <- list(wdpa_laea, wdpa_buf_5km, unp)
  lin_regs <- lapply(lin_regs, FUN=function(x) 
    crop(mask(lin_reg_ndvi, x), x))
  
  # Boxplot of Slope and AdjR2 for protected areas, 
  # buffer zone and non-protected areas
  lin_regs <- lapply(lin_regs, FUN=function(x) 
    as.data.frame(rasterToPoints(x)))
  lin_regs_df <- do.call("rbind", lin_regs)
  status <- c("Protected", "5km Buffer", "Unprotected")
  lin_regs_df$status <- do.call("c", Map(x=lin_regs, 
                                         y=status, f=function(x,y) 
                                           rep(y, nrow(x))))
  lin_regs_df <- tidyr::gather(lin_regs_df, var, value, -c(x,y,status))
  lin_regs_df$var <- factor(lin_regs_df$var, 
                            labels=c("Slope", "Adjusted R2"))
  lin_regs_df$status <- factor(lin_regs_df$status)
  ggplot() + geom_boxplot(data=lin_regs_df, aes(x=status, 
                                                y=value, fill=status)) + 
    scale_x_discrete(limits=status) + 
    scale_fill_manual(values=c("darkseagreen2", "forestgreen", "snow2"), 
                      breaks=status) + 
    facet_wrap(~ var, scales="free_y") + labs(x="Status", y="") + 
    theme_bw() + theme(legend.position="none")
  ggsave(paste0("Figures/Gimms_Status_NDVI_", isos[i], ".pdf"), 
         width=10, height=5)

  # Mask data by protected and non-protected areas
  ndvi_pro <- crop(mask(ndvi_sa, wdpa_laea), wdpa_laea)
  # Create shapefile of unprotected areas
  unp <- gDifference(studyarea_laea, wdpa_laea)
  ndvi_unp <- crop(mask(ndvi_sa, unp), unp)
  
  # Calculate spatial average of NDVI and save in data frame
  ndvi_pro_sp <- as.data.frame(cbind(c(rep(1982:2013, each=24)), 
                                     c(rep(rep(1:12, each=2), times=32)) , 
                                     c(rep(seq(1, 16, by=15))), 
                                     cellStats(ndvi_pro, mean), 
                                     cellStats(ndvi_pro, sd)))
  names(ndvi_pro_sp) <- c("year", "month", "day","mn_ndvi", "sd_ndvi")
  ndvi_pro_sp$date <- as.Date(paste(c(rep(1982:2013, each=24)), 
                                    c(rep(rep(1:12, each=2), times=32)) , 
                                    c(rep(seq(1, 16, by=15)))), format="%Y%m%d")
  ndvi_unp_sp <- as.data.frame(cbind(c(rep(1982:2013, each=24)), 
                                     c(rep(rep(1:12, each=2), times=32)), 
                                     c(rep(seq(1, 16, by=15))), 
                                     cellStats(ndvi_unp, mean), 
                                     cellStats(ndvi_unp, sd)))
  names(ndvi_unp_sp) <- c("year", "month", "day","mn_ndvi", "sd_ndvi")
  ndvi_unp_sp$date <- as.Date(paste(c(rep(1982:2013, each=24)), 
                                    c(rep(rep(1:12, each=2), times=32)), 
                                    c(rep(seq(1, 16, by=15)))), format="%Y%m%d")
  
  #Combine protected and unprotected data to one dataframe
  ndvi_pro_sp$status <- "Protected"
  ndvi_unp_sp$status <- "Unprotected"
  ndvi_sp <- rbind(ndvi_pro_sp, ndvi_unp_sp)

  #Calculate 1 year running mean & sd of NDVI
  ndvi_pro_sp$mn_sym_1 <- stats::filter(ndvi_pro_sp$mn_ndvi, 
                                        rep(1/25, 25), sides=2)
  ndvi_pro_sp$sd_sym_1 <- stats::filter(ndvi_pro_sp$sd_ndvi, 
                                        rep(1/25, 25), sides=2)
  ndvi_unp_sp$mn_sym_1 <- stats::filter(ndvi_unp_sp$mn_ndvi, 
                                        rep(1/25, 25), sides=2)
  ndvi_unp_sp$sd_sym_1 <- stats::filter(ndvi_unp_sp$sd_ndvi, 
                                        rep(1/25, 25), sides=2)
  
  # Create plot of NDVI time series
  ggplot() + geom_ribbon(data=ndvi_pro_sp, aes(x=date, 
                                               ymin=mn_sym_1-sd_sym_1, 
                                               ymax=mn_sym_1+sd_sym_1), 
                         fill="forestgreen", alpha=0.5) + 
    geom_ribbon(data=ndvi_unp_sp, aes(x=date, ymin=mn_sym_1-sd_sym_1, 
                                      ymax=mn_sym_1+sd_sym_1), 
                fill="lightgrey", alpha=0.5) + 
    geom_line(data=ndvi_pro_sp, aes(x=date, y=mn_sym_1, colour="ProMean")) + 
    geom_line(data=ndvi_unp_sp, aes(x=date, y=mn_sym_1, colour="UnpMean")) + 
    scale_colour_manual(name="Status", breaks=c("ProMean", "UnpMean"), 
                        values=c("green3", "darkgrey"), 
                        labels=c("Protected", "Unprotected")) + 
    labs(x="Date", y="NDVI (1 year running Mean ± SD)")
  ggsave(paste0(workdir, "Figures/NDVI_1yr_mean_", isos[i], ".pdf"), 
         width=10, height=6, bg="transparent")
  
  #Calculate trend in NDVI for each protected area
  ndvi_wdpa <- list()
  for(j in 1:nrow(wdpa_laea)){
    ndvi_pro_sp <- data.frame(mn_ndvi=
                                apply((t(data.frame(
                                  raster::extract(ndvi_sa, wdpa_laea[j,])))), 
                                  1, mean, na.rm=TRUE), 
                              sd_ndvi= apply((t(
                                data.frame(raster::extract(ndvi_sa, 
                                                           wdpa_laea[j,])))), 
                                1, sd, na.rm=TRUE))
    ndvi_pro_sp$date <- as.Date(paste(c(rep(1982:2013, each=24)), 
                                      c(rep(rep(1:12, each=2), times=32)) , 
                                      c(rep(seq(1, 16, by=15)))), 
                                format="%Y%m%d")
    ndvi_pro_sp$status <- rep(wdpa_laea[j,]$NAME, length=nrow(ndvi_pro_sp))
    ndvi_wdpa[[j]] <- ndvi_pro_sp
  }
  ndvi_wdpa[[nrow(wdpa_laea)+1]] <- ndvi_unp_sp[,c(4,5,6,7)]
  ndvi_wdpa_df <- do.call("rbind", ndvi_wdpa)
  ggplot(data=ndvi_wdpa_df, aes(x=date, y=mn_ndvi)) + geom_line() + 
    labs(x="Date",y="Mean NDVI") + geom_smooth(method="lm", formula = y ~ x) + 
    facet_wrap(~ status) + 
    stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..rr.label.., 
                                                    sep = "~~~")), 
                 label.x.npc="right", label.y.npc="top", parse = TRUE) + 
    facet_wrap(~ status, ncol=2)
  ggsave(paste0("Figures/Mn_NDVI_WDPA_", isos[i], ".pdf"), width=7, height=9)
}

## ---- srtm30 ----

# Specify file names
srtm_v2 <- c("20160919155849_370254493", "20161205141802_228173302")

# Load & Process 30m SRTM DEM Data
for(i in 1:2){
  # Read study area shapefile 
  studyarea <- readOGR(dsn = "Results/", 
                       layer = paste0("studyarea_", 
                                      isos[i]), verbose=FALSE)
  
  if(file.exists(paste0(workdir, "Results/dem_", isos[i], ".tif"))){
    # Read file
    dem <- raster(paste0(workdir, "Results/dem_", isos[i], ".tif"))
  }else{
    # Read file
    dem <- raster(paste0(workdir, "Data/", srtm_v2[i], ".tif"))
    
    # Crop raster
    dem <- mask(crop(dem, studyarea), studyarea)
    dem <- projectRaster(dem, crs=crs.laea, 
                         filename=paste0(workdir, 
                                         "Results/dem_", isos[i], ".tif"), 
                         format="GTiff", overwrite=TRUE)
  }
  if(file.exists(paste0(workdir, "Results/slope_", isos[i], ".tif"))){
    slope <- raster(paste0(workdir, "Results/slope_", isos[i], ".tif"))
    aspect <- raster(paste0(workdir, "Results/aspect_", isos[i], ".tif"))
  }else{
    # Calculate slope from altitude data
    slope <- terrain(dem, opt='slope', unit='tangent', neighbors=8, 
                     filename=paste0(workdir, "Results/slope_", 
                                     isos[i], ".tif"), 
                     format="GTiff", overwrite=TRUE)
    aspect <- terrain(dem, opt='aspect', unit='degrees', neighbors=8, 
                      filename=paste0(workdir, "Results/aspect_", 
                                      isos[i], ".tif"), 
                      format="GTiff", overwrite=TRUE)
  }
  
  # Read zebra movement data
  zebra_data_laea <- readRDS(paste0(workdir, "Data/zebra_data_laea_", 
                                    isos[i], ".rds"))
  if(class(zebra_data_laea) != "MoveStack"){
    zebra_data_laea <- moveStack(zebra_data_laea)
  }
  
  # Read study area shapefile 
  studyarea_laea <- readOGR(dsn = "Results/", layer = 
                              paste0("studyarea_laea_", isos[i]))
  
  # Extract the raster information for each zebra location
  zebra_data_laea$dem <- raster::extract(dem, zebra_data_laea)
  
  # Plot dem
  dem_df <- fortify(dem, maxpixels=maxpixels); rm(dem)
  names(dem_df) <- c("x","y","dem")
  
  # Extract the raster information for each zebra location
  zebra_data_laea$slope <- raster::extract(slope, zebra_data_laea)
  zebra_data_laea$aspect <- raster::extract(aspect, zebra_data_laea)
  
  # Save zebra movement data
  saveRDS(zebra_data_laea, paste0(workdir, "Data/zebra_data_laea_", 
                                  isos[i], ".rds"))
  
  # Plot slope
  slope_df <- fortify(slope, maxpixels=maxpixels); rm(slope)
  names(slope_df) <- c("x","y","slope")
  
  # Plot aspect
  aspect_df <- fortify(aspect,maxpixels=maxpixels); rm(aspect)
  names(aspect_df) <- c("x","y","aspect")
  p1 <- ggplot() + geom_raster(data=dem_df, aes(x,y, fill = dem)) + 
    scale_fill_gradientn(name="DEM (m)", na.value="transparent", 
                         colours=terrain.colors(255)) + 
    geom_polygon(data=studyarea_laea, aes(x=long, y=lat, group=group), 
                 fill="transparent", color="red") + 
    theme(legend.background = element_rect(fill="transparent", colour=NA), 
          axis.text.y = rotatedAxisElementText(90,"y"), 
          legend.position="bottom") + 
    labs(x="Eastings (m)", y="Northings (m)") + coord_equal(expand=FALSE)
  rm(dem_df)
  p2 <- ggplot() + geom_raster(data=slope_df, aes(x,y, fill = slope)) + 
    scale_fill_gradientn(name="Slope (°)", na.value="transparent", 
                         colours=rev(terrain.colors(255))) + 
    geom_polygon(data=studyarea_laea, aes(x=long, y=lat, group=group), 
                 fill="transparent", color="red") + 
    theme(legend.background = element_rect(fill="transparent", colour=NA), 
          axis.text.y = rotatedAxisElementText(90,"y"), 
          legend.position="bottom") + 
    labs(x="Eastings (m)", y="Northings (m)") + coord_equal(expand=FALSE)
  rm(slope_df)
  p3 <- ggplot() + geom_raster(data=aspect_df, aes(x,y, fill = aspect)) + 
    scale_fill_gradientn(name="Aspect", na.value="transparent",
                         colours=rev(terrain.colors(255))) + 
    geom_polygon(data=studyarea_laea, aes(x=long, y=lat, group=group), 
                 fill="transparent", color="red") + 
    theme(legend.background = element_rect(fill="transparent", colour=NA), 
          axis.text.y = rotatedAxisElementText(90,"y"), 
          legend.position="bottom") + 
    labs(x="Eastings (m)", y="Northings (m)") + coord_equal(expand=FALSE)
  rm(aspect_df)
  cowplot::plot_grid(p1,p2,p3,ncol=3, labels="auto")
  if(i == 1){
    ggsave(paste0("Figures/terrain_", isos[i], ".pdf"), width=15, height=5)
  } else if(i == 2){
    ggsave(paste0("Figures/terrain_", isos[i], ".pdf"), width=10, height=6)
  }
}

## ---- terrain ----
for(i in 1:length(isos)){
  # Read zebra movement data
  zebra_data_laea <- readRDS(paste0(workdir, "Data/zebra_data_laea_", 
                                    isos[i], ".rds"))
  if(class(zebra_data_laea) != "MoveStack"){
    zebra_data_laea <- moveStack(zebra_data_laea)
  }
  
  # Plot zebras with elevation over time
  if(i == 1){
    ggplot(data=data.frame(zebra_data_laea), 
           aes(x = timestamps(zebra_data_laea),
               y = dem, colour=individual.local.identifier)) + 
      geom_line() + scale_colour_discrete(name="ID") + 
      theme_bw() + labs(x="Date", y="Altitude (m)")
  } else if(i == 2){
    df <- data.frame(zebra_data_laea[
      order(zebra_data_laea$individual.local.identifier, 
            timestamps(zebra_data_laea)),])
    df$date <- timestamps(zebra_data_laea)
    df$split <- 2007
    df$split[df$individual.local.identifier %in% 
               c("PZ6", "PZ8", "PZ10", "PZ14")] <- 2005
    ggplot(data=df, aes(x = date, y = dem, 
                        colour=individual.local.identifier)) + 
      geom_line() + scale_colour_discrete(name="ID") + theme_bw() + 
      labs(x="Date", y="Altitude (m)") + facet_wrap(~ split, scales="free_x")
  }
  ggsave(paste0("Figures/Zebra_dem_", isos[i], ".pdf"), width=10, height=5)
  
  # Create boxplot of Altitude per month zebras experience
  zebra_data_laea$month <- lubridate::month(timestamps(zebra_data_laea))
  ggplot(data=data.frame(zebra_data_laea), aes(factor(month), dem)) + 
    geom_boxplot() + 
    scale_x_discrete(breaks=c(1:12), 
                     labels=c("J","F","M","A","M","J",
                              "J","A","S","O","N","D")) + 
    labs(x="Month", y="Altitude (m)") + theme_bw()
  ggsave(paste0("Figures/dem_month_", isos[i], ".pdf"), width=9, height=7)
}

## ---- landsat_pathrow ----

# Read study area shapefile 
studyarea <- lapply(isos, FUN=function(x) readOGR(dsn = paste0(workdir, "/Results/"), layer = paste0("studyarea_", x), verbose=FALSE))

# Read zebra data
zebra_data <- lapply(isos, FUN=function(x) readRDS(paste0(workdir, "Data/zebra_data_", x, ".rds")))

# Retrieve path and row for required Landsat scenes
library(wrspathrow)
pr_sa <- lapply(studyarea, pathrow_num)

# Get polygons of required Landsat scenes
pr_sa_sp <- lapply(studyarea, pathrow_num, as_polys=TRUE)

for(i in 1:length(isos)){
  # Create plot of study area and the according Landsat scenes
  ggplot() + geom_polygon(data=pr_sa_sp[[i]], aes(x=long, y=lat, group=group), fill="gray", color="black") + 
    geom_text(data=pr_sa_sp[[i]]@data, aes(x = centroid_x, y = centroid_y, label = PATH), size=4) + 
    geom_text(data=pr_sa_sp[[i]]@data, aes(x = centroid_x, y = centroid_y, label = ROW), 
              size=4, nudge_x=0.07, nudge_y=-0.3) + 
    geom_polygon(data=studyarea[[i]], aes(x=long, y=lat, group=group), 
                 fill="transparent", colour="red", lwd=1.5) + 
    geom_path(data=data.frame(zebra_data[[i]]), aes(x=coords.x1, y=coords.x2, group=individual.local.identifier)) + 
    labs(x=expression(paste("Longitude (",degree,")")), y=expression(paste("Latitude (",degree,")"))) + 
    coord_cartesian(expand=FALSE)
  ggsave(paste0("Figures/Pathrow_", isos[i], ".pdf"), width=6, height=6)
}

## ---- landsat_times ----

# Specify path and rows
pathrows <- list(BWA=c("173073", "173074", "173075", "174073", "174074"), 
                 KEN=c("168059", "168060", "169059", "169060"))

# List Landsat directories
#landsat_dirs <- list.dirs(paste0(filedir, "EarthExplorer"))
#saveRDS(landsat_dirs, file="Results/Landsat_Dirs")
landsat_dirs <- readRDS("Results/Landsat_Dirs")
subset_dirs <- lapply(pathrows, FUN=function(x) grep(paste0(x, collapse="|"), 
                                                     landsat_dirs, value=TRUE))

# List doy and year of landsat images
landsat_alltimes <- lapply(subset_dirs, FUN=function(y)
  sort(unique(unlist(
    lapply(y, FUN=function(x)
      substring(unlist(strsplit(unlist(
        strsplit(x, split="/"))[8], split="-"))[1], 
        first=10, last=16))))))

# List all Landsat images that are available
landsat_scenes <- lapply(landsat_alltimes, FUN=function(x) {
  unlist(lapply(x, FUN= function(x) c(
    list.files(path = paste0(filedir, "EarthExplorer"), 
               pattern = paste0(x), full.names=TRUE))))})

# Get Scene Names
scene_names <- lapply(subset_dirs, FUN=function(y) 
  unlist(lapply(y, FUN=function(x) 
    substring(unlist(strsplit(unlist(
      strsplit(x, split="/"))[8], split="-"))[1], 
      first=1, last=16))))

# Check for Duplicates
print(scene_names[[1]][duplicated(scene_names[[1]])])
print(scene_names[[2]][duplicated(scene_names[[2]])])

# Get Dates
landsat_alldates <- lapply(landsat_alltimes, function(x) 
  as.Date(as.numeric(substring(x, first=5)) - 1, 
          origin = paste0(substring(x, first=1, last=4), "-01-01")))

## ---- landsat ----
for(i in 1:2){
  for(k in 1:length(landsat_scenes[[i]])){
    # Get name of scene and check if files exist
    scene_name <- strsplit(strsplit(
      list.files(paste0(landsat_scenes[[i]][k]), pattern=".xml", 
                 full.names=TRUE), split="/")[[1]][8], split="-")[[1]][1]
    sensor <- substring(scene_name, first=1, last=3)
    if(file.exists(paste0(workdir, "Results/", sensor, "/", scene_name ,".tif"))){
      if(!file.exists(paste0(workdir, "Figures/", sensor, "/", scene_name, ".pdf"))){
        ls_pr_cdr <- stack(paste0("/home/mabi/GitHub/ZebraMove/Results/", sensor, "/", scene_name,".tif"))
        ggRGB(ls_pr_cdr, r = 1, g = 2, b = 3, stretch = "lin")
        ggsave(paste0("Figures/", sensor, "/", scene_name, ".pdf"), width=9, height=7)
      }
    }else{
      # Load list of Landsat images for one time frame
      xml_meta <- readMeta(list.files(paste0(landsat_scenes[[i]][k]), pattern=".xml", full.names=TRUE))
      
      # Load surface reflectance
      ls_pr_cdr <- stackMeta(xml_meta, quantity = c("sre"))
      
      # Read study area shapefile 
      studyarea <- readOGR(dsn = "Results/", layer = paste0("studyarea_", isos[i]), verbose=FALSE)
      
      # Create studyarea shapefile with same projection as landsat image
      studyarea_lsat <- spTransform(studyarea, crs(ls_pr_cdr)); rm(studyarea)
      
      # Get scale Factor
      scaleF <- getMeta(ls_pr_cdr, xml_meta, what = "SCALE_FACTOR")
      
      # Crop Landsat image by extent of studyarea
      ls_pr_cdr <- crop(ls_pr_cdr, studyarea_lsat); rm(studyarea_lsat)
      
      # Multiple data by scaling factor
      ls_pr_cdr <- ls_pr_cdr * scaleF; rm(scaleF)
      
      # Project and mask dem according to landsat scene
      dem <- raster(paste0(workdir, "Results/dem_", isos[i], ".tif"))
      beginCluster(n)
      dem <- projectRaster(dem, crs=projection(ls_pr_cdr), 
                           filename="dem_lsat.tif", format="GTiff", 
                           overwrite=TRUE)
      endCluster()
      beginCluster(n)
      dem <- resample(dem, ls_pr_cdr[[1]])
      endCluster() 
      dem <- crop(dem, ls_pr_cdr[[1]])
      
      # Apply topographic illumination correction
      beginCluster(n)
      ls_pr_cdr <- topCor(ls_pr_cdr, dem = dem, metaData = xml_meta, method = "C")
      endCluster(); rm(dem)
      
      # Obtain cloud product and crop it by extent of the study area
      ls_pr_qa <- stackMeta(xml_meta, category = "qa")
      ls_pr_qa <- crop(ls_pr_qa, ls_pr_cdr); rm(xml_meta)
      
      # Mask clouds and shadows
      ls_pr_cdr <- mask(ls_pr_cdr, ls_pr_qa[["QA_cloud"]], 
                        maskvalue=255)
      ls_pr_cdr <- mask(ls_pr_cdr, ls_pr_qa[["QA_cloud_shadow"]], 
                        maskvalue=255)
      ls_pr_cdr <- mask(ls_pr_cdr, ls_pr_qa[["QA_adjacent_cloud"]], 
                        maskvalue=255); rm(ls_pr_qa)
      
      # Project raster files to country laea
      beginCluster(n)
      ls_pr_cdr <- projectRaster(ls_pr_cdr, crs=crs.laea)
      endCluster()
      
      # Save modified landsat images to file
      scene_name <- strsplit(strsplit(list.files(paste0(landsat_scenes[[i]][k]), 
                                                 pattern=".xml", full.names=TRUE), 
                                      split="/")[[1]][8], split="-")[[1]][1]
      if(!file.exists(paste0(workdir, "Figures/", sensor, "/", scene_name, ".pdf"))){
        ggRGB(ls_pr_cdr, r = 1, g = 2, b = 3, stretch = "lin")
        ggsave(paste0("Figures/", sensor, "/", scene_name, ".pdf"), width=9, height=7)
      }
      
      ## Write raster
      beginCluster(n)
      ls_pr_cdr <- projectRaster(ls_pr_cdr, crs=crs.laea, 
                                 filename=paste0("/home/mabi/GitHub/ZebraMove/Results/", 
                                                 sensor, "/", scene_name,".tif"), 
                                 format = "GTiff", overwrite = TRUE)
      endCluster(); rm(ls_pr_cdr)
      removeTmpFiles(h=0.2)
    }
  }
}

## ---- landsat_mosaick ----

# List doy and year of landsat images
landsat_times <- list(zb1=c(2007282, 2007298, 2007314, 2007346, 2007362, 
                            2008013, 2008029, 2008037, 2008053, 2008069, 
                            2008085, 2008101, 2008117, 2008157, 2008173, 
                            2008189, 2008205, 2008221, 2008237, 2008253, 
                            2008269, 2008285, 2008301, 2008333, 2009015, 
                            2009031, 2009047, 2009079, 2009095, 2009111, 
                            2009127, 2009143), 
                      zb2=c(2005145,2005161,2005193, 2007151, 2007183, 2007199, 
                            2007215, 2007231, 2007247, 2007263, 2007279, 2007295, 
                            2007311, 2007327, 2007343, 2007359))
landsat_dates <- lapply(landsat_times, function(x) 
  as.Date(as.numeric(substring(x, first=5)) - 1, 
          origin = paste0(substring(x, first=1, last=4), "-01-01")))

# List Landsat images
for(i in 1:2){
  # Read study area shapefile 
  studyarea_laea <- readOGR(dsn = "Results/", layer = paste0("studyarea_laea_", isos[i]), verbose=FALSE)
  # Loop for running through each date 
  for(k in 1:length(unique(landsat_times[[i]]))){
    sensor <- c("LT5", "LE7")
    for(l in 1:length(sensor)){
      # Get name of scene and load file
      (landsat_images <- c(list.files(path = paste0(workdir, "Results/", sensor[l]), 
                                      pattern = paste0(strftime(landsat_dates[[i]][k], format = "%Y%j")), 
                                      full.names=TRUE), 
                           list.files(path = paste0(workdir, "Results/", sensor[l]), 
                                      pattern = paste0(strftime(landsat_dates[[i]][k]+7, format = "%Y%j")), 
                                      full.names=TRUE), 
                           list.files(path = paste0(workdir, "Results/", sensor[l]), 
                                      pattern = paste0(strftime(landsat_dates[[i]][k]+8, format = "%Y%j")), 
                                      full.names=TRUE), 
                           list.files(path = paste0(workdir, "Results/", sensor[l]), 
                                      pattern = paste0(strftime(landsat_dates[[i]][k]-1, format = "%Y%j")), 
                                      full.names=TRUE)))
      
      if(length(landsat_images) > 0){
        if(!file.exists(paste0("/home/mabi/GitHub/ZebraMove/Results/LS_Merged/", 
                               sensor[l], "_", isos[i], "_", landsat_times[[i]][k], ".grd"))){
          # read the individual rasters into a list and resample them
          merge_list <- lapply(landsat_images, stack)
          if(length(landsat_images) > 1){
            ls_stack <- raster(resolution = c(30.7, 29.3), extent(studyarea_laea), crs=crs.laea)
            merge_list <- lapply(merge_list, resample, ls_stack); rm(ls_stack)
            
            # Mosaic all the rasters together
            merge_list$tolerance <- 0.3
            merge_list$fun <- "mean"
            beginCluster(n)
            merged <- do.call(mosaic, merge_list)
            endCluster()
            rm(merge_list)
          }else if(length(landsat_images) == 1){
            ls_stack <- raster(resolution = c(30.7, 29.3), extent(studyarea_laea), crs=crs.laea)
            merge_list <- lapply(merge_list, resample, ls_stack); rm(ls_stack)
            merged <- merge_list[[1]]; rm(merge_list)
          }
          # Save modified landsat images to file
          writeRaster(merged, filename=paste0("/home/mabi/GitHub/ZebraMove/Results/LS_Merged/", 
                                              sensor[l], "_", isos[i], "_", landsat_times[[i]][k], ".grd"), 
                      bandorder = "BIL", options = c("COMPRESS=NONE"), dataType = "INT2U", overwrite=TRUE)
          rm(merged)
          removeTmpFiles(h = 0.01)
        }
        if(!file.exists(paste0("Figures/LS_Merged/", sensor[l], "_", 
                               isos[i], "_", landsat_times[[i]][k], ".pdf"))){
          merged <- stack(paste0("/home/mabi/GitHub/ZebraMove/Results/LS_Merged/", 
                                 sensor[l], "_", isos[i], "_", landsat_times[[i]][k], ".grd"))
        ggRGB(merged, r = 1, g = 2, b = 3, stretch = "lin") + 
          labs(x=expression(paste("Longitude (",degree,")")), 
               y=expression(paste("Latitude (",degree,")"))) + 
          geom_polygon(data=studyarea_laea, aes(x=long, y=lat, group=group), 
                       fill="transparent", color="red")
        coord_cartesian(expand=FALSE)
        ggsave(paste0("Figures/LS_Merged/", sensor[l], "_", isos[i], "_", 
                      landsat_times[[i]][k], ".pdf"), width=9, height=7)
        }
      }
    }
  }
}

## ---- ls_timeseries ----
bands <- c("b01", "b02", "b03", "b04", "b05", "b07")
for(i in 1:2){
  if(unique(sapply(bands, FUN=function(x) !file.exists(paste0(workdir, "Results/LS_", x, "_", isos[i], ".grd")))) == TRUE){
    # Read study area shapefile 
    studyarea_laea <- readOGR(dsn = "Results/", layer = paste0("studyarea_laea_", isos[i]), verbose=FALSE)
    
    # Get file names
    ls_list <- list.files(path = paste0(workdir, "Results/LS_Merged"), pattern = isos[i], full.names=TRUE)
    ls_list <- ls_list[c(grep(x=ls_list, pattern=".grd"))]
    
    # Get dates
    ls_dates <- lapply(ls_list, FUN= function(x) unlist(strsplit(unlist(strsplit(x[1], split="_"))[4], split=".", fixed=TRUE))[1])
    ls_dates <- lapply(ls_dates, function(x) as.Date(as.numeric(substring(x, first=5)) - 1, origin = paste0(substring(x, first=1, last=4), "-01-01")))
    ls_dates <- do.call("c", ls_dates)
    
    # Load raster stacks
    ls_list <- lapply(ls_list, stack)
    
    for(j in 1:length(bands)){
      if(!file.exists(paste0(workdir, "Results/LS_", bands[j], "_", isos[i], ".grd"))){
        # Set names
        ls_band <- stack(lapply(ls_list, function(x) x[[j]]))
        ls_band <- setZ(ls_band, z=ls_dates, name="date")
        names(ls_band) <- getZ(ls_band)
        
        # Save Bands to file
        writeRaster(ls_band, filename = paste0(workdir, "Results/LS_", bands[j], "_", isos[i], ".grd"), bandorder = "BIL", options = c("COMPRESS=NONE"), dataType = "INT2U", overwrite=TRUE); rm(ls_band)
        removeTmpFiles(h=0.01)
      }
    };rm(ls_dates, ls_list)
  }
}

## ---- modis_tiles ----
# install.packages("MODIS", repos = "http://R-Forge.R-project.org")
library(MODIS)
for(i in 1:2){
  # Read study area shapefile 
  studyarea <- readOGR(dsn = "Results/", layer = paste0("studyarea_", isos[i]), verbose=FALSE)
  
  # Read zebra data
  zebra_data <- readRDS(paste0(workdir, "Data/zebra_data_", isos[i], ".rds"))
  
  # Specify tiles
  tileH <- getTile(extent=studyarea)$tileH
  tileV <- getTile(extent=studyarea)$tileV
  
  # Get extent of tiles
  modis_tiles <- genTile(tileSize=10)
  sa_tiles <- subset(modis_tiles, iv %in% tileV & ih == tileH)
  # sp_tiles <- apply(data.frame(sa_tiles[3], sa_tiles[5], sa_tiles[6], sa_tiles[4]), MARGIN=1, FUN=extent)
  sa_tiles$centroid_x <- rowMeans(sa_tiles[,c(3,5)])
  sa_tiles$centroid_y <- rowMeans(sa_tiles[,c(4,6)])
  
  # Plot MODIS tiles and extent study area
  ggplot() + geom_rect(data=sa_tiles, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
                       fill="transparent", colour="black") + 
    geom_text(data=sa_tiles, aes(x = centroid_x, y = centroid_y, label = iv), size=4) + 
    geom_text(data=sa_tiles, aes(x = centroid_x, y = centroid_y, label = ih), 
              size=4, nudge_x=0.05, nudge_y=-1) + 
    labs(x=expression(paste("Longitude (",degree,")")), y=expression(paste("Latitude (",degree,")"))) + 
    geom_polygon(data=studyarea, aes(x=long, y=lat, group=group), fill="transparent", color="red") + 
    geom_path(data=data.frame(zebra_data), aes(x=coords.x1, y=coords.x2, group=individual.local.identifier)) + 
    coord_map()
  ggsave(paste0("Figures/tileH_V_", isos[i], ".pdf"), width=7, height=7)
}

## ---- mlc ----
MODISoptions(localArcPath = paste0(filedir, "MODIS"), 
             outDirPath = paste0(filedir, "MODIS/Processed"))
product <- "MCD12Q1" # LandCoverType, 500m, Yearly

for(i in 1:2){
  if(!file.exists(paste0(workdir,"Results/lcov_", product, "_", isos[i], "_brick.grd"))){
    # Read study area shapefile 
    studyarea <- readOGR(dsn = "Results/", layer = paste0("studyarea_", isos[i]))
    studyarea_sin <- spTransform(studyarea, 
                                 CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 
                                     +a=6371007.181 +b=6371007.181 +units=m +no_defs"))
    
    # Specify tiles and time period
    tileH <- getTile(extent=studyarea)$tileH
    tileV <- getTile(extent=studyarea)$tileV
    
    # Create empty list
    merge_list <- list()
    
    # Loop through all tiles
    for(j in 1:length(tileV)){
      # Download data
      modis.hdf <- getHdf(product = product, tileH = tileH, tileV = tileV[j], checkIntegrity=TRUE)
      lcov_all <- stack()
      for(k in 1:length(modis.hdf[[1]])){
        ifl <- strsplit(basename(modis.hdf[[1]][k]), ".hdf")
        fn <- paste0("lcov_", ifl, ".tif")
        if(is.na(ifl) | file.exists(paste0(filedir, "MODIS/Processed/", fn))){
          print("file exists or is not available on the server")
        }else{
          # specify layer you want to extract (NDVI)
          sds <- get_subdatasets(modis.hdf[[1]][k]) 
          #Convert to Tiff
          tmp <- rasterTmpFile()
          extension(tmp) <- "tif"
          gdal_translate(sds[1], dst_dataset = tmp) # Landcover one out of 5
          lcov <- crop(x = raster(tmp), y = studyarea_sin, 
                       filename=paste0(filedir, "MODIS/Processed/", fn), 
                       format="GTiff", overwrite=TRUE); rm(lcov)
        }
        lcov_all <- addLayer(lcov_all, raster(paste0(filedir, "MODIS/Processed/", fn)))
      }
      # Create list of scenes to merge
      merge_list  <- append(merge_list, lcov_all);rm(lcov_all)
    }
    if(length(merge_list)>1){
      lcov_all <- stack()
      for(l in 1:nlayers(merge_list[[1]])){
        lcov_ind <- merge(merge_list[[1]][[l]], merge_list[[2]][[l]], tolerance=0.05)
        beginCluster(n)
        lcov_ind <- projectRaster(lcov_ind, method = "ngb", crs = crs.laea)
        endCluster()  
        lcov_all <- stack(lcov_all, lcov_ind); rm(lcov_ind)
      }
    } else{
      lcov_all <- merge_list[[1]]; rm(merge_list)
      beginCluster(n)
      lcov_all <- projectRaster(lcov_all, method = "ngb", crs = crs.laea)
      endCluster()
    }
    # Give name to layers  
    names(lcov_all) <- extractDate(modis.hdf[[1]], asDate=TRUE)[[1]]; rm(modis.hdf)
    
    # Save as grid, as this will keep the layer names
    writeRaster(lcov_all, filename = paste0(workdir,"Results/lcov_", product, "_", 
                                            isos[i], "_brick.grd"), 
                bandorder = "BIL", overwrite=TRUE, options = c("COMPRESS=NONE"), 
                dataType = "INT2U"); rm(lcov_all)
  }
}

## ---- mlc_plot ----
# Land cover classes
classes <- data.frame(ID = c(0:16, 254, 255), 
                      Class = factor(c("Water", "Evergreen Needleleaf forest", 
                                "Evergreen Broadleaf forest", 
                                "Deciduous Needleleaf forest", 
                                "Deciduous Broadleaf forest", 
                                "Mixed forest", "Closed shrublands", 
                                "Open shrublands", "Woody savannas", 
                                "Savannas", "Grasslands", 
                                "Permanent wetlands", "Croplands", 
                                "Urban and built-up", "Cropland/Natural vegetation mosaic", 
                                "Snow and ice", "Barren or sparsely vegetated", 
                                "Unclassified", "Fill Value"), ordered=TRUE), 
                      colours=factor(c("blue", NA, "darkgreen", NA, 
                                "green2", "green3", "brown", 
                                "coral1", "lightseagreen", 
                                "seagreen2", "green", "lightblue", 
                                "yellow", "grey", "orange", "white", 
                                "brown2", "transparent", "transparent"), ordered=TRUE))

# Create table of Landcover classes
print(xtable(classes[c(1:17),c(1:2)], 
             caption=c("Land cover classes of the IGBP scheme.", "Land cover classes of the IGBP scheme"), 
             label="table:mlc_classes"), include.rownames=FALSE, type="latex", 
      file = paste0("Tables/MLC_classes.tex"), 
      caption.placement="top", booktabs=TRUE, table.placement="H")

# Re-set colnames
colnames(classes) <- c("id", "class", "colours")

# Plot Landcover data
for(i in 1:2){
  # Read study area shapefile 
  studyarea_laea <- readOGR(dsn = paste0(workdir, "Results/"), 
                            layer = paste0("studyarea_laea_", isos[i]))
  
  # Read Landcover data
  lcov <- stack(paste0(workdir,"Results/lcov_MCD12Q1_", isos[i], "_brick.grd"))
  lcov <- crop(lcov, studyarea_laea)
  
  # Calculate area of each land cover class
  lcov_freq <- Reduce(function(x,y) merge(x, y, by="value", all=TRUE), 
                      freq(lcov, useNA = "no"))
  years <- 2001:2013
  names(lcov_freq) <- c("id", years)
  lcov_area <- lcov_freq
  lcov_perc <- lcov_freq
  for(j in 1:length(years)){
    lcov_area[,j+1] <- lcov_freq[,j+1]*res(lcov)[1]*res(lcov)[2]*1e-06
    lcov_perc[,j+1] <- sapply(lcov_area[,j+1], FUN=function(x) 
      round(x/sum(lcov_area[,j+1], na.rm=TRUE)*100, 2))
  }
  lcov_area <- plyr::join(lcov_area, classes, by = "id", type="left", match="first")
  lcov_perc <- plyr::join(lcov_perc, classes, by = "id", type="left", match="first")
  
  # Plot time-series steamgraph of land cover area
  lcov_ts <- tidyr::gather(lcov_area[,1:15], year, area, -c(id, class))
  lcov_ts$year <- as.numeric(lcov_ts$year)
  lcov_ts <- lcov_ts[complete.cases(lcov_ts),]
  ggplot(data = lcov_ts, aes(x = year, y = area, group = class, fill = factor(id))) + 
    geom_area() + labs(list(x = "Year", y = "Area (km2)")) + 
    scale_x_continuous(breaks = seq(min(lcov_ts$year), max(lcov_ts$year), by = 1)) +
    scale_fill_manual(values=as.character(classes$colours[classes$id %in% levels(factor(lcov_ts$id))]), 
                      labels=as.character(classes$class[classes$id %in% levels(factor(lcov_ts$id))]), 
                      name="Land Cover") + theme_bw() + theme(legend.position="bottom")
  ggsave(paste0("Figures/MLC_Timeseries_", isos[i], ".pdf"), width=9, height=6)
  
  # Create table
  years_sub <- list(BWA=c("2007", "2008", "2009"), KEN=c("2005", "2007", "2008"))
  lcov_sa <- cbind(lcov_perc$class, lcov_perc[,c(paste0(years_sub[[i]]))])
  colnames(lcov_sa)[1] <- c("Landcover Class")
  print(xtable(lcov_sa, caption=c(paste0("Percentage of each land cover class of the 
                                       studyarea per year over the course of the study."), 
                                  paste0("Percentage of each LC class per year, ", region[i])), 
               label=paste0("table:mlc_percentage_", isos[i]), digits=2), 
        type="latex", file = paste0("Tables/MLC_percentage_", isos[i], ".tex"), 
        caption.placement="top", include.rownames=FALSE, 
        booktabs=TRUE, table.placement="H")
  
  # Read study area shapefile 
  wdpa_laea <- readOGR(dsn = paste0(workdir, "Results/wdpa_laea_", isos[i]), layer = "wdpa")
  
  # Create table of mean percentage of land cover per PA
  lcov_pas <- list()
  mn_lcov <- calc(lcov, mean)
  for(j in 1:nrow(wdpa_laea)){
    lcov_pa <- crop(mn_lcov, wdpa_laea[j,])
    lcov_pa_freq <- data.frame(freq(lcov_pa, useNA = "no"))
    lcov_pa_freq[,2] <- lcov_pa_freq[,2]* res(lcov_pa)[1]* res(lcov_pa)[2] * 1e-06
    names(lcov_pa_freq) <- c("id", "area")
    lcov_pa_freq$percentage <- 
      sapply(lcov_pa_freq[,2], FUN=function(x) 
        round(x/sum(lcov_pa_freq[,2])*100, 2))
    lcov_pa_freq <- 
      plyr::join(lcov_pa_freq, classes, by = "id", type="left", match="first")
    lcov_pa_freq$name <- 
      rep(wdpa_laea[j,]$NAME, length=nrow(lcov_pa_freq))
    lcov_pas[[j]] <- lcov_pa_freq 
  }
  lcov_unp <- mask(mn_lcov, wdpa_laea, inverse=TRUE)
  lcov_unp_freq <- data.frame(freq(lcov_unp, useNA="no"))
  lcov_unp_freq[,2] <- lcov_unp_freq[,2]* res(lcov_unp)[1]* res(lcov_unp)[2] * 1e-06
  names(lcov_unp_freq) <- c("id", "area")
  lcov_unp_freq$percentage <- 
    sapply(lcov_unp_freq[,2], FUN=function(x) 
      round(x/sum(lcov_unp_freq[,2])*100, 2))
  lcov_unp_freq <- plyr::join(lcov_unp_freq, classes, by = "id", type="left", match="first")
  lcov_unp_freq$name <- rep("Non-protected", length=nrow(lcov_unp_freq))
  lcov_pas[[nrow(wdpa_laea)+1]] <- lcov_unp_freq
  lcov_pas_df <- do.call("rbind", lcov_pas)
  lcov_pas_df <- lcov_pas_df[c("class", "name", "percentage")]
  lcov_pas_df <- tidyr::spread(lcov_pas_df, key="name", value="percentage")
  colnames(lcov_pas_df)[1] <- ""
  print(xtable(lcov_pas_df, caption=c("Percentage of land cover 
                                           for protected and non-protected 
                                           areas.", paste0("Percentage of LC 
                                           for protected and non-protected 
                                           areas, ", region[i])), 
               label=paste0("table:mlc_percentage_wdpa_", isos[i]), digits=1), 
        type="latex", file = paste0("Tables/MLC_Percentage_WDPA_", isos[i], ".tex"), 
        include.rownames=FALSE, floating=TRUE, table.placement="H", 
        caption.placement="top", booktabs=TRUE, 
        rotate.colnames=TRUE) #floating.environment = "sidewaystable",
  
  # Fortify Landcover
  lcov_df <- as.data.frame(rasterToPoints(lcov)); rm(lcov)
  years <- 2001:2013
  colnames(lcov_df) <- c("x", "y", years)
  lcov_df <- tidyr::gather(lcov_df, "year", "lcov",3:ncol(lcov_df))
  colnames(lcov_df) <- c("x", "y", "year", "id")
  
  # Merge lcov with classes
  lcov_df <- plyr::join(lcov_df, classes, by = "id", type="left", match="first")
  
  # Create land cover map
  ggplot() + geom_raster(data=lcov_df, 
                         aes(x,y,fill=factor(id))) + 
    labs(x="Eastings (m)", y="Northings (m)") + 
    scale_fill_manual(values=as.character(classes$colours[classes$id %in% levels(factor(lcov_df$id))]), 
                      labels=as.character(classes$class[classes$id %in% levels(factor(lcov_df$id))]), 
                      name="Land Cover") + 
    geom_polygon(data=studyarea_laea, aes(long,lat, group=group), 
                 fill="transparent", colour="red") + 
    theme(axis.text.y = rotatedAxisElementText(90,"y")) + coord_equal(expand=FALSE) + 
    facet_wrap(~ year)
  if(i == 1){
    ggsave(paste0("Figures/Landcover_", isos[i], ".pdf"), width=10, height=7)
  } else if(i == 2){
    ggsave(paste0("Figures/Landcover_", isos[i], ".pdf"), width=8, height=10)
  }
}

## ---- modis_vi ----
MODISoptions(localArcPath = paste0(filedir, "MODIS"), 
             outDirPath = paste0(filedir, "MODIS/Processed"))

# Function to discard NDVI time series pixels that do not fulfil a certain minimum criterion for quality
f_cleanvi <- function(vi, rel){
  i <- (rel <=1)
  res <- matrix(NA, length(i), 1)
  if(sum(i, na.rm = TRUE) > 0){
    i <- which(i)
    res[i] <- vi[i]
  }
  res
}

for(i in 1:length(isos)){
  # Read study area shapefile 
  studyarea <- readOGR(dsn = "Results/", layer = paste0("studyarea_", isos[i]))
  
  if(!file.exists(paste0(workdir, "Results/M.D13_NDVI_", isos[i], ".grd"))){
    # Transform studyarea shapefile to sinusoidal projection
    studyarea_sin <- spTransform(studyarea, CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 
                                                +b=6371007.181 +units=m +no_defs")); rm(studyarea)
    
    # Specify start and end period
    begin <- c("2007.10.01", "2005.06.01")
    end <- c("2009.05.02", "2008.07.01")
    
    # Download data
    if(!dir.exists(paste0(options("MODIS_outDirPath"), isos[i], "_VI"))){
      modis.hdf <- getHdf(product = "M.D13Q1", begin = begin[i], end = end[i], 
                          extent=studyarea_sin, checkIntegrity=TRUE)
      
      # Extract GeoTiff for each time stamp and SDS Layer specified
      runGdal(job=paste0(isos[i], "_VI"), product = "M.D13Q1", 
              begin = begin[i], end = end[i], extent=studyarea_sin, 
              SDSstring="111000000000", checkIntegrity=TRUE)
    }
    
    # Specify path of files
    path <- paste0(options("MODIS_outDirPath"), isos[i], "_VI")
    
    # Load band stacks and apply conversion factor
    vi <- c("NDVI", "EVI")
    preStack_list <- lapply(vi, FUN=function(x) 
      preStack(path=path, pattern=paste0(x, ".tif")))
    stack_list <- lapply(preStack_list, stack)
    beginCluster(n)
    stack_list <- lapply(stack_list, FUN=function(x){
      clusterR(x, fun=function(y) y/10000)})
    endCluster()
    crop_list <- lapply(stack_list, FUN=function(x) 
      crop(x, y = studyarea_sin))
    
    # Load quality layer and extract Bits
    quality <- stack(preStack(path=path, pattern="*VI_Quality.tif$"))
    beginCluster(n)
    quality <- clusterR(quality, fun=function(y) y/10000)
    endCluster()
    quality <- crop(quality, studyarea_sin)
    
    # Clean according to quality
    stack_clean <- lapply(crop_list, FUN=function(x) 
      overlay(x, y=quality, fun = f_cleanvi, unstack=FALSE))
    
    # Give name to layers  
    dates <- lapply(vi, FUN=function(x) 
      extractDate(preStack(path=path, pattern=paste0(x, ".tif$")), 
                  asDate=TRUE)[[1]])
    stack_clean <- mapply(setZ, stack_clean, dates)
    
    # Project list of raster stacks
    beginCluster(n)
    stack_laea <- lapply(stack_clean, projectRaster, method="bilinear", crs=crs.laea)
    endCluster()
    names <- lapply(vi, FUN=function(x) paste0(workdir, "Results/M.D13_", x, "_", isos[i], ".grd"))
    # And save to file
    mapply(FUN=function(x,y) writeRaster(x=x, filename=y, bandorder = "BIL", 
                                         options = c("COMPRESS=NONE"), dataType = "INT2U", 
                                         overwrite=TRUE), x=stack_laea, y = names)
  }
}

## ---- modis_sur_refl ----
for(i in 1:2){
  # Read study area shapefile 
  studyarea <- readOGR(dsn = "Results/", layer = paste0("studyarea_", isos[i]))
  
  if(!file.exists(paste0(workdir, "Results/M.D09_b01_", isos[i], "_brick.grd"))){
    # Transform studyarea shapefile to sinusoidal projection
    studyarea_sin <- spTransform(studyarea, CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 
                                                +b=6371007.181 +units=m +no_defs")); rm(studyarea)
    
    # Specify start and end period
    begin <- c("2007.10.01", "2005.06.01")
    end <- c("2009.05.02", "2008.06.25")
    
    if(!dir.exists(paste0(isos[i], "_Sur_Refl"))){
      # Download data
      modis.hdf <- getHdf(product = "M.D09A1", begin = begin[i], end = end[i], 
                          extent=studyarea_sin, checkIntegrity=TRUE) #Extract MOD09A1 and MYD09A1 products
      
      # Extract GeoTiff for each time stamp and SDS Layer specified
      runGdal(job=paste0(isos[i], "_Sur_Refl"), product = "M.D09A1", begin = begin[i], end = end[i], 
              extent=studyarea_sin, SDSstring="1111111000010", checkIntegrity=TRUE)
    }
    
    # Specify path of files
    path <- paste0(options("MODIS_outDirPath"), isos[i], "_Sur_Refl")
    
    # Load band stacks and convert to reflectance (?)
    bands <- c("b01", "b02", "b03", "b04", "b05", "b06", "b07")
    preStack_list <- lapply(bands, FUN=function(x) preStack(path=path, pattern=paste0(x, ".tif$")))
    stack_list <- lapply(preStack_list, stack)
    beginCluster(n)
    stack_list <- lapply(stack_list, FUN=function(x){clusterR(x, fun=function(y) y/10000)})
    endCluster()

    # Give name to layers
    dates <- lapply(c(1:7), FUN=function(x) extractDate(preStack(path=path, pattern=paste0(bands, ".tif$")), asDate=TRUE)[[1]])
    stack_list <- mapply(setZ, stack_list, dates)
    
    # Project list of raster stacks
    beginCluster(n)
    bands_laea <- lapply(stack_list, projectRaster, method="bilinear", crs=crs.laea)
    endCluster()
    names <- lapply(bands, FUN=function(x) paste0(workdir, "Results/M.D09_", x, "_", isos[i], "_brick.grd"))
    # And save to file
    mapply(FUN=function(x,y) writeRaster(x=x, filename=y, bandorder = "BIL", options = c("COMPRESS=NONE"), dataType = "INT2U", overwrite=TRUE), x=bands_laea, y = names)
  }
}

## ---- modis_vi_wi ----
for(i in 1:length(isos)){
  # Read study area shapefile 
  studyarea_laea <- readOGR(dsn = "Results/", layer = paste0("studyarea_laea_", isos[i]), verbose=FALSE)
  
  if(!file.exists(paste0(workdir, "Results/M.D09_ndvi_", isos[i], ".grd"))){
    band1 <- stack(paste0(workdir, "Results/M.D09_b01_", isos[i], "_brick.grd"))
    band2 <- stack(paste0(workdir, "Results/M.D09_b02_", isos[i], "_brick.grd"))
    ndvi <- overlay(band1, band2, fun = function(b1, b2){return((b2-b1)/(b2+b1))}, unstack=FALSE, 
                    filename=paste0("Results/M.D09_ndvi_", isos[i], ".grd"), bandorder = "BIL", 
                    options = c("COMPRESS=NONE"), dataType = "INT2U"); rm(band1, band2) # Rouse 1974
  }
  
  if(!file.exists(paste0(workdir, "Results/M.D09_evi_", isos[i], ".grd"))){
    band1 <- stack(paste0(workdir, "Results/M.D09_b01_", isos[i], "_brick.grd"))
    band2 <- stack(paste0(workdir, "Results/M.D09_b02_", isos[i], "_brick.grd"))
    band3 <- stack(paste0(workdir, "Results/M.D09_b03_", isos[i], "_brick.grd"))
    evi <- overlay(band3, band1, band2, fun = function(b3, b1, b2){return(2.5*((b2-b1)/(b2 + 6*b1 - 7.5*b3 + 1)))}, 
                   unstack=FALSE, filename=paste0("Results/M.D09_evi_", isos[i], ".grd"), 
                   bandorder = "BIL", options = c("COMPRESS=NONE"), dataType = "INT2U")
    rm(band1, band2, band3) # Huete 1999
  }
  
  if(!file.exists(paste0(workdir, "Results/M.D09_ndwi_mcf_", isos[i], ".grd"))){
    band4 <- stack(paste0(workdir, "Results/M.D09_b04_", isos[i], "_brick.grd"))
    band2 <- stack(paste0(workdir, "Results/M.D09_b02_", isos[i], "_brick.grd"))
    ndwi_mcf <- overlay(band4, band2, fun = function(b4, b2){(b4-b2)/(b4+b2)}, 
                        unstack=FALSE, filename=paste0("Results/M.D09_ndwi_mcf_", isos[i], ".grd"), 
                        bandorder = "BIL", options = c("COMPRESS=NONE"), dataType = "INT2U")
    rm(band2, band4) # Gao 1996, McFeeters
  }
  
  if(!file.exists(paste0(workdir, "Results/M.D09_ndwi_xu_", isos[i], ".grd"))){
    band6 <- stack(paste0(workdir, "Results/M.D09_b06_", isos[i], "_brick.grd"))
    band4 <- stack(paste0(workdir, "Results/M.D09_b04_", isos[i], "_brick.grd"))
    ndwi_xu <- overlay(band4, band6, fun = function(b4, b6){(b4-b6)/(b4+b6)}, 
                       unstack=FALSE, filename=paste0("Results/M.D09_ndwi_xu_", isos[i], ".grd"), 
                       bandorder = "BIL", options = c("COMPRESS=NONE"), dataType = "INT2U")
    rm(band4, band6) # Xu 2006
  }
  
  if(!file.exists(paste0(workdir, "Results/M.D09_ndmi_", isos[i], ".grd"))){
    band6 <- stack(paste0(workdir, "Results/M.D09_b06_", isos[i], "_brick.grd"))
    band2 <- stack(paste0(workdir, "Results/M.D09_b02_", isos[i], "_brick.grd"))
    ndmi <- overlay(band2, band6, fun = function(b2, b6){(b2-b6)/(b2+b6)}, 
                    filename=paste0("Results/M.D09_ndmi_", isos[i], ".grd"), 
                    bandorder = "BIL", options = c("COMPRESS=NONE"), 
                    dataType = "INT2U")
    rm(band2, band6) # 
  }
  
  if(!file.exists(paste0(workdir, "Results/M.D09_awei_ns_", isos[i], ".grd"))){
    band6 <- stack(paste0(workdir, "Results/M.D09_b06_", isos[i], "_brick.grd"))
    band2 <- stack(paste0(workdir, "Results/M.D09_b02_", isos[i], "_brick.grd"))
    band4 <- stack(paste0(workdir, "Results/M.D09_b04_", isos[i], "_brick.grd"))
    awei_ns <- overlay(band4, band2, band6, fun = 
                         function(b4, b2, b6){4*(b4 - b6)- (0.25*b2+2.75*b6)}, 
                       unstack=F, filename=paste0("Results/M.D09_awei_ns_", isos[i], ".grd"), 
                       bandorder = "BIL", options = c("COMPRESS=NONE"), dataType = "INT2U")
  }
  
  if(!file.exists(paste0(workdir, "Results/M.D09_awei_sh_", isos[i], ".grd"))){
    band3 <- stack(paste0(workdir, "Results/M.D09_b03_", isos[i], "_brick.grd"))
    band7 <- stack(paste0(workdir, "Results/M.D09_b07_", isos[i] , "_brick.grd"))
    band2 <- stack(paste0(workdir, "Results/M.D09_b02_", isos[i], "_brick.grd"))
    band4 <- stack(paste0(workdir, "Results/M.D09_b04_", isos[i] , "_brick.grd"))
    band6 <- stack(paste0(workdir, "Results/M.D09_b06_", isos[i], "_brick.grd"))
    awei_sh <- overlay(band3, band4, band2, band6, band7, fun= 
                         function(b3, b4, b2, b6, b7){b3 + 2.5*b4 - 1.5*(b2 + b6) - 0.25*b7}, 
                       unstack=F, filename=paste0("Results/M.D09_awei_sh_", isos[i], ".grd"), 
                       bandorder = "BIL", options = c("COMPRESS=NONE"), dataType = "INT2U"); rm(band3, band7)
  }
  
  if(!file.exists(paste0(workdir, "Results/M.D09_wi_2015_", isos[i], ".grd"))){
    band1 <- stack(paste0(workdir, "Results/M.D09_b01_", isos[i], "_brick.grd"))
    band7 <- stack(paste0(workdir, "Results/M.D09_b07_", isos[i], "_brick.grd"))
    band2 <- stack(paste0(workdir, "Results/M.D09_b02_", isos[i], "_brick.grd"))
    band4 <- stack(paste0(workdir, "Results/M.D09_b04_", isos[i] , "_brick.grd"))
    band6 <- stack(paste0(workdir, "Results/M.D09_b06_", isos[i], "_brick.grd"))
    wi_2015 <- overlay(band4, band1, band2, band6, band7, fun = 
                         function(b4, b1, b2, b6, b7){1.7204 + 171*b4 + 3*b1 - 70*b2 - 45*b6 - 71*b7}, 
                       unstack=F, filename=paste0("Results/M.D09_wi_2015_", isos[i], ".grd"), 
                       bandorder = "BIL", options = c("COMPRESS=NONE"), dataType = "INT2U"); 
    rm(band1, band2, band4, band6, band7)
  }
}

## ---- ls_vi_wi ----
for(i in 1:length(isos)){
  # Read study area shapefile 
  studyarea_laea <- readOGR(dsn = "Results/", layer = paste0("studyarea_laea_", isos[i]), verbose=FALSE)
  
  if(!file.exists(paste0(workdir, "Results/LS_ndvi_", isos[i], ".grd"))){
    band3 <- stack(paste0(workdir, "Results/LS_b03_", isos[i], "_brick.grd"))
    band4 <- stack(paste0(workdir, "Results/LS_b04_", isos[i], "_brick.grd"))
    ndvi <- overlay(band3, band4, fun = function(b3, b4){return((b4-b3)/(b4+b3))}, 
                    unstack=FALSE, filename=paste0("Results/LS_ndvi_", isos[i], ".grd"), 
                    bandorder = "BIL", options = c("COMPRESS=NONE"), dataType = "INT2U") # Rouse 1974
  }
  
  if(!file.exists(paste0(workdir, "Results/LS_evi_", isos[i], ".grd"))){
    band1 <- stack(paste0(workdir, "Results/LS_b01_", isos[i], "_brick.grd"))
    band3 <- stack(paste0(workdir, "Results/LS_b03_", isos[i], "_brick.grd"))
    band4 <- stack(paste0(workdir, "Results/LS_b04_", isos[i], "_brick.grd"))
    evi <- overlay(band1, band3, band4, fun = function(b1, b3, b4){return(2.5*((b4-b3)/(b4 + 6*b3 - 7.5*b1 + 1)))}, unstack=FALSE, filename=paste0("Results/LS_evi_", isos[i], ".grd"), bandorder = "BIL", options = c("COMPRESS=NONE"), dataType = "INT2U") # Huete 1999
  }
  
  if(!file.exists(paste0(workdir, "Results/LS_ndwi_mcf_", isos[i], ".grd"))){
    band2 <- stack(paste0(workdir, "Results/LS_b02_", isos[i], "_brick.grd"))
    band4 <- stack(paste0(workdir, "Results/LS_b04_", isos[i], "_brick.grd"))
    ndwi_mcf <- overlay(band2, band4, fun = function(b2, b4){(b2-b4)/(b2+b4)}, unstack=FALSE, filename=paste0("Results/LS_ndwi_mcf_", isos[i], ".grd"), bandorder = "BIL", options = c("COMPRESS=NONE"), dataType = "INT2U") # Gao 1996
  }
  
  if(!file.exists(paste0(workdir, "Results/LS_ndwi_xu_", isos[i], ".grd"))){
    band2 <- stack(paste0(workdir, "Results/LS_b02_", isos[i], "_brick.grd"))
    band5 <- stack(paste0(workdir, "Results/LS_b05_", isos[i], "_brick.grd"))
    ndwi_xu <- overlay(band2, band5, fun = function(b2, b5){(b2-b5)/(b2+b5)}, unstack=FALSE, filename=paste0("Results/LS_ndwi_xu_", isos[i], ".grd"), bandorder = "BIL", options = c("COMPRESS=NONE"), dataType = "INT2U")
  }
  
  if(!file.exists(paste0(workdir, "Results/LS_ndmi_", isos[i], ".grd"))){
    band4 <- stack(paste0(workdir, "Results/LS_b04_", isos[i], "_brick.grd"))
    band5 <- stack(paste0(workdir, "Results/LS_b05_", isos[i], "_brick.grd"))
    ndmi <- overlay(band4, band5, fun = function(b4, b5){(b4-b5)/(b4+b5)}, 
                    filename=paste0("Results/LS_ndmi_", isos[i], ".grd"), bandorder = "BIL", options = c("COMPRESS=NONE"), dataType = "INT2U")
  }
  
  if(!file.exists(paste0(workdir, "Results/LS_awei_ns_", isos[i], ".grd"))){
    band2 <- stack(paste0(workdir, "Results/LS_b02_", isos[i], "_brick.grd"))
    band4 <- stack(paste0(workdir, "Results/LS_b04_", isos[i], "_brick.grd"))
    band5 <- stack(paste0(workdir, "Results/LS_b05_", isos[i], "_brick.grd"))
    awei_ns <- overlay(band2, band4, band5, fun = function(b2, b4, b5){4*(b2 - b5)- (0.25*b4+2.75*b5)}, unstack=F, filename=paste0("Results/LS_awei_ns_", isos[i], ".grd"), bandorder = "BIL", options = c("COMPRESS=NONE"), dataType = "INT2U")
  }
  
  if(!file.exists(paste0(workdir, "Results/LS_awei_sh_", isos[i], ".grd"))){
    band1 <- stack(paste0(workdir, "Results/LS_b01_", isos[i], "_brick.grd"))
    band2 <- stack(paste0(workdir, "Results/LS_b02_", isos[i], "_brick.grd"))
    band4 <- stack(paste0(workdir, "Results/LS_b04_", isos[i], "_brick.grd"))
    band5 <- stack(paste0(workdir, "Results/LS_b05_", isos[i], "_brick.grd"))
    band7 <- stack(paste0(workdir, "Results/LS_b07_", isos[i] , "_brick.grd"))
    awei_sh <- overlay(band1, band2, band4, band5, band7, fun= function(b1, b2, b4, b5, b7){b1 + 2.5*b2 - 1.5*(b4 + b5) - 0.25*b7}, unstack=F, filename=paste0("Results/LS_awei_sh_", isos[i], ".grd"), bandorder = "BIL", options = c("COMPRESS=NONE"), dataType = "INT2U"); rm(band1, band7)
  }
  
  if(!file.exists(paste0(workdir, "Results/LS_wi_2015_", isos[i], ".grd"))){
    band2 <- stack(paste0(workdir, "Results/LS_b02_", isos[i], "_brick.grd"))
    band3 <- stack(paste0(workdir, "Results/LS_b03_", isos[i], "_brick.grd"))
    band4 <- stack(paste0(workdir, "Results/LS_b04_", isos[i], "_brick.grd"))
    band5 <- stack(paste0(workdir, "Results/LS_b05_", isos[i], "_brick.grd"))
    band7 <- stack(paste0(workdir, "Results/LS_b07_", isos[i], "_brick.grd"))
    wi_2015 <- overlay(band2, band3, band4, band5, band7, fun = function(b2, b3, b4, b5, b7){1.7204 + 171*b2 + 3*b3 - 70*b4 - 45*b5 - 71*b7}, unstack=F, filename=paste0("Results/LS_wi_2015_", isos[i], ".grd"), bandorder = "BIL", options = c("COMPRESS=NONE"), dataType = "INT2U"); rm(band2, band3, band4, band5, band7)
  }
}

## ---- LS_tasseledCap ----
for(i in 1:2){
  if(!file.exists(paste0(workdir, "Results/LS_brightness_", isos[i], ".grd")) | !file.exists(paste0(workdir, "Results/LS_greenness_", isos[i], ".grd")) | !file.exists(paste0(workdir, "Results/LS_wetness_", isos[i], ".grd"))){
    
    # Read study area shapefile 
    studyarea_laea <- readOGR(dsn = "Results/", layer = paste0("studyarea_laea_", isos[i]), verbose=FALSE)
    
    # Get file names
    ls_list <- list.files(path = paste0(workdir, "Results/LS_Merged"), pattern = isos[i], full.names=TRUE)
    ls_list <- ls_list[c(grep(x=ls_list, pattern=".grd"))]
    
    # Get dates
    ls_dates <- lapply(ls_list, FUN= function(x) unlist(strsplit(unlist(strsplit(x, split="_"))[4], split=".", fixed=TRUE))[1])
    ls_dates <- lapply(ls_dates, function(x) as.Date(as.numeric(substring(x, first=5)) - 1, origin = paste0(substring(x, first=1, last=4), "-01-01")))
    ls_dates <- do.call("c", ls_dates)
    
    # Load raster stacks
    ls_stack <- lapply(ls_list, stack)
    
    # Split list according to sensor
    ls_le7 <- ls_stack[c(grep(x=ls_list, pattern="LE7"))]
    ls_lt5 <- ls_stack[c(grep(x=ls_list, pattern="LT5"))]
    
    # Run Tassel Cap Transformation
    beginCluster(n)
    le7_tc <- lapply(ls_le7, tasseledCap, sat="Landsat7ETM")
    lt5_tc <- lapply(ls_lt5, tasseledCap, sat="Landsat5TM")
    endCluster()
    ls_stack_tc <- c(le7_tc, lt5_tc)
    
    tc <- c("brightness", "greenness", "wetness")
    for(j in 1:length(tc)){
      if(!file.exists(paste0(workdir, "Results/LS_", tc[j], "_", isos[i], ".grd"))){
        # Transpose stack from bands to time series
        ls_tc <- stack(lapply(ls_stack_tc, function(x) x[[j]]))
        # Set names
        ls_tc <- setZ(ls_tc, z=ls_dates, name="date")
        names(ls_tc) <- getZ(ls_tc)
        # Save Bands to file
        writeRaster(ls_tc, filename = paste0(workdir, "Results/LS_", tc[j], "_", isos[i], ".grd"), bandorder = "BIL", options = c("COMPRESS=NONE"), dataType = "INT2U", overwrite=TRUE)
        removeTmpFiles(h=0.01)
      }
    }
  }
}

## ---- MODIS_tasselCap ----
for(i in 1:length(isos)){
  if(!file.exists(paste0(workdir, "Results/M.D09_brightness_", isos[i], ".grd")) | 
     !file.exists(paste0(workdir, "Results/M.D09_greenness_", isos[i], ".grd")) | 
     !file.exists(paste0(workdir, "Results/M.D09_wetness_", isos[i], ".grd"))){
    # Read study area shapefile 
    studyarea_laea <- readOGR(dsn = "Results/", layer = 
                                paste0("studyarea_laea_", isos[i]))
    
    # Specify path of files
    path <- paste0(options("MODIS_outDirPath"), isos[i], "_Sur_Refl")
    
    # Load band stacks and convert to reflectance (?)
    bands <- c("b01", "b02", "b03", "b04", "b05", "b06", "b07")
    preStack_list <- lapply(bands, FUN=function(x) 
      preStack(path=path, pattern=paste0(x, ".tif$")))
    
    # Get dates
    modis_dates <- lapply(preStack_list[[1]], FUN= function(x) 
      unlist(strsplit(unlist(strsplit(x, split=".", fixed=TRUE))[2], 
                      split=".", fixed=TRUE))[1])
    modis_dates <- lapply(modis_dates, function(x) 
      as.Date(as.numeric(substring(x, first=6)) - 1, origin = 
                paste0(substring(x, first=2, last=5), "-01-01")))
    modis_dates <- do.call("c", modis_dates)
    
    # Save modis dates to file
    saveRDS(modis_dates, paste0(workdir, "/Results/modis_dates_", isos[i], ".rds"))
    
    # Transpose stack & apply conversion factor
    preStack_list_t <- purrr::transpose(preStack_list)
    stack_list <- lapply(preStack_list_t, stack)
    beginCluster(n)
    stack_list <- lapply(stack_list, FUN=
                           function(x){clusterR(x, fun=function(y) y/10000)})
    endCluster()
    
    # Project list of raster stacks
    beginCluster(n)
    stack_list <- lapply(stack_list, projectRaster, method="bilinear", crs=crs.laea)
    endCluster()
    
    # Run Tassel Cap Transformation
    beginCluster(n)
    modis_stack_tc <- lapply(stack_list, tasseledCap, sat="MODIS")
    endCluster()
    
    tc <- c("brightness", "greenness", "wetness")
    for(j in 1:length(tc)){
      if(!file.exists(paste0("Results/M.D09_", tc[j], "_", isos[i], ".grd"))){
        # Transpose stack from bands to time series
        modis_tc <- stack(lapply(modis_stack_tc, FUN=function(x){x[[j]]}))
        # Set names
        modis_tc <- setZ(modis_tc, z=modis_dates, name="date")
        names(modis_tc) <- getZ(modis_tc)
        # Save Bands to file
        writeRaster(modis_tc, filename = paste0(workdir, "Results/M.D09_", tc[j], "_", isos[i], ".grd"), bandorder = "BIL", options = c("COMPRESS=NONE"), dataType = "INT2U", overwrite=TRUE)
        removeTmpFiles(h=0.01)
      }
    }
  }
}


## ---- modis_vi_plot ----
# Load VI satellist list
vi_names <- c("M.D09_evi", "M.D09_ndvi", "M.D09_ndwi_mcf", "M.D09_ndmi", 
              "M.D09_ndwi_xu", "M.D09_awei_ns", "M.D09_awei_sh", "M.D09_wi_2015", 
              "M.D09_brightness", "M.D09_greenness", "M.D09_wetness", 
              "M.D13_NDVI", "M.D13_EVI")
m.d09_names <- c("M.D09_evi", "M.D09_ndvi", "M.D09_ndwi_mcf", "M.D09_ndmi", 
                 "M.D09_ndwi_xu", "M.D09_awei_ns", "M.D09_awei_sh", "M.D09_wi_2015", 
                 "M.D09_brightness", "M.D09_greenness", "M.D09_wetness")
m.d13_names <- c("M.D13_NDVI", "M.D13_EVI")

for(i in 1:length(isos)){
  # Read study area shapefile 
  studyarea_laea <- readOGR(dsn = "Results/", 
                            layer = paste0("studyarea_laea_", isos[i]))
  
  # Load files 
  m.d09 <- lapply(m.d09_names, FUN=function(x) 
    paste0(workdir, "Results/", x, "_", isos[i], ".grd"))
  m.d09_stack <- lapply(m.d09, stack)
  
  # Read modis dates from file
  m.d09_dates <- readRDS(paste0(workdir, "/Results/modis_dates_", isos[i], ".rds"))
  
  # Set date and time
  m.d09_stack <- lapply(m.d09_stack, FUN=function(x) setZ(x, z=m.d09_dates, name="date"))
  
  # Load M.13 Data
  m.d13 <- lapply(m.d13_names, FUN=function(x) 
    paste0(workdir, "Results/", x, "_", isos[i], ".grd"))
  m.d13_stack <- lapply(m.d13, stack)
  
  # Specify path of files
  library(MODIS)
  MODISoptions(localArcPath = paste0(filedir, "MODIS"), 
               outDirPath = paste0(filedir, "MODIS/Processed"))
  path <- paste0(options("MODIS_outDirPath"), isos[i], "_VI")
  # Give name to layers 
  m.d13_dates <- extractDate(preStack(path=path, 
                                      pattern=paste0("NDVI.tif$")), 
                             asDate=TRUE)[[1]]
  m.d13_stack <- lapply(m.d13_stack, FUN=function(x) 
    setZ(x, z=m.d13_dates, name="date"))
  
  # Merge M.D09 and M.D13
  vi_stack <- c(m.d09_stack, m.d13_stack)
  
  # Calculate mean and sd over time
  if(!file.exists(paste0("Figures/MODIS_sum_", isos[i], ".pdf"))){
    vi_sum <- lapply(vi_stack, FUN=function(x) 
      data.frame(date = getZ(x), mean = 
                   cellStats(x, mean), sd = 
                   cellStats(x, sd)))
    vi_sum <- do.call("rbind", vi_sum)
    vi_sum$vi <- do.call("c", 
                         mapply(x=c(1:length(vi_names)), y=
                                  lapply(vi_stack, FUN=function(x) 
                                    length(getZ(x))), FUN=function(x,y) 
                                      rep(x, time= y)))
    vi_sum$vi <- factor(vi_sum$vi, labels=vi_names)
    ggplot(data=vi_sum, aes(x=date)) + geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), fill = "grey70") + 
      geom_line(aes(y=mean)) + facet_wrap(~ vi, scales="free_y", ncol=2) + 
      labs(x="Date", y="Mean ± SD")
    ggsave(paste0("Figures/MODIS_sum_", isos[i], ".pdf"), width=8, height=12)
  }
  
  # Plot VI maps for selected dates
  names <- c("evi", "ndvi", "ndwi_mcf", "ndmi", "ndwi_xu", 
             "awei_ns", "awei_sh", "wi_2015", "brightness", 
             "greenness", "wetness", "NDVI", "EVI")
  for(k in 1:length(vi_stack)){
    vi_ind <- vi_stack[[k]]
    names(vi_ind) <- getZ(vi_ind)
    vi_sub <- as.data.frame(rasterToPoints(vi_ind[[round(seq(1, nlayers(vi_ind), length=6))]]))
    vi_dates <- rep(getZ(vi_ind[[round(seq(1, nlayers(vi_ind), length=6))]]), each=nrow(vi_sub))
    vi_sub <- tidyr::gather(vi_sub, date, value, -c(x,y))
    vi_sub$date <- as.Date(vi_dates)
    if(!file.exists(paste0("Figures/MODIS_", names[k], "_", isos[i], ".pdf"))){
      ggplot(data = vi_sub, aes(x = x, y = y), maxpixels=maxpixels) + 
        geom_raster(aes(fill=value)) + facet_wrap(~ date) + 
        scale_fill_gradientn(name=names[k], colors = rev(terrain.colors(255)), 
                             na.value="transparent")  + 
        labs(x="Eastings (m)", y="Northings (m)") + 
        geom_polygon(data=studyarea_laea, aes(x=long, y=lat, group=group), 
                     fill="transparent", color="red") + 
        theme(legend.background = element_rect(fill="transparent", colour=NA), 
              axis.text.y = rotatedAxisElementText(90,"y")) + 
        coord_equal(expand=FALSE)
      if(i == 1){
        ggsave(paste0("Figures/MODIS_", names[k], "_", isos[i], ".pdf"), width=12, height=6)
      } else if(i == 2){
        ggsave(paste0("Figures/MODIS_", names[k], "_", isos[i], ".pdf"), width=7, height=7)
      }
    }
  }
}    

## ---- ls_vi_plot ----
# Load VI list
vi_names <- c("evi", "ndvi", "ndwi_mcf", "ndmi", "ndwi_xu", 
              "awei_ns", "awei_sh", "wi_2015", 
              "brightness", "greenness", "wetness")

ls_dates <- list(BWA = c("2007-10-09", "2007-10-25", "2007-11-10", 
                         "2007-12-12", "2007-12-28", "2008-01-13", 
                         "2008-02-06", "2008-02-22", "2008-03-09", 
                         "2008-03-25", "2008-04-26", "2008-06-05", 
                         "2008-06-21", "2008-07-07", "2008-07-23", 
                         "2008-08-08", "2008-08-24", "2008-09-09", 
                         "2008-09-25", "2008-10-11", "2008-10-27", 
                         "2008-11-28", "2009-01-15", "2009-01-31", 
                         "2009-02-16", "2009-03-20", "2009-04-05", 
                         "2009-04-21", "2008-01-29", "2008-02-06", 
                         "2008-02-22", "2008-03-09", "2008-03-25", 
                         "2008-04-10", "2008-04-26", "2008-08-08", 
                         "2008-08-24", "2008-09-09", "2008-09-25", 
                         "2009-03-20", "2009-04-21", "2009-05-07", 
                         "2009-05-23"), 
                 KEN = c("2005-05-25", "2005-06-10", "2005-07-12", "2007-05-31", 
                         "2007-07-02", "2007-07-18", "2007-08-03", "2007-08-19", 
                         "2007-09-04", "2007-09-20", "2007-10-06", "2007-10-22", 
                         "2007-11-07", "2007-11-23", "2007-12-09", 
                         "2007-12-25"))
ls_dates <- lapply(ls_dates, FUN=function(x) as.Date(x))

for(i in 1:length(isos)){

# Read study area shapefile 
studyarea_laea <- readOGR(dsn = "Results/", layer = paste0("studyarea_laea_", isos[i]))

# Calculate monthly vi
library(zoo)
if(unique(sapply(vi_names, FUN=function(x) 
  !file.exists(paste0(workdir, "Results/LS_", x, "_monthly_", 
                      isos[i], "_brick.grd"))))==TRUE){
  # Load files 
  vi <- lapply(vi_names, FUN=function(x) paste0(workdir, 
                                                "Results/LS_", x, "_", 
                                                isos[i], ".grd"))
  vi_stack <- lapply(vi, stack)
  
  # Set dates
  vi_stack <- lapply(vi_stack, FUN=function(x) 
    setZ(x, z=ls_dates[[i]], name="date"))
  
  # Calculate monthly index
  vi_month <- lapply(vi_stack, FUN=function(x) 
    zApply(x, by=as.yearmon, fun=mean, name="yearmon"))
  for(k in 1:length(vi_names)){
    if(!file.exists(paste0(workdir, "Results/LS_", 
                           vi_names[k], "_monthly_", 
                           isos[i], "_brick.grd"))){
      writeRaster(x=vi_month[[k]], 
                  filename=paste0(workdir, "Results/LS_", vi_names[k], 
                                  "_monthly_", isos[i], "_brick.grd"), 
                  bandorder = "BIL", options = c("COMPRESS=NONE"), 
                  dataType = "INT2U", overwrite=TRUE)
    }
  }
} else{
  vi <- lapply(vi_names, FUN=function(x) 
    paste0(workdir, "Results/LS_", x, "_monthly_", 
           isos[i], "_brick.grd"))
  vi_month <- lapply(vi, stack)
  vi_dates <- unique(as.yearmon(ls_dates[[i]]))
  vi_month <- lapply(vi_month, FUN=function(x) 
    setZ(x, z=vi_dates, name="date"))
}

# Calculate monthly mean and sd
if(!file.exists(paste0("Figures/LS_vi_monthly_sum_", isos[i], ".pdf"))){
  vi_sum <- lapply(vi_month, FUN=function(x) data.frame(month = getZ(x), 
                                                        mean = cellStats(x, mean), 
                                                        sd = cellStats(x, sd)))
  # cellStats does not work for large files!
  vi_sum <- do.call("rbind", vi_sum)
  vi_sum$vi <- do.call("c", lapply(c(1:length(vi_names)), FUN=function(x) 
    rep(x, length(getZ(vi_month[[1]])))))
  vi_sum$vi <- factor(vi_sum$vi, labels=vi_names)
  vi_sum$month <- as.Date(vi_sum$month)
  ggplot(data=vi_sum, aes(month, mean, ymin = mean-sd, ymax = mean+sd)) + 
    geom_pointrange() + facet_wrap(~ vi, scales="free") + labs(x="Date", y="Mean ± SD")
  ggsave(paste0("Figures/LS_vi_monthly_sum_", isos[i], ".pdf"), width=9, height=7)
}

# Plot monthly vi
for(k in 1:length(vi_month)){
  vi_date <- getZ(vi_month[[k]])
  vi_ind <- fortify(vi_month[[k]], maxpixels=maxpixels)
  vi_ind <- tidyr::gather(vi_ind, month, value, -c(x,y))
  vi_ind$month <- as.Date(vi_date)
  if(!file.exists(paste0("Figures/LS_", vi_names[k], "_monthly_", isos[i], ".pdf"))){
    ggplot(data = vi_ind, aes(x = x, y = y), maxpixels=maxpixels) + 
      geom_raster(aes(fill=value)) + facet_wrap(~ month) + 
      scale_fill_gradientn(name=vi_names[k], colors = rev(terrain.colors(255)), na.value="transparent") + 
      labs(x="Eastings (m)", y="Northings (m)") + 
      geom_polygon(data=studyarea_laea, aes(x=long, y=lat, group=group), 
                   fill="transparent", color="red") + 
      theme(legend.background = element_rect(fill="transparent", colour=NA), 
            axis.text.y = rotatedAxisElementText(90,"y")) + 
      coord_equal(expand=FALSE)
    ggsave(paste0("Figures/LS_", vi_names[k], 
                  "_monthly_", isos[i], ".pdf"), 
           width=9, height=7)
  }
}
}

## ---- landuse_ken ----
# Load land use data
landuse <- readOGR(dsn = "Data/landuse_shp", layer = "Ewaso_LandUse_")
#landuse <- readOGR(dsn = "Data/landuse_shp", layer = "gz_land-use_surface_Project")

# Load study area shapefile
studyarea_laea <- readOGR(dsn = "Results/", layer = "studyarea_laea_KEN")

# Transform into crs.laea projection
landuse <- spTransform(landuse, crs.laea)

# Crop landuse by study area
landuse <- raster::intersect(studyarea_laea, landuse)

# Convert shapefiles into dataframes
landuse@data$id <- rownames(landuse@data)
df_landuse <- fortify(landuse)
df_landuse <- plyr::join(df_landuse,landuse@data, by="id")

# Create table of Landuse classes
print(xtable(data.frame(ID = c(1:nlevels(df_landuse$LANDUSE)), 
                        Class = levels(df_landuse$LANDUSE)), 
             caption=c(paste0("Land use classes of Laikipia, Kenya."), 
                       paste0("Land use classes of Laikipia")), 
             label="table:lu_classes"), include.rownames=FALSE, type="latex", 
      file = paste0("Tables/LU_classes.tex"), 
      caption.placement="top", booktabs=TRUE, table.placement="H")

# Calculate area of each landuse class
luse <- landuse@data
luse$area <- gArea(landuse, byid=TRUE)
library(dplyr)
luse_area <- data.frame(luse %>% group_by(LANDUSE) %>% summarise(total = sum(area)))
luse_area$perc <- (luse_area[,2]/sum(luse_area[,2])*100)
names(luse_area) <- c("Land use", "Area", "Percentage (%)")
luse_area <- na.omit(luse_area)

# Create table
print(xtable(luse_area[,c("Land use", "Percentage (%)")], 
             caption=c(paste0("Percentage area of each land use class of the study area in Kenya."),
                       paste0("Percentage area of each LU class, ", region[2])), 
             label=paste0("table:lu_percentage_", isos[2]), digits=2), 
      type="latex", file = paste0("Tables/LU_percentage_", isos[2], ".tex"), 
      caption.placement="top", include.rownames=FALSE, 
      booktabs=TRUE, table.placement="H")

# Map of Landuse, colour coded by land use type
ggplot() + geom_polygon(data=df_landuse, aes(x=long, y=lat, group=group, fill=factor(LANDUSE)), color="black") + 
  geom_polygon(data=studyarea_laea, aes(x=long, y=lat, group=group), fill="transparent", color="red") + 
  scale_fill_manual(values=viridis::viridis(nlevels(factor(df_landuse$LANDUSE))), name="Landuse") + 
  labs(x="Eastings (m)", y="Northings (m)") + 
  theme(legend.background = element_rect(fill="transparent", colour=NA),
        axis.text.y = rotatedAxisElementText(90,"y")) + coord_equal(expand=FALSE)
ggsave(paste0("Figures/KEN_Landuse.pdf"), width=6, height=6)

# Rasterize Land use to 30m resolution
ls_raster <- raster(resolution = c(30.7, 29.3), extent(studyarea_laea), vals=NA, crs=crs.laea)
landuse_30m <- rasterize(landuse, ls_raster, field=landuse$LANDUSE)
writeRaster(landuse_30m, paste0(workdir, "/Results/landuse_30m_", isos[2], ".tif"), overwrite=TRUE)

# Read zebra movement data
zebra_data_laea <- readRDS(paste0(workdir, "Data/zebra_data_laea_", isos[2], ".rds"))

# Extract land use information for Kenya
zebra_data_laea$landuse <- raster::extract(landuse_30m, zebra_data_laea)

# Save zebra data
saveRDS(zebra_data_laea, paste0(workdir, "Data/zebra_data_laea_", isos[2], ".rds"))

## ---- expl_move ----
for(i in 1:length(isos)){
  
  # Read zebra data
  zebra_data_laea <- readRDS(paste0(workdir, "Data/zebra_data_laea_", isos[i],".rds"))
  if(class(zebra_data_laea) != "MoveStack") {
    zebra_data_laea <- moveStack(zebra_data_laea)
  }
  
  if(i == 1){
    source("Functions/ssf_fun.R")
    step_duration <- 180
    thinned_data <- moveStack(lapply(split(zebra_data_laea), thin_split, step_duration, unit = "mins"))
  }
  
  # Calculate median speed, aver distance, aver duration, averAzimuth
  individualData <- data.frame(unique(zebra_data_laea$individual.local.identifier), 
                               unlist(lapply(speed(zebra_data_laea), mean)), 
                               unlist(lapply(distance(zebra_data_laea), mean)), 
                               unlist(lapply(angle(zebra_data_laea), mean)))
  colnames(individualData) <- c("ID", "Speed", "Step length", "Direction")
  print(xtable(individualData, caption=c("Mean speed (m/s), step length (m) and direction 
               for each individual zebra.", paste0("Movement summary of each individual, ", region[i])), 
        label=paste0("table:summary_zebras_", isos[i])), 
        type="latex", file=paste0("Tables/Summary_zebras_", isos[i], ".tex"), 
        caption.placement="top", include.rownames=FALSE, 
        bookstab=TRUE, table.placement="H")
  
  # Calculate speed
  zebra_data_laea$speed <- unlist(lapply(speed(zebra_data_laea), function(x) c(x, NA)))
  ggplot() + geom_histogram(data = data.frame(zebra_data_laea), aes(x=speed), breaks=seq(0,1.5,0.05)) + 
    labs(x="Speed (m/s)", y="Frequency") + facet_wrap(~ individual.local.identifier) + theme_bw()
  ggsave(paste0(workdir, "Figures/Speed_", isos[i], ".pdf"), width=9, height=7)
  
  # Calculate distance and angle
  zebra_data_laea$angle <- unlist(lapply(angle(zebra_data_laea), function(x) c(x, NA)))
  zebra_data_laea$distance <- unlist(lapply(distance(zebra_data_laea), function(x) c(x, NA)))
  
  # Polar diagram of directions
  ggplot(data = data.frame(angle = unlist(lapply(angle(zebra_data_laea), function(x) c(x, NA))), 
                           individual.local.identifier = trackId(zebra_data_laea)), aes(angle)) + 
    geom_histogram(binwidth=30) + facet_wrap(~ individual.local.identifier) + theme_bw() + 
    coord_polar(theta = "x", start=pi, direction=1) + labs(x="Direction (°)", y="")
  ggsave(paste0("Figures/Zebra_", isos[i], "_Angle_Polar.pdf"), width=6, height=7)
  
  # Plot boxplot of speed within and outside protected area
  status <- c("Protected", "5km Buffer", "Unprotected")
  ggplot(data = data.frame(speed = unlist(lapply(speed(zebra_data_laea), function(x) c(x, NA))), 
                           individual.local.identifier = trackId(zebra_data_laea), 
                           status = zebra_data_laea$status), 
         aes(x=factor(status), y=speed, fill=factor(status))) + 
    geom_boxplot() + scale_x_discrete(labels=status) + scale_y_continuous(limits=c(0,2)) + 
    scale_fill_manual(values=c("forestgreen", "darkseagreen2", "snow2"), breaks=c(0,1,2), labels=status) + 
    labs(x="Status", y="Speed (m/s)") + theme(legend.position = "none") + theme_bw()
  ggsave(paste0("Figures/wdpa_speed_", isos[i], ".pdf"), width=7, height=5)

  # Plot histogram of step length and turning angle
  # Create control locations for each case
  if(i == 1){
    angle_dist <- prepare_angle_dist(thinned_data); summary(angle_dist)
  } else if(i == 2){
    angle_dist <- prepare_angle_dist(zebra_data_laea); summary(angle_dist)
    angle_dist <- data.frame(angle_dist)
    angle_dist <- angle_dist[angle_dist$dist>0,]
  }
  
  # Fit distributions to distances
  library(MASS)
  fexp <- fitdistr(angle_dist[, "dist"], "exponential")
  flogn <- fitdistr(angle_dist[, "dist"], "lognormal")
  theta <- 1/mean(angle_dist[, "dist"]) # Parameter for halfnorm
  
  p1 <- ggplot() + geom_histogram(data=data.frame(angle_dist), 
                                  aes(dist, y=..density..), bins=50) + 
    stat_function(aes(x = 0:5000, color="red"), fun = 
                    function(x) dexp(x, rate = fexp$estimate)) + 
    stat_function(aes(x = 0:5000, color="blue"), fun = 
                    function(x) dlnorm(x, meanlog = flogn$estimate["meanlog"], 
                                       sdlog = flogn$estimate["sdlog"])) + 
    stat_function(aes(x = 0:5000, color="green"), fun = 
                    function(x) fdrtool::dhalfnorm(x, theta = theta)) + 
    ylim(c(0,1e-3)) + xlim(c(0,15000)) + labs(x="Step length (m)", y="Frequency") + 
    scale_colour_manual(name="Distribution", labels = c("exp", "lnorm", "halfnorm"), 
                        values = c("red", "blue", "green")) + theme_bw() +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
  
  # Fit distribution to turning angles
  library(CircStats)
  fkappa <- est.kappa(angle_dist[,"rel.angle"])
  p2 <- ggplot() + geom_histogram(data=data.frame(angle_dist), aes(rel.angle, y=..density..), bins=15) + 
    stat_function(aes(x = -pi:pi), fun = function(x) dvonmises(x, circular(0), kappa=fkappa), color="red") + 
    labs(x="Turning angle", y="Frequency") + theme_bw()
  cowplot::plot_grid(p1, p2, ncol=2, labels = "auto"); rm(p1, p2)
  ggsave(paste0(workdir, "Figures/Angle_Dist_", isos[i], ".pdf"), width=10, height=6, bg="transparent")
  
  # Speed against month
  df <- data.frame(time=month(unlist(timestamps(zebra_data_laea))), 
                   speed = unlist(lapply(speed(zebra_data_laea), function(x) c(x, NA))), 
                   individual.local.identifier = trackId(zebra_data_laea))
  ggplot(data= df, aes(x=factor(time), y=speed)) + 
    geom_boxplot() + 
    scale_y_continuous(limits = quantile(df$speed, c(0.1, 0.9), na.rm=TRUE)) + 
    scale_fill_manual(name="ID", values=colourtheme[1:length(trackId(zebra_data_laea))]) + 
    labs(x="Month", y="Speed (m/s)")
  ggsave(paste0(workdir, "Figures/Speed_Month_", isos[i], ".pdf"), width=6, height=5)
  
  # Define daytime (0600 - 1800) and nighttime (1800 - 0600)
  zebra_data_laea$hour <- hour(timestamps(zebra_data_laea)) + (minute(timestamps(zebra_data_laea))/ 60)
  zebra_data_laea$time <- NA
  zebra_data_laea$time[floor(zebra_data_laea$hour) >= 6 & floor(zebra_data_laea$hour) < 18] <- "Daytime"
  zebra_data_laea$time[floor(zebra_data_laea$hour) < 6] <- "Nighttime"
  zebra_data_laea$time[floor(zebra_data_laea$hour) >= 18] <- "Nighttime"
  zebra_data_laea$time <- as.factor(zebra_data_laea$time)
  
  # Plot speed against daytime and nighttime
  df <- data.frame(time=zebra_data_laea$time, speed = 
                     unlist(lapply(speed(zebra_data_laea), function(x) c(x, NA))), 
             individual.local.identifier = trackId(zebra_data_laea))
  ggplot(data=df, aes(x=factor(time), y=speed)) + 
    geom_boxplot() + 
    scale_y_continuous(limits = quantile(df$speed, c(0.1, 0.9), na.rm=TRUE)) + 
    labs(x="", y="Speed (m/s)")
  ggsave(paste0(workdir, "Figures/Speed_Hour_", isos[i], ".pdf"), width=6, height=5)
  
  # Calculate turning angle
  zebra_data<- readRDS(paste0(workdir, "Data/zebra_data_", isos[i],".rds"))
  zebra_data_laea$turnAngleGc <- unlist(lapply(turnAngleGc(zebra_data), function(x) c(NA, x, NA)))
  zebra_data_laea$abs_turnAngle <- abs(zebra_data_laea$turnAngleGc)
  
  # Plot speed, distance and turning angle over time
  ggplot() + geom_line(data=data.frame(time=timestamps(zebra_data_laea), 
                              speed = unlist(lapply(speed(zebra_data_laea), 
                                                    function(x) c(x, NA))), 
                              individual.local.identifier = trackId(zebra_data_laea)), 
              aes(x=time, y=speed)) + scale_y_continuous(limits=c(0, 2.5)) + 
    facet_wrap(~ individual.local.identifier, ncol=2, scales="free_x") + 
    labs(x="Time", y="Speed (m/s)")
  ggsave(paste0(workdir, "Figures/Time_Speed_", isos[i], ".pdf"), width=7, height=7, bg="transparent")

  # Define dry and rainy season
  if(i == 1){
    zebra_data_laea$cat <- NA
    zebra_data_laea$cat[month(timestamps(zebra_data_laea)) %in% c(4,5,6,7,8,9,10,11)] <- "Dry season"
    zebra_data_laea$cat[month(timestamps(zebra_data_laea)) %in% c(12,1,2,3)] <- "Wet season"
  } else if(i == 2){
    zebra_data_laea$cat <- NA
    zebra_data_laea$cat[month(timestamps(zebra_data_laea)) %in% c(1,2,6,7,8,9)] <- "Dry season"
    zebra_data_laea$cat[month(timestamps(zebra_data_laea)) %in% c(3,4,5,10,11,12)] <- "Wet season"
  }
  
  # Calculate mean NDVI for season
  ndvi_stack <- stack(paste0(workdir, "Results/M.D13_NDVI_", isos[i], ".grd"))
  # Specify path of files
  MODISoptions(localArcPath = paste0(filedir, "MODIS"), outDirPath = paste0(filedir, "MODIS/Processed"))
  path <- paste0(options("MODIS_outDirPath"), isos[i], "_VI")
  # Give name to layers 
  m.d13_dates <- extractDate(preStack(path=path, pattern=paste0("NDVI.tif$")), asDate=TRUE)[[1]]
  ndvi_stack <- setZ(ndvi_stack, z=m.d13_dates, name="date")
  library(lubridate)
  ndvi_month <- zApply(ndvi_stack, by=month, fun=mean, name='month')
  if(i == 1){
    ndvi_season <- overlay(stack(ndvi_month[[c(4,5,6,7,8,9,10,11)]]), fun=mean) # Dry season
    ndvi_season <- addLayer(ndvi_season, overlay(stack(ndvi_month[[c(12,1,2,3)]]), fun=mean)) # Wet season
    names(ndvi_season) <- c("Dry season", "Wet season")
  }else if(i == 2){
    ndvi_season <- overlay(stack(ndvi_month[[c(1,2,6,7,8,9)]]), fun=mean) # Dry season
    ndvi_season <- addLayer(ndvi_season, overlay(stack(ndvi_month[[c(3,4,5,10,11,12)]]), fun=mean)) # Wet season
    names(ndvi_season) <- c("Dry season", "Wet season")
  }
  ndvi_season <- fortify(ndvi_season, maxpixels=maxpixels)
  ndvi_season <- tidyr::gather(ndvi_season, cat, ndvi, -c(x,y))
  ndvi_season$cat <- factor(ndvi_season$cat, labels=c("Dry season", "Wet season"))
  
  # Plot zebra data dry vs. rainy season 
  ggplot() + geom_raster(data=ndvi_season, aes(x=x, y=y, fill=ndvi)) + 
    scale_fill_gradientn(name="NDVI", colors = rev(terrain.colors(255)), na.value="transparent") + 
      geom_point(data = data.frame(zebra_data_laea), aes(x=coords.x1,y=coords.x2), shape=1) + 
      facet_wrap(~ cat) + labs(x="Eastings (m)", y="Northings (m)") + 
    theme(axis.text.y = rotatedAxisElementText(90,"y")) + coord_equal(expand=FALSE)
  if(i == 1){
    ggsave(paste0("Figures/Zebra_", isos[i], "_Dry_vs_Wet.pdf"), width=10, height=4)
  } else if(i == 2){
    ggsave(paste0("Figures/Zebra_", isos[i], "_Dry_vs_Wet.pdf"), width=9, height=7)
  }
 
  # Save move object to file
  saveRDS(zebra_data_laea, paste0(workdir, "Data/zebra_data_laea_", isos[i], ".rds"))
  
  if(i == 1){
    df <- data.frame(zebra_data_laea)
    df$mig <- NA
    df$mig[df$state %in% c(1,3)] <- 0
    df$mig[df$state %in% c(2,4)] <- 1
    library(dplyr)
    df %>% group_by(mig) %>% 
      summarise(mn_speed = mean(speed, na.rm=TRUE), 
                sd_speed = sd(speed, na.rm=TRUE),
                mn_angle = mean(abs_turnAngle, na.rm=TRUE),
                sd_angle = sd(abs_turnAngle, na.rm=TRUE))
    
    df %>% filter(state == 1) %>% # northern range
      summarise(mn_dem = mean(dem, na.rm=TRUE),
                sd_dem = sd(dem, na.rm=TRUE))
    df %>% filter(state == 3) %>% # southern range
      summarise(mn_dem = mean(dem, na.rm=TRUE),
                sd_dem = sd(dem, na.rm=TRUE))
    
    df <- df[,c("state", "speed", "abs_turnAngle")]
    df <- subset(df, speed < 2)
    df <- tidyr::gather(df, var, value, -state)
    df$var <- factor(df$var, labels=c("Absolute Turning Angle", "Speed"))
    ggplot() + geom_boxplot(data=df, aes(x=factor(state), y=value, 
                                         fill=factor(state))) + 
      labs(x="State", y="") + 
      scale_x_discrete(breaks=c(1:4), labels=c("Northern Range", "Southern Migration", 
                                               "Southern Range", "Northern Migration")) +
      facet_wrap(~ var, scales = "free_y") + 
      theme_bw() + theme(legend.position = "none")
    ggsave(paste0(workdir, "Figures/State_Move_", isos[i], ".pdf"), 
           width=10, height=5, bg="transparent")
  }
}

## ---- variogram ----
library(snow)
library(ctmm)
for(i in 1:length(isos)){
  # Read data
  zebra_data_laea <- readRDS(paste0(workdir, "Data/zebra_data_laea_", isos[i],".rds"))
  
  # Convert to telemetry data
  zebra_tele <- lapply(split(zebra_data_laea), as.telemetry)
  
  # Calculate variogram
  if(i == 1){
    zebra_vg <- list()
    # Plot variogram
    cairo_pdf(paste0(workdir, "Figures/Zebra_Vario_6h_", isos[i], ".pdf"), width=6, height=12, bg="transparent")
    par(mfrow=c(4,2), mar=c(4.1,4.1,0.6,0.6))
    for(j in 1:length(zebra_tele)){
      zebra_vg[[j]] <- variogram(zebra_tele[[j]], dt=c(960, rep(3600,6))[j])
      plot(zebra_vg[[j]], fraction= 21600/(nrow(zebra_vg[[j]])*(c(960, rep(3600,6))[j])))
    }
    dev.off()
    # Plot variogram
    cairo_pdf(paste0(workdir, "Figures/Zebra_Vario_", isos[i], ".pdf"), width=6, height=12, bg="transparent")
    par(mfrow=c(4,2), mar=c(4.1,4.1,0.6,0.6))
    for(j in 1:length(zebra_vg)){
      plot(zebra_vg[[j]])
    }
    dev.off()
  } else if(i == 2){
    zebra_vg <- list()
    cairo_pdf(paste0(workdir, "Figures/Zebra_Vario_6h_", isos[i], ".pdf"), width=6, height=12, bg="transparent")
    par(mfrow=c(5,3), mar=c(4.1,4.1,0.6,0.6))
    for(j in 1:length(zebra_tele)){
      zebra_vg[[j]] <- variogram(zebra_tele[[j]], dt=c(rep(3600, 9), rep(480,4))[j])
      # Plot variogram
      plot(zebra_vg[[j]], fraction=21600/(nrow(zebra_vg[[j]])*(c(rep(3600, 9), rep(480,4))[j])))
    }
    dev.off()
    # Plot variogram
    cairo_pdf(paste0(workdir, "Figures/Zebra_Vario_", isos[i], ".pdf"), width=6, height=12, bg="transparent")
    par(mfrow=c(5,3), mar=c(4.1,4.1,0.6,0.6))
    for(j in 1:length(zebra_vg)){
      plot(zebra_vg[[j]])
    }
    dev.off()
  }
  
  # Save variogram file
  saveRDS(zebra_vg, paste0(workdir, "Results/zebra_vg_", isos[i], ".rds"))
}

## ---- periodogram ----
for(i in 1:length(isos)){
  # Read data
  zebra_data_laea <- readRDS(paste0(workdir, "Data/zebra_data_laea_", isos[i],".rds"))
  
  # Convert to telemetry data
  zebra_tele <- lapply(split(zebra_data_laea), as.telemetry)
  
  # If you want to look at these periodicities, use periodogram
  zebra_pg <- lapply(zebra_tele, periodogram)
  
  # Plot periodogram
  cairo_pdf(paste0(workdir, "Figures/Zebra_Perio_", isos[i], ".pdf"), width=4, height=8, bg="transparent")
  if(i == 1){par(mfrow=c(4,2))}else{par(mfrow=c(5,3))}
  for(j in 1:length(zebra_pg)){
    par(mar=c(4.1,4.1,1.1,1.1))
    plot(zebra_pg[[j]])
  }
  dev.off()
}  

## ---- akde_ken ----
library(ctmm)
# Read data
zebra_data_laea <- readRDS(paste0(workdir, "Data/zebra_data_laea_", isos[2],".rds"))
if(class(zebra_data_laea) != "MoveStack") {
  zebra_data_laea <- moveStack(zebra_data_laea)
}

# Convert to telemetry data
zebra_tele <- lapply(split(zebra_data_laea), as.telemetry, projection=crs.laea)

# Facilitate a good starting guess for model fitting
zebra_guess <- lapply(zebra_tele, ctmm.guess, variogram=NULL, interactive=FALSE)

if(!file.exists(paste0(workdir, "Results/zebra_ud_", isos[2], ".rds"))){
  # Select best model
  zebra_select <- rep(list(NULL),length(zebra_guess))
  for(j in 1:length(zebra_guess)){
    zebra_select[[j]] <- ctmm.select(zebra_tele[[j]], zebra_guess[[j]])
  }
  
  # Calculate akde
  ud_zebra <- akde(zebra_tele, zebra_select)
  
  # Read WDPA shapefile
  wdpa_laea <- readOGR(dsn=paste0(workdir, "Results/wdpa_laea_", isos[2]), 
                       layer = "wdpa", verbose=FALSE)
  
  # Plot UD
  ud_zebra_df <- lapply(ud_zebra, FUN=function(x) 
    spTransform(SpatialPolygonsDataFrame.UD(x, level.UD=0.95, level=0.95), 
                crs.laea))
  ud_zebra_df <- do.call("rbind", ud_zebra_df)
  ud_zebra_df <- fortify(ud_zebra_df)
  ud_zebra_df$id <- as.factor(substring(ud_zebra_df$id, first=1, last=4))
  ggplot() + geom_polygon(data=wdpa_laea, aes(x=long, y=lat, group=group), fill="forestgreen") + 
    geom_polygon(data=ud_zebra_df, aes(x=long, y=lat, group=group, colour=id), fill=NA) + 
    labs(x="Eastings (m)", y="Northings (m)") + coord_equal(expand=FALSE) + scale_colour_discrete(name="ID")
  ggsave(paste0(workdir, "Figures/AKDE_", isos[2], ".pdf"), width=6, height=8, bg="transparent")
  
  # Save Akde Result
  saveRDS(ud_zebra, paste0(workdir, "Results/zebra_ud_", isos[2], ".rds"))
  
  # Summary of AKDE
  summary_ud <- data.frame(ID=levels(trackId(zebra_data_laea)), 
                           t(sapply(ud_zebra, FUN=function(x) summary(x)$DOF)), 
                           t(sapply(ud_zebra, FUN=function(x) summary(x)$CI)))
  colnames(summary_ud)[2:6] <- c("DOF.area", "DOF.bandwidth", "CI.low", "CI.ML", "CI.high") 
  # area (square kilometers)
  print(xtable(summary_ud, caption=c("Summary of autocorrelated kernel density 
               estimation for each individual in Kenya.", "Summary of AKDE for each individual"), 
               label=paste0("table:summary_akde_", isos[2])), type="latex", 
        file=paste0("Tables/summary_akde_", isos[2], ".tex"), caption.placement="top", 
        include.rownames=FALSE, bookstab=TRUE, table.placement="H")
  
  # Mean akde per species
  ids <- as.character(unique(zebra_data_laea$individual.local.identifier))
  mean_gz_ud <-mean(ud_zebra[1:7])
  mean_pz_ud <- mean(ud_zebra[8:13])
  
  # Size comparsion
  summary(mean_gz_ud)
  summary(mean_pz_ud)
  
  # Calculate species overlap
  ov <- overlap(ud_zebra); ov[,,"ML"]
  ov_sp <- overlap(list(mean_gz_ud, mean_pz_ud));ov_sp[,,"ML"]
  
  # Plot
  mean_gz_ud_df <- spTransform(
    SpatialPolygonsDataFrame.UD(mean_gz_ud, level.UD=0.95, 
                                level=0.95), crs.laea)
  mean_pz_ud_df <- spTransform(
    SpatialPolygonsDataFrame.UD(mean_pz_ud, level.UD=0.95, 
                                level=0.95), crs.laea)
  ggplot() + geom_polygon(data=wdpa_laea, aes(x=long, y=lat, group=group), fill="forestgreen") + 
    geom_polygon(data=mean_gz_ud_df, aes(x=long, y=lat, group=group, colour="gz"), fill=NA) + 
    geom_polygon(data=mean_pz_ud_df, aes(x=long, y=lat, group=group, colour="pz"), fill=NA) + 
    labs(x="Eastings (m)", y="Northings (m)") + coord_equal(expand=FALSE) + 
    scale_colour_manual(name="Species", breaks=c("gz", "pz"), values=c("#F8766D", "#00BFC4"), 
                          labels=c("Equus grevyi", "Equus burchelli"))
  ggsave(paste0(workdir, "Figures/Mean_UD_Sp_", isos[2], ".pdf"), width=6, height=7, bg="transparent")
}

## ---- overlap_ken ----

# Read Akde Result
ud_zebra <- readRDS(paste0(workdir, "Results/zebra_ud_", isos[2], ".rds"))

# Get mean UD
mean_zebra_ud <- mean(ud_zebra[1:13])
mean_gz_ud <-mean(ud_zebra[1:7])
mean_pz_ud <- mean(ud_zebra[8:13])

# Read WDPA shapefile
wdpa_laea <- readOGR(dsn=paste0(workdir, "Results/wdpa_laea_", isos[2]), 
                     layer = "wdpa", verbose=FALSE)

# Calculate overlap of shapefiles
mean_zebra_ud_95 <-SpatialPolygonsDataFrame.UD(mean_zebra_ud, level.UD=0.95, level=0.95)
mean_gz_ud_95 <-SpatialPolygonsDataFrame.UD(mean_gz_ud, level.UD=0.95, level=0.95)
mean_pz_ud_95 <-SpatialPolygonsDataFrame.UD(mean_pz_ud, level.UD=0.95, level=0.95)

pa_ud95 <- gIntersection(mean_zebra_ud_95[2,], wdpa_laea)
pa_gz_ud95 <- gIntersection(mean_gz_ud_95[2,], wdpa_laea)
pa_pz_ud95 <- gIntersection(mean_pz_ud_95[2,], wdpa_laea)
nonpa_ud95 <- gDifference(mean_zebra_ud_95[2,], wdpa_laea)

# Proportion of UD protected
(gArea(pa_ud95)/gArea(mean_zebra_ud_95[2,]))*100 #meters
#gArea(nonpa_ud95)/gArea(mean_zebra_ud_95[2,])*100 #meters

# Proportion of UD protected by species
(gArea(pa_gz_ud95)/gArea(mean_gz_ud_95[2,]))*100 #meters
(gArea(pa_pz_ud95)/gArea(mean_pz_ud_95[2,]))*100 #meters


# Read data
zebra_data_laea <- readRDS(paste0(workdir, "Data/zebra_data_laea_", isos[2],".rds"))
if(class(zebra_data_laea) != "MoveStack") {
  zebra_data_laea <- moveStack(zebra_data_laea)
}

# Plot map of difference 
ggplot() + 
  geom_polygon(data=wdpa_laea, aes(x=long, y=lat, group=group), fill="forestgreen") + 
  geom_polygon(data=mean_zebra_ud_95, aes(x=long, y=lat, group=group), colour="blue", fill="NA") + 
  geom_point(data=data.frame(zebra_data_laea), aes(x=coords.x1, y=coords.x2), shape=1) + 
  labs(x="Eastings (m)", y="Northings (m)") + coord_equal(expand=FALSE)
ggsave(paste0(workdir, "Figures/UD_WDPA_Overlap_", isos[2], ".pdf"), 
       width=6, height=8, bg="transparent")

## ---- akde_bwa ----

# Read zebra movement data
zebra_data_laea <- readRDS(paste0(workdir, "Data/zebra_data_laea_", isos[1], ".rds"))
zebra_data_laea <- split(zebra_data_laea)

# Burst zebra data
zebra_burst <- mapply(x=zebra_data_laea, f=states, FUN=function(x,f) burst(x=x, f=f[-1]))

# Split bursts
zebra_bursts <- list()
for(j in 1:length(zebra_burst)){
  burst <- split(zebra_burst[[j]])
  zebra_bursts <- append(zebra_bursts, burst)
}

# Convert to telemetry data
zebra_tele <- lapply(zebra_bursts, as.telemetry)

# Facilitate a good starting guess for model fitting
zebra_guess <- lapply(zebra_tele, ctmm.guess, variogram=NULL, interactive=FALSE)

if(!file.exists(paste0(workdir, "Results/zebra_ud_", isos[1], ".rds"))){
  # Select best model
  zebra_select <- rep(list(NULL),length(zebra_guess))
  for(j in 1:length(zebra_guess)){
    zebra_select[[j]] <- ctmm.select(zebra_tele[[j]], zebra_guess[[j]])
  }
  
  # Calculate akde
  ud_zebra <- akde(zebra_tele, zebra_select)
  
  # Save Akde Result
  saveRDS(ud_zebra, paste0(workdir, "Results/zebra_ud_", isos[1], ".rds"))
}

# Read Akde Result
ud_zebra <- readRDS(paste0(workdir, "Results/zebra_ud_", isos[1], ".rds"))

# Mean akde per burst
#mean_states_ud <-mean(ud_zebra[1:7])

# Plot segment UD
#cairo_pdf(paste0(workdir, "Figures/Mean_States_UD_", isos[1], ".pdf"), width=800, height=800, bg="transparent")
#ctmm::plot(mean_zebra_ud, level.UD=c(0.50, 0.95), col.DF='orange', col.level='orange', col.grid = NA)
#dev.off()

# Load wdpa shapefile

# Calculate overlap of shapefiles
?gIntersection()
?gArea()
?gDifference()

## ---- absence ----
# Load split_thin function
source("Functions/ssf_fun.R")

# Read zebra movement data
zebra_data_laea <- readRDS(paste0(workdir, "Data/zebra_data_laea_", isos[1], ".rds"))  

# Create absence points
if(!file.exists(paste0(workdir, "Data/ssf_data_exp.rds"))){
  step_duration <- 180
  thinned_data <- moveStack(lapply(split(zebra_data_laea), thin_split, step_duration, unit = "mins"))
  my_K <- 100
  set.seed(310)
  theta <- 1/mean(prepare_angle_dist(thinned_data)[, "dist"]) # Parameter for halfnorm
  ssf_data_exp <- prepare_ssf_steps(thinned_data, method = "halfnormal", 
                                    K = my_K, theta = theta/2, crs = crs.laea)
} else{
  ssf_data_exp <- readRDS(paste0(workdir, "Data/ssf_data_exp.rds"))
}

# Plot the choice set for a given step
if(!file.exists(paste0("Figures/SSF_Absence_Test.pdf"))){
  dat <- data.frame(subset(ssf_data_exp, stratum == "Z3864_2_9"))
  dat0 <- data.frame(subset(ssf_data_exp, stratum == "Z3864_2_7" & used == 1))
  dat2 <- data.frame(subset(ssf_data_exp, stratum %in% c("Z3864_2_7", "Z3864_2_8") & used == 1))
  ggplot(data=rbind(dat, dat0), aes(x,y)) + geom_point(shape=1) + 
    geom_point(data=data.frame(
      subset(ssf_data_exp, stratum == "Z3864_2_9" & used == 1)), 
      aes(x,y, colour="t+1")) + 
    geom_point(data=data.frame(
      subset(ssf_data_exp, stratum == "Z3864_2_8" & used == 1)), 
      aes(x,y, colour ="t")) + 
    geom_point(data=dat0, aes(x,y, colour="t-1")) + 
    geom_line(data=dat2, aes(x,y), colour="green") + 
    scale_colour_manual(name="Time step", 
                        values=c("green", "blue", "red"), 
                        labels = c("t-1", "t", "t+1")) + 
    theme_bw() + theme(legend.justification=c(0,1), legend.position=c(0,1)) + 
    labs(x="Eastings (m)", y="Northings (m)") + coord_equal(expand=FALSE)
  ggsave(paste0("Figures/SSF_Absence_Test.pdf"), width=10, height=6, bg="transparent")
}
saveRDS(ssf_data_exp, paste0(workdir, "Data/ssf_data_exp.rds"))

## ---- spt_modis_bwa ----

# Read SSF data
ssf_data_exp <- readRDS(paste0(workdir, "Data/ssf_data_exp.rds"))

# Read different VI stacks
satellites <- c("M.D09", "M.D13")
vi_names <- list(c("evi", "ndvi", "ndwi_mcf", "ndmi", "ndwi_xu", "awei_ns", "awei_sh", 
                   "wi_2015", "brightness", "greenness", "wetness"), c("EVI", "NDVI"))
vis <- mapply(x=satellites, y=vi_names, FUN=function(x,y) 
  lapply(y, FUN=function(y) paste(x, y, sep="_")))
vi <- lapply(vis, FUN=function(x) paste0(workdir, "Results/", x, "_", isos[1], ".grd"))
vi_stack_list <- lapply(vi, FUN=function(x) lapply(x, stack))

# Set time and get time interval
timestamp <- readRDS(paste0(workdir, "Results/modis_dates_", isos[1], ".rds"))
dif <- as.numeric(difftime(timestamp[2], timestamp[1], units="days"), units="days")

# Extract data for step selection function
for(j in 1:length(vi_stack_list)){
  for(k in 1:length(vi_stack_list[[j]])){
    for(l in 1:length(timestamp)){
      ssf_data_exp[[vi_names[[j]][[k]]]][ssf_data_exp$date >= 
                                           timestamp[l]-dif/2 & 
                                           ssf_data_exp$date < timestamp[l]+dif/2] <- 
        extract(vi_stack_list[[j]][[k]], subset(ssf_data_exp, ssf_data_exp$date >= 
                                                  timestamp[l]-dif/2 & 
                                                  ssf_data_exp$date < timestamp[l]+dif/2))
    }
  }  
}

# Save SSF data
saveRDS(ssf_data_exp, paste0(workdir, "Data/ssf_data_exp.rds"))

## ---- spt_mlc_bwa ----
library(lubridate)

# Read SSF data
ssf_data_exp <- readRDS(paste0(workdir, "Data/ssf_data_exp.rds"))

# Read studyarea shapefile
studyarea_laea <- readOGR(dsn = "Results/", layer = paste0("studyarea_laea_", isos[1]))

# Read Landcover data
lcov <- stack(paste0(workdir,"Results/lcov_MCD12Q1_", isos[1], "_brick.grd"))
lcov <- crop(lcov, studyarea_laea)

# Get time and set time interval
timestamp <- do.call("c", lapply(names(lcov), FUN=function(x) as.Date(paste0(substring(x, first=2, last=5), "-01-01"))))
lcov <- setZ(lcov, timestamp, name="date")
timestamp <- year(timestamp)

# Extract spatio-temporal land cover per year
for(j in 1:length(timestamp)){
  ssf_data_exp$lcov[year(ssf_data_exp$date) == timestamp[j]] <- 
    extract(lcov[[j]], subset(ssf_data_exp, year(ssf_data_exp$date) == timestamp[j]))
}

# Save SSF data
saveRDS(ssf_data_exp, paste0(workdir, "Data/ssf_data_exp.rds"))

## ---- spt_ls_bwa ----
# Read SSF data
ssf_data_exp <- readRDS(paste0(workdir, "Data/ssf_data_exp.rds"))

# Load monthly Landsat VI data
vi_names <- c("evi", "ndvi", "ndwi_mcf", "ndmi", "ndwi_xu", "awei_ns", "awei_sh", "wi_2015", "brightness", "greenness", "wetness")
vi <- lapply(vi_names, FUN=function(x) paste0(workdir, "Results/LS_", x, "_monthly_", isos[i], "_brick.grd"))
vi_stack_list <- lapply(vi, stack)

# Specify dates
ls_dates <- list(BWA = c("2007-10-09", "2007-10-25", "2007-11-10", 
                         "2007-12-12", "2007-12-28", "2008-01-13", 
                         "2008-02-06", "2008-02-22", "2008-03-09", 
                         "2008-03-25", "2008-04-26", "2008-06-05", 
                         "2008-06-21", "2008-07-07", "2008-07-23", 
                         "2008-08-08", "2008-08-24", "2008-09-09", 
                         "2008-09-25", "2008-10-11", "2008-10-27", 
                         "2008-11-28", "2009-01-15", "2009-01-31", 
                         "2009-02-16", "2009-03-20", "2009-04-05", 
                         "2009-04-21", "2008-01-29", "2008-02-06", 
                         "2008-02-22", "2008-03-09", "2008-03-25", 
                         "2008-04-10", "2008-04-26", "2008-08-08", 
                         "2008-08-24", "2008-09-09", "2008-09-25", 
                         "2009-03-20", "2009-04-21", "2009-05-07", 
                         "2009-05-23"), 
                 KEN = c("2005-05-25", "2005-06-10", "2005-07-12", "2007-05-31", 
                         "2007-07-02", "2007-07-18", "2007-08-03", "2007-08-19", 
                         "2007-09-04", "2007-09-20", "2007-10-06", "2007-10-22", 
                         "2007-11-07", "2007-11-23", "2007-12-09", 
                         "2007-12-25"))

# Set time stamp
library(zoo)
timestamp <- unique(as.yearmon(ls_dates[[1]]))

# Column names
vi_names <- c("LS_evi", "LS_ndvi", "LS_ndwi_mcf", "LS_ndmi", "LS_ndwi_xu", "LS_awei_ns", 
              "LS_awei_sh", "LS_wi_2015", "LS_brightness", "LS_greenness", "LS_wetness")

# Extract monthly data
for(k in 1:length(vi_stack_list)){
  for(l in 1:length(timestamp)){
    ssf_data_exp[[vi_names[[k]]]][as.yearmon(ssf_data_exp$date) == timestamp[l]] <- 
      extract(vi_stack_list[[k]], subset(ssf_data_exp, as.yearmon(ssf_data_exp$date) == timestamp[l]))
  }
}

# Save SSF data
saveRDS(ssf_data_exp, paste0(workdir, "Data/ssf_data_exp.rds"))
removeTmpFiles(h=0.01)

## ---- spt_modis_ken ----
i <- 2

# Read zebra movement data
zebra_data_laea <- readRDS(paste0(workdir, "Data/zebra_data_laea_", isos[i], ".rds"))
zebra_data_laea <- split(zebra_data_laea)

# Read different VI stacks
satellites <- c("M.D09", "M.D13")
vi_names <- list(c("evi", "ndvi", "ndwi_mcf", "ndmi", "ndwi_xu", "awei_ns", 
                   "awei_sh", "wi_2015", "brightness", "greenness", "wetness"), c("EVI", "NDVI"))
vis <- mapply(x=satellites, y=vi_names, FUN=function(x,y) lapply(y, FUN=function(y) paste(x, y, sep="_")))
vi <- lapply(vis, FUN=function(x) paste0(workdir, "Results/", x, "_", isos[i], ".grd"))
vi_stack_list <- lapply(vi, FUN=function(x) lapply(x, stack))

# Set time and get time interval
timestamp <- readRDS(paste0(workdir, "Results/modis_dates_", isos[i], ".rds"))
dif <- as.numeric(difftime(timestamp[2], timestamp[1], units="days"), units="days")

# Extract data for movement data
for(j in 1:length(vi_stack_list)){
  for(k in 1:length(vi_stack_list[[j]])){
    for(l in 1:length(timestamp)){
      for(m in 1:length(zebra_data_laea)){
        zebra_data_laea[[m]][[vi_names[[j]][[k]]]][timestamps(zebra_data_laea[[m]]) >= 
                                                     timestamp[l]-dif/2 & timestamps(zebra_data_laea[[m]]) < 
                                                     timestamp[l]+dif/2] <- 
          extract(vi_stack_list[[j]][[k]], subset(zebra_data_laea[[m]], timestamps(zebra_data_laea[[m]]) >=
                                                    timestamp[l]-dif/2 & timestamps(zebra_data_laea[[m]]) < 
                                                    timestamp[l]+dif/2))
      }
    }
  }
}

# Save zebra data
saveRDS(zebra_data_laea, paste0(workdir, "Data/zebra_data_laea_", isos[i], ".rds"))

## ---- spt_mlc_ken ----
i <- 2

# Read zebra movement data
zebra_data_laea <- readRDS(paste0(workdir, "Data/zebra_data_laea_", isos[i], ".rds"))

# Read studyarea shapefile
studyarea_laea <- readOGR(dsn = "Results/", layer = paste0("studyarea_laea_", isos[i]))

# Read Landcover data
lcov <- stack(paste0(workdir,"Results/lcov_MCD12Q1_", isos[i], "_brick.grd"))
lcov <- crop(lcov, studyarea_laea)

# Get time and set time interval
timestamp <- do.call("c", lapply(names(lcov), FUN=function(x) 
  as.Date(paste0(substring(x, first=2, last=5), "-01-01"))))
lcov <- setZ(lcov, timestamp, name="date")
timestamp <- year(timestamp)

# Extract land cover per year
for(j in 1:length(timestamp)){
  for(k in 1:length(zebra_data_laea)){
    zebra_data_laea[[k]]$lcov[year(timestamps(zebra_data_laea[[k]])) == timestamp[j]] <- 
      extract(lcov[[j]], subset(zebra_data_laea[[k]], year(timestamps(zebra_data_laea[[k]])) == timestamp[j]))
  }
}

# Save zebra data
saveRDS(zebra_data_laea, paste0(workdir, "Data/zebra_data_laea_", isos[i], ".rds"))

## ---- spt_ls_ken ----
i <- 2

# Read zebra movement data
zebra_data_laea <- readRDS(paste0(workdir, "Data/zebra_data_laea_", isos[i], ".rds"))

# Load monthly Landsat VI data
vi_names <- c("evi", "ndvi", "ndwi_mcf", "ndmi", "ndwi_xu", "awei_ns", "awei_sh", 
              "wi_2015", "brightness", "greenness", "wetness")
vi <- lapply(vi_names, FUN=function(x) paste0(workdir, "Results/LS_", x, "_monthly_", isos[i], "_brick.grd"))
vi_stack_list <- lapply(vi, stack)

# Get dates
ls_dates <- list(BWA = c("2007-10-09", "2007-10-25", "2007-11-10", 
                         "2007-12-12", "2007-12-28", "2008-01-13", 
                         "2008-02-06", "2008-02-22", "2008-03-09", 
                         "2008-03-25", "2008-04-26", "2008-06-05", 
                         "2008-06-21", "2008-07-07", "2008-07-23", 
                         "2008-08-08", "2008-08-24", "2008-09-09", 
                         "2008-09-25", "2008-10-11", "2008-10-27", 
                         "2008-11-28", "2009-01-15", "2009-01-31", 
                         "2009-02-16", "2009-03-20", "2009-04-05", 
                         "2009-04-21", "2008-01-29", "2008-02-06", 
                         "2008-02-22", "2008-03-09", "2008-03-25", 
                         "2008-04-10", "2008-04-26", "2008-08-08", 
                         "2008-08-24", "2008-09-09", "2008-09-25", 
                         "2009-03-20", "2009-04-21", "2009-05-07", 
                         "2009-05-23"), 
                 KEN = c("2005-05-25", "2005-06-10", "2005-07-12", "2007-05-31", 
                         "2007-07-02", "2007-07-18", "2007-08-03", "2007-08-19", 
                         "2007-09-04", "2007-09-20", "2007-10-06", "2007-10-22", 
                         "2007-11-07", "2007-11-23", "2007-12-09", 
                         "2007-12-25"))

# Set time and time interval
library(zoo)
timestamp <- unique(as.yearmon(ls_dates[[2]]))

# Column names
vi_names <- c("LS_evi", "LS_ndvi", "LS_ndwi_mcf", "LS_ndmi", "LS_ndwi_xu", "LS_awei_ns", 
              "LS_awei_sh", "LS_wi_2015", "LS_brightness", "LS_greenness", "LS_wetness")

# Extract monthly data
for(k in 1:length(vi_stack_list)){
  for(l in 1:length(timestamp)){
    for(m in 1:length(zebra_data_laea)){
      zebra_data_laea[[m]][[vi_names[[k]]]][as.yearmon(timestamps(zebra_data_laea[[m]])) == timestamp[l]] <- 
        extract(vi_stack_list[[k]], 
                subset(zebra_data_laea[[m]], as.yearmon(timestamps(zebra_data_laea[[m]])) == timestamp[l]))
    }
  }
}  

# Save zebra data
saveRDS(zebra_data_laea, paste0(workdir, "Data/zebra_data_laea_", isos[i], ".rds"))

## ---- spaceuse_ken ----
# Read data
zebra_data_laea <- readRDS(paste0(workdir, "Data/zebra_data_laea_", isos[2],".rds"))
if(class(zebra_data_laea) != "MoveStack") {
  zebra_data_laea <- moveStack(zebra_data_laea)
}

# Stacked bar chart of land use
landuse <- data.frame(id = c(1:14), 
                      class = factor(c("CCA", "Forest Reserves", "Government Land", 
                                       "Large Scale Farms", "Mukogodo Group Ranches", 
                                       "National Park", "National Reserves", 
             "PWS", "Ranches", "Rhino Sanctuary", "Settlements", "Swamp", 
             "Trust Land", "Urban Settlements"), ordered=TRUE))
df_lu <- lapply(split(zebra_data_laea), FUN=function(x) plyr::count(x$landuse))
df_lu <- do.call("rbind", df_lu)
df_lu$id <- sapply(rownames(df_lu), FUN=function(x) unlist(strsplit(x, split=".", fixed=TRUE))[1])
df_lu$species <- NA
df_lu$species[grep(df_lu$id, pattern="GZ")] <- "GZ"
df_lu$species[grep(df_lu$id, pattern="PZ")] <- "PZ"
df_lu$species <- as.factor(df_lu$species)
df_lu_id <- dplyr::right_join(df_lu, dplyr::count(df_lu, id, wt = freq, by="id"))
(df_lu_id$perc <- df_lu_id$freq/df_lu_id$n*100)

p1 <- ggplot(data=df_lu_id, aes(x=id, y=perc)) + 
  geom_col(aes(fill=factor(x))) + labs(x="ID", y="Percentage (%)") + 
  scale_fill_manual(values=viridis::viridis(nlevels(factor(df_lu_id$x))), 
                    labels=as.character(landuse$class[landuse$id %in% levels(factor(df_lu_id$x))]), 
                    name="Land Use") + theme_bw() + theme(legend.position = "bottom")

classes <- data.frame(id = c(0:16, 254, 255), 
                      class = factor(c("Water", "Evergreen Needleleaf forest", 
                                       "Evergreen Broadleaf forest", 
                                       "Deciduous Needleleaf forest", 
                                       "Deciduous Broadleaf forest", 
                                       "Mixed forest", "Closed shrublands", 
                                       "Open shrublands", "Woody savannas", 
                                       "Savannas", "Grasslands", 
                                       "Permanent wetlands", "Croplands", 
                                       "Urban and built-up", "Natural vegetation mosaic", 
                                       "Snow and ice", "Barren or sparsely vegetated", 
                                       "Unclassified", "Fill Value"), ordered=TRUE), 
                      colours=factor(c("blue", NA, "darkgreen", NA, 
                                       "green2", "green3", "brown", 
                                       "coral1", "lightseagreen", 
                                       "seagreen2", "green", "lightblue", 
                                       "yellow", "grey", "orange", "white", 
                                       "brown2", "transparent", "transparent"), ordered=TRUE))

# Stacked bar chart of land cover
df_lc <- lapply(split(zebra_data_laea), FUN=function(x) plyr::count(x$lcov))
df_lc <- do.call("rbind", df_lc)
df_lc <- df_lc[!is.na(df_lc$x),]
df_lc$id <- sapply(rownames(df_lc), FUN=function(x) unlist(strsplit(x, split=".", fixed=TRUE))[1])
df_lc$species <- NA
df_lc$species[grep(df_lc$id, pattern="GZ")] <- "GZ"
df_lc$species[grep(df_lc$id, pattern="PZ")] <- "PZ"
df_lc$species <- as.factor(df_lc$species)
df_lc_id <- dplyr::right_join(df_lc, dplyr::count(df_lc, id, wt = freq, by="id"))
df_lc_id$perc <- df_lc_id$freq/df_lc_id$n*100
p2 <- ggplot(data=df_lc_id, aes(x=id, y=perc)) + 
  geom_col(aes(fill=factor(x))) + labs(x="ID", y="Percentage (%)") +
  scale_fill_manual(values=as.character(classes$colours[classes$id %in% levels(factor(df_lc_id$x))]), 
                    labels=as.character(classes$class[classes$id %in% levels(factor(df_lc_id$x))]), 
                    name="Land Cover") + theme_bw() + 
  theme(legend.position = "bottom") + guides(fill=guide_legend(nrow=2,byrow=TRUE))
cowplot::plot_grid(p1,p2, labels="auto", ncol=2)
ggsave(paste0("Figures/Zebra_", isos[2], "_LULC.pdf"), width=14, height=6)

# Species comparison
df_lu_sp <- dplyr::right_join(df_lu, dplyr::count(df_lu, species, wt = freq, by="id"))
(df_lu_sp$perc <- df_lu_sp$freq/df_lu_sp$n*100)
p1 <- ggplot(data=df_lu_sp, aes(x=species, y=perc)) + 
  geom_col(aes(fill=factor(x))) + labs(x="Species", y="Percentage (%)") +
  scale_fill_manual(values=viridis::viridis(nlevels(factor(df_lu_id$x))), 
                    labels=as.character(landuse$class[landuse$id %in% levels(factor(df_lu_id$x))]), 
                    name="Land Use") + theme_bw() + theme(legend.position = "bottom")
# Stacked bar chart of land cover per status and species
df_lc_sp <- dplyr::right_join(df_lc, dplyr::count(df_lc, species, wt = freq, by="id"))
df_lc_sp$perc <- df_lc_sp$freq/df_lc_sp$n*100
p2 <- ggplot(data=df_lc_sp, aes(x=species, y=perc)) + 
  geom_col(aes(fill = factor(x))) + labs(x="Species", y="Percentage (%)") + 
  scale_fill_manual(values=as.character(classes$colours[classes$id %in% levels(factor(df_lc_sp$x))]), 
                    labels=as.character(classes$class[classes$id %in% levels(factor(df_lc_sp$x))]), 
                    name="Land Cover") + theme_bw() + 
  theme(legend.position = "bottom") + guides(fill=guide_legend(nrow=2,byrow=TRUE))
cowplot::plot_grid(p1,p2, ncol=2, labels="auto")
ggsave(paste0("Figures/Zebra_Sp_", isos[2], "_LULC.pdf"), width=14, height=6)

# Compare speed and absolute turning angle among land cover and land use classes
df <- data.frame(zebra_data_laea)
df <- df[c("landuse", "speed", "species", "lcov", "abs_turnAngle")]
df <- df[!is.na(df$lcov),]
df <- subset(df, speed < 1)
df$landuse <- factor(df$landuse, labels=c(
  as.character(landuse$class[landuse$id %in% levels(factor(df$landuse))])))
df$lcov <- factor(df$lcov, labels=c(
  as.character(classes$class[classes$id %in% levels(factor(df$lcov))])))
df <- tidyr::gather(df, var, value, -c(landuse, lcov, species))
df <- tidyr::gather(df, lulc, class, -c(var, value, species))
df$lulc <- factor(df$lulc, labels = c("Land Use", "Land Cover"))
df$var <- factor(df$var, labels = c("Absolute turning angle (°)", "Speed (m/s)"))
ggplot() + 
  geom_boxplot(data=df, aes(x=factor(class), y=value, fill=species)) + 
  labs(x="", y="") + 
  scale_fill_discrete(name="Species", breaks=c("GZ", "PZ"), 
                      labels=c("Equus grevyi", "Equus burchelli")) +
  facet_grid(var ~ lulc, scales= "free", drop=TRUE) + theme_bw() + 
  theme(strip.background = element_blank(), legend.position="bottom", 
        axis.text.x = element_text(angle = 45, vjust = 0.5))
ggsave(paste0("Figures/LULC_angle_speed_", isos[2], ".pdf"), width=8, height=6)

## ---- spaceuse_bwa ----
# Read data
ssf_data_exp <- readRDS(paste0(workdir, "Data/ssf_data_exp.rds"))
ssf_data_exp <- subset(ssf_data_exp, used==1)
ssf_data_list <- lapply(levels(ssf_data_exp$id), FUN=function(x) 
  subset(ssf_data_exp, ssf_data_exp$id == x))
names(ssf_data_list) <- levels(ssf_data_exp$id)

classes <- data.frame(id = c(0:16, 254, 255), 
                      class = factor(c("Water", "Evergreen Needleleaf forest", 
                                       "Evergreen Broadleaf forest", 
                                       "Deciduous Needleleaf forest", 
                                       "Deciduous Broadleaf forest", 
                                       "Mixed forest", "Closed shrublands", 
                                       "Open shrublands", "Woody savannas", 
                                       "Savannas", "Grasslands", 
                                       "Permanent wetlands", "Croplands", 
                                       "Urban and built-up", "Natural vegetation mosaic", 
                                       "Snow and ice", "Barren or sparsely vegetated", 
                                       "Unclassified", "Fill Value"), ordered=TRUE), 
                      colours=factor(c("blue", NA, "darkgreen", NA, 
                                       "green2", "green3", "brown", 
                                       "coral1", "lightseagreen", 
                                       "seagreen2", "green", "lightblue", 
                                       "yellow", "grey", "orange", "white", 
                                       "brown2", "transparent", "transparent"), ordered=TRUE))

# Stacked bar chart of land cover
df_lc <- lapply(ssf_data_list, FUN=function(x) plyr::count(x$lcov))
df_lc <- do.call("rbind", df_lc)
df_lc$id <- sapply(rownames(df_lc), FUN=function(x) unlist(strsplit(x, split=".", fixed=TRUE))[1])
df_lc_id <- dplyr::right_join(df_lc, dplyr::count(df_lc, id, wt = freq, by="id"))
df_lc_id$perc <- df_lc_id$freq/df_lc_id$n*100
ggplot(data=df_lc_id, aes(x=id, y=perc)) + 
  geom_col(aes(fill=factor(x))) + labs(x="ID", y="Percentage (%)") +
  scale_fill_manual(values=as.character(classes$colours[classes$id %in% levels(factor(df_lc_id$x))]), 
                    labels=as.character(classes$class[classes$id %in% levels(factor(df_lc_id$x))]), 
                    name="Land Cover") + theme_bw() + 
  theme(legend.position = "bottom")
ggsave(paste0("Figures/Zebra_", isos[1], "_LC.pdf"), width=10, height=5)

# Compare step length and absolute turning angle among land cover classes
df <- data.frame(ssf_data_exp)
df <- df[!is.na(df$lcov),]
p1 <- ggplot() + geom_boxplot(data=df, aes(x=factor(lcov), 
                                           y=dist)) + 
  labs(x = "Land Cover", y = "Step length (m)") + 
  scale_x_discrete(breaks=classes$id[classes$id %in% 
                                       levels(factor(df$lcov))], 
                   labels=as.character(
                     classes$class[classes$id %in% 
                                     levels(factor(df$lcov))])) + 
  theme_bw() +  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 45, vjust = 0.5))

p2 <- ggplot() + 
  geom_boxplot(data=df, aes(x=factor(lcov), 
                            y=abs(radians2degrees(rel.angle)))) + 
  labs(x="Land Cover", y="Absolute Turning angle (°)") + 
  scale_x_discrete(breaks=classes$id[classes$id %in% 
                                       levels(factor(df$lcov))], 
                   labels=as.character(
                     classes$class[classes$id %in% 
                                     levels(factor(df$lcov))])) + 
  theme_bw() + theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 45, vjust = 0.5))
cowplot::plot_grid(p1, p2, ncol=2, labels="auto")
ggsave(paste0("Figures/lcov_angle_dist_", isos[1], ".pdf"), 
       width=12, height=6)

## ---- move_env_ken ----

# Read zebra movement data
zebra_data_laea <- readRDS(paste0(workdir, "Data/zebra_data_laea_", isos[2], ".rds"))
if(class(zebra_data_laea) != "MoveStack"){
  zebra_data_laea <- moveStack(zebra_data_laea)
}

# Plot zebras with NDVI over time
df <- data.frame(zebra_data_laea[order(
  zebra_data_laea$individual.local.identifier, timestamps(zebra_data_laea)),])
df$date <- timestamps(zebra_data_laea)
df$split <- 2007
df$split[df$individual.local.identifier %in% c("PZ6", "PZ8", "PZ10", "PZ14")] <- 2005
ggplot(data=df, aes(x = date, y = ndvi, colour=individual.local.identifier)) + geom_line() + 
  scale_colour_discrete(name="ID") + labs(x="Date", y="NDVI") + theme_bw() + 
  facet_wrap(~ split, scales="free_x")
ggsave(paste0("Figures/Zebra_ndvi_", isos[2], ".pdf"), width=10, height=5)

# Plot boxplot of NDVI for dry vs. rainy season
ggplot(data=data.frame(zebra_data_laea), 
       aes(x=individual.local.identifier, y=ndvi, fill=cat)) + 
  geom_boxplot() + labs(x="ID", y="NDVI") + theme_bw() + 
  scale_fill_manual(name="Season", labels=c("Dry", "Wet"), 
                      values=c("#F8766D", "#00BFC4")) 
ggsave(paste0("Figures/NDVI_", isos[2], "_Dry_vs_Wet.pdf"), width=9, height=7)

zebra_data_laea$month <- lubridate::month(timestamps(zebra_data_laea))
zebra_data_laea$species <- NA
zebra_data_laea$species[grep(zebra_data_laea$individual.local.identifier, pattern="GZ")] <- "GZ"
zebra_data_laea$species[grep(zebra_data_laea$individual.local.identifier, pattern="PZ")] <- "PZ"
zebra_data_laea$species <- as.factor(zebra_data_laea$species)

df <- data.frame(zebra_data_laea)
df <- df[c("month", "species", "ndvi", "evi")] # select any index that is available
# Result looks different with NDVI, EVI
df$month <- factor(df$month)
df <- tidyr::gather(df, var, value, -c(month, species))

ggplot(data=df, aes(x=factor(month), y=value, fill=species)) + 
  geom_boxplot() + 
  scale_x_discrete(breaks=c(1:12),
                   labels=c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) + 
  labs(x="", y="") + 
  scale_fill_discrete(name="Species", breaks=c("GZ", "PZ"), 
                      labels=c("Equus grevyi", "Equus burchelli")) +
  facet_wrap(~ var, scales= "free") + theme_bw() + 
  theme(strip.background = element_blank(), legend.position="bottom")
ggsave(paste0("Figures/vi_month_sp_", isos[2], ".pdf"), width=8, height=5)
saveRDS(zebra_data_laea, paste0(workdir, "Data/zebra_data_laea_", isos[2], ".rds"))

## ---- move_env_ssf ----

# Read SSF data
ssf_data_exp <- readRDS(paste0(workdir, "Data/ssf_data_exp.rds"))

# Plot boxplot of NDVI for dry vs. rainy season
ssf_data_exp$cat <- NA
ssf_data_exp$cat[month(ssf_data_exp$date) %in% c(4,5,6,7,8,9,10,11)] <- "Dry season"
ssf_data_exp$cat[month(ssf_data_exp$date) %in% c(12,1,2,3)] <- "Wet season"

ggplot(data=data.frame(subset(ssf_data_exp, used==1)), 
       aes(x=id, y=ndvi, fill=factor(cat))) + geom_boxplot() + 
  labs(x="ID", y="NDVI") + theme_bw() + 
  scale_fill_manual(name="Season", labels=c("Dry", "Wet"), values=c("#F8766D", "#00BFC4")) 
ggsave(paste0("Figures/NDVI_", isos[1], "_Dry_vs_Wet.pdf"), width=9, height=7)

df %>% group_by(month) %>%
  summarise(mn_ndvi = mean(ndvi, na.rm=TRUE), 
            sd_ndvi = sd(ndvi, na.rm=TRUE))

# Plot zebra NDVI over time
ggplot(data=data.frame(subset(ssf_data_exp, used==1)), 
       aes(x = date, y = ndvi, colour=id)) + geom_line() + theme_bw() + 
  scale_colour_discrete(name="ID") + labs(x="Date", y="NDVI")
ggsave(paste0("Figures/Zebra_ndvi_", isos[1], ".pdf"), width=10, height=5)

# Plot NDVI against month
ssf_data_exp$month <- lubridate::month(ssf_data_exp$date)
df <- data.frame(ssf_data_exp); rm(ssf_data_exp)
df<- df[c("month", "used", "awei_ns", "brightness", "goal", 
          "LS_ndmi", "LS_ndvi", "ndmi", "ndvi", "slope")]
p1 <- ggplot() + geom_boxplot(data=df, aes(x=factor(month), y=awei_ns, fill=factor(used)),
                              outlier.shape=NA) +
  scale_y_continuous(limits = quantile(df$awei_ns, c(0.1, 0.9), na.rm=TRUE)) + 
  scale_x_discrete(breaks=c(1:12), 
                   labels=c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) + 
  labs(x="Month of the Year", y="awei_ns") + 
  scale_fill_manual(name="", breaks=c(0,1), values=c("#F8766D", "#00BFC4"), 
                    labels=c("Absence", "Presence")) + 
  theme_bw() + theme(legend.position = "bottom")
p2 <- ggplot() + geom_boxplot(data=df, aes(x=factor(month), y=brightness, fill=factor(used)),
                              outlier.shape=NA) +
  scale_y_continuous(limits = quantile(df$brightness, c(0.1, 0.9), na.rm=TRUE)) + 
  scale_x_discrete(breaks=c(1:12), 
                   labels=c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) + 
  labs(x="Month of the Year", y="brightness") + 
  scale_fill_manual(name="", breaks=c(0,1), values=c("#F8766D", "#00BFC4"), 
                    labels=c("Absence", "Presence")) + 
  theme_bw() + theme(legend.position = "bottom")
p3 <- ggplot() + geom_boxplot(data=df, aes(x=factor(month), y=goal, fill=factor(used)),
                              outlier.shape=NA) +
  scale_y_continuous(limits = quantile(df$goal, c(0.1, 0.9), na.rm=TRUE)) + 
  scale_x_discrete(breaks=c(1:12), 
                   labels=c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) + 
  labs(x="Month of the Year", y="Goal") + 
  scale_fill_manual(name="", breaks=c(0,1), values=c("#F8766D", "#00BFC4"), 
                    labels=c("Absence", "Presence")) + 
  theme_bw() + theme(legend.position = "bottom")
p4 <- ggplot() + geom_boxplot(data=df, aes(x=factor(month), y=LS_ndmi, fill=factor(used)),
                              outlier.shape=NA) +
  scale_y_continuous(limits = quantile(df$LS_ndmi, c(0.1, 0.9), na.rm=TRUE)) + 
  scale_x_discrete(breaks=c(1:12), 
                   labels=c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) + 
  labs(x="Month of the Year", y="LS_ndmi") + 
  scale_fill_manual(name="", breaks=c(0,1), values=c("#F8766D", "#00BFC4"), 
                    labels=c("Absence", "Presence")) + 
  theme_bw() + theme(legend.position = "bottom")
p5 <- ggplot() + geom_boxplot(data=df, aes(x=factor(month), y=LS_ndvi, fill=factor(used)),
                              outlier.shape=NA) +
  scale_y_continuous(limits = quantile(df$LS_ndvi, c(0.1, 0.9), na.rm=TRUE)) + 
  scale_x_discrete(breaks=c(1:12), 
                   labels=c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) + 
  labs(x="Month of the Year", y="LS_ndvi") + 
  scale_fill_manual(name="", breaks=c(0,1), values=c("#F8766D", "#00BFC4"), 
                    labels=c("Absence", "Presence")) + 
  theme_bw() + theme(legend.position = "bottom")
p6 <- ggplot() + geom_boxplot(data=df, aes(x=factor(month), y=ndmi, fill=factor(used)),
                              outlier.shape=NA) +
  scale_y_continuous(limits = quantile(df$ndmi, c(0.1, 0.9), na.rm=TRUE)) + 
  scale_x_discrete(breaks=c(1:12), 
                   labels=c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) + 
  labs(x="Month of the Year", y="ndmi") + 
  scale_fill_manual(name="", breaks=c(0,1), values=c("#F8766D", "#00BFC4"), 
                    labels=c("Absence", "Presence")) + 
  theme_bw() + theme(legend.position = "bottom")
p7 <- ggplot() + geom_boxplot(data=df, aes(x=factor(month), y=ndvi, fill=factor(used)),
                              outlier.shape=NA) +
  scale_y_continuous(limits = quantile(df$ndvi, c(0.1, 0.9), na.rm=TRUE)) + 
  scale_x_discrete(breaks=c(1:12), 
                   labels=c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) + 
  labs(x="Month of the Year", y="ndvi") + theme(legend.position = "bottom") + 
  scale_fill_manual(name="", breaks=c(0,1), values=c("#F8766D", "#00BFC4"), 
                    labels=c("Absence", "Presence")) + 
  theme_bw() + theme(legend.position = "bottom")
p8 <- ggplot() + geom_boxplot(data=df, aes(x=factor(month), y=slope, fill=factor(used)),
                              outlier.shape=NA) +
  scale_y_continuous(limits = quantile(df$slope, c(0.1, 0.9), na.rm=TRUE)) + 
  scale_x_discrete(breaks=c(1:12), 
                   labels=c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) + 
  labs(x="Month of the Year", y="slope") + theme(legend.position = "bottom") + 
  scale_fill_manual(name="", breaks=c(0,1), values=c("#F8766D", "#00BFC4"), 
                    labels=c("Absence", "Presence")) + 
  theme_bw() + theme(legend.position = "bottom")
cowplot::plot_grid(p1,p2,p3,p4,p5,p6,p7,p8, ncol=2, labels = "auto")
ggsave(paste0("Figures/vi_month_", isos[1], ".pdf"), width=6, height=12)

## ---- envdata ----

# Read SSF data
ssf_data_exp <- readRDS(paste0(workdir, "Data/ssf_data_exp.rds"))

# Stack env data
env_data <- list(paste0(workdir, "Results/dem_", isos[1], ".tif"), 
                 paste0(workdir, "Results/slope_", isos[1], ".tif"), 
                 paste0(workdir, "Results/aspect_", isos[1], ".tif"))
env_data <- stack(env_data)

# Read goal data
goal <- raster(paste0(workdir, "Results/goal_30m_BWA.tif"))

# Extract env data for ssf_data_exp
env_extract <- raster::extract(env_data, ssf_data_exp)
colnames(env_extract) <- c("dem", "slope", "aspect")
ssf_data_exp@data <- cbind(ssf_data_exp@data, env_extract)

# Extract goal data
ssf_data_exp$goal <- extract(goal, ssf_data_exp)

# Save ssf_data_exp
saveRDS(ssf_data_exp, paste0(workdir, "Data/ssf_data_exp.rds"))

## ---- ssf ----

# Read SSF data
ssf_data_exp <- readRDS(paste0(workdir, "Data/ssf_data_exp.rds"))
(colnames(ssf_data_exp@data))
summary(ssf_data_exp)

# Define variables
vars <- c("dem", "slope", "aspect", "ndvi", "evi", "ndwi_mcf", 
          "ndwi_xu", "ndmi", "awei_ns", "awei_sh", "wi_2015", "brightness", 
          "greenness", "wetness", "NDVI", "EVI", "LS_evi", "LS_ndvi", 
          "LS_ndwi_xu", "LS_ndwi_mcf", "LS_ndmi", "LS_awei_ns", 
          "LS_awei_sh", "LS_wi_2015", "LS_brightness", "LS_greenness", "LS_wetness", "goal")

# Remove NA's
ssf_data_sub <- ssf_data_exp[colnames(ssf_data_exp@data) %in% vars]
ssf_data_sub <- ssf_data_sub[complete.cases(ssf_data_sub@data),]

# Collinearity
(cm <- round(cor(ssf_data_sub@data), 2))
cairo_pdf(paste0(workdir, "Figures/Collinearity_Matrix.pdf"), 
          width=6, height=6, bg="transparent")
ellipse::plotcorr(cm, col=ifelse(abs(cm) > 0.7, "red", "grey"), 
                  main="", mar=c(0.1,0.1,0.1,0.1))
dev.off()

# Obtain cross-combinations of uncollinear variables
cm <- data.frame(cm)
cm_na <- cm
cm_na[cm_na == 1] <- NA
(vars_uncor <- colnames(cm)[unlist(apply(cm_na, FUN=function(x) 
  all(abs(x) < 0.7, na.rm=TRUE), MARGIN=2))])
(vars_cor <- colnames(cm)[!colnames(cm) %in% vars_uncor])
(cor <- cm_na[colnames(cm_na) %in% vars_cor, 
              rownames(cm_na) %in% vars_cor,])

# Obtain different combinations of uncorrelated variables
var_uncor_comb <- list()
for(j in 1:ncol(cor)){
  (var_uncor_comb[[j]] <- c(vars_uncor, vars_cor[j]))
  (vars_uncor_ind <- c(na.omit(colnames(cor)[abs(cor[,j]) < 0.7])))
  (cm_sub <- cm[colnames(cm) %in% c(vars_uncor_ind), 
                rownames(cm) %in% c(vars_uncor_ind)])
  vars_uncor_ind_new <- vars_uncor_ind
  for(k in 1:length(vars_uncor_ind)){
    if(!vars_uncor_ind[[k]] %in% vars_uncor_ind_new){
      next
    }else{
      (vars_uncor_ind_l2 <- c(na.omit(colnames(cm_sub)[abs(cm_sub[,vars_uncor_ind[k]]) < 0.7])))
      (vars_uncor_ind_new <- vars_uncor_ind_new[vars_uncor_ind_new %in% vars_uncor_ind_l2])
      (var_uncor_comb[[j]] <- 
          append(var_uncor_comb[[j]], vars_uncor_ind[k]))
    }
  }
}

# Remove identical combinations
var_uncor_comb <- lapply(var_uncor_comb, 
                         FUN=function(x) x[order(x)])
var_uncor_comb <- lapply(var_uncor_comb, as.factor)
var_uncor_comb <- var_uncor_comb[!duplicated(var_uncor_comb)]
var_uncor_comb <- lapply(var_uncor_comb, as.character)

# Create uncorrelated model formulas
(uncor_fmlas <- lapply(var_uncor_comb, FUN=function(x) 
  as.formula(paste("used ~ cos(rel.angle) + strata(stratum) + ", 
                   paste(x, collapse= "+")))))

# Convert data into list split by ID
ssf_data_list <- lapply(levels(ssf_data_exp$id), FUN=function(x) 
  subset(ssf_data_exp, ssf_data_exp$id == x))

# Run logistic regression model with all non-correlated variables
library(survival)
library(pbs)
m <- lapply(uncor_fmlas, FUN=function(formula) 
  lapply(ssf_data_list, FUN=function(data) 
    clogit(formula=formula, data = data)))

# Calculate AIC for all models and all individuals
m_AIC <- data.frame(t(sapply(m, FUN=function(x) sapply(x, AIC))))

# Create table with Model formulas and AIC
colnames(m_AIC) <- levels(ssf_data_exp$id)
m_AIC$Mean <- apply(m_AIC[,c(1:7)], MARGIN=1, FUN="mean")
m_AIC$Variables <- sapply(var_uncor_comb, FUN=function(x) as.character(paste(x, collapse=" + ")))
m_AIC$Formula <- as.character(uncor_fmlas)
m_AIC$Model <- sapply(1:nrow(m_AIC), FUN=function(x) paste0("M", x))

#' Table with different models
library(xtable)
print(xtable(m_AIC[c("Model", "Variables")], 
             caption=c("Explanatory variables of each model.", "Explanatory variables of each model"), 
             label=paste0("table:summary_models_", isos[1])), 
      type="latex", file=paste0("Tables/summary_models_", isos[1], ".tex"), 
      caption.placement="top", include.rownames=FALSE, 
      bookstab=TRUE, table.placement="H", scalebox="0.85")

#' Table with AIC of the different models
print(xtable(m_AIC[c("Model", levels(ssf_data_exp$id), "Mean")], 
             caption=c("Akaike Information Criterion (AIC) for each 
             model and individual and 
             mean AIC for each model.", "AIC for each 
             model and individual and 
             mean AIC for each model"), 
             label=paste0("table:summary_AIC_", isos[1])), 
      type="latex", file=paste0("Tables/summary_AIC_", isos[1], ".tex"), 
      caption.placement="top", include.rownames=FALSE, 
      bookstab=TRUE, table.placement="H")

#' Drop collinear variables according to its explanatory value
fmla_final <- as.formula(m_AIC$Formula[m_AIC$Mean == min(m_AIC$Mean)])

#' Calculate step selection function
m_final <- list()
p_final <- list()
#' Backwise model selection until only significant variables remain
for(j in 1:length(ssf_data_list)){
  #' Run final model
  m <- clogit(fmla_final, data = ssf_data_list[[j]])
  
  # Get variables  
  (var_final <- c(
    strsplit(as.character(fmla_final), 
             split=" + ", fixed=TRUE)[[3]][-c(1,2,3)]))
  
  #' Get summary
  summary(m)
  
  #' Get p-values
  (pValues <- summary(m)$coefficients[,5])
  
  #' Drop variable with highest p-value from variables
  while(max(pValues, na.rm=TRUE) > 0.05){
    (var_final <- var_final[var_final != names(which.max(pValues))])
    (fmla_sub <- as.formula(paste0("used ~ cos(rel.angle) + strata(stratum) + ", 
                                     paste(var_final, collapse= "+"))))
    m <- clogit(fmla_sub, data = ssf_data_list[[j]])
    (pValues <- summary(m)$coefficients[,5])
  }
  
  # Save model and p-values
  m_final[[j]] <- m; rm(m)
  print(j)
}

# Create table with p-Values of final models
p_final <- sapply(m_final, FUN=function(x) summary(x)$coefficients[,5])
names(p_final) <-  levels(ssf_data_exp$id)
pValues <- data.frame(Variable = unique(unlist(sapply(p_final, FUN=function(x) names(x)))), 
                      t(do.call(rbind, lapply(lapply(p_final, unlist), "[", 
                                            unique(unlist(c(sapply(p_final,names))))))))
print(xtable(pValues, caption=c("p-Values for each component of the final model and each individual.", 
                                "p-Values for each variable and each individual"), 
             label=paste0("table:summary_pvalue_", isos[1]), 
             digits=3, display=c("s", rep("g", length=ncol(pValues)))), 
      type="latex", file=paste0("Tables/summary_pvalue_", isos[1], ".tex"), 
      caption.placement="top", include.rownames=FALSE, 
      bookstab=TRUE, table.placement="H")

# Create table with p-Values and parameter estimates
pValues <- tidyr::gather(pValues, ID, pValue, -Variable)
coef_final <- sapply(m_final, FUN=function(x) summary(x)$coefficients[,1])
names(coef_final) <-  levels(ssf_data_exp$id)
se_final <- sapply(m_final, FUN=function(x) summary(x)$coefficients[,3])
names(se_final) <-  levels(ssf_data_exp$id)
coefs <- data.frame(Variable = unique(unlist(sapply(coef_final, FUN=function(x) names(x)))), 
                    t(do.call(rbind, lapply(lapply(coef_final, unlist), "[", 
                                            unique(unlist(c(sapply(coef_final,names))))))))
print(xtable(coefs, caption=c("Coefficients for each component of the final model and each individual.", 
                              "Coefficients for each variable and each individual"), 
             label=paste0("table:summary_coef_", isos[1]), 
             digits=3, display=c("s", rep("g", length=ncol(coefs)))), 
      type="latex", file=paste0("Tables/summary_coef_", isos[1], ".tex"), 
      caption.placement="top", include.rownames=FALSE, 
      bookstab=TRUE, table.placement="H")
coefs <- tidyr::gather(coefs, ID, Coefficient, -Variable)
se_coefs <- data.frame(Variable = unique(unlist(sapply(se_final, FUN=function(x) names(x)))), 
                       t(do.call(rbind, lapply(lapply(se_final, unlist), "[", 
                                               unique(unlist(c(sapply(se_final,names))))))))
se_coefs <- tidyr::gather(se_coefs, ID, SE, -Variable)

model_sum <- na.omit(cbind(coefs, se_coefs, pValues)[,c(1,2,3,6,9)])
print(xtable(model_sum, caption=c("Model coefficient, SE and p-Value for each component of the 
             final model and each individual.", "Model summary for each variable and each individual"), 
             label=paste0("table:summary_model_", isos[1]), 
             digits=3, display=c("s", rep("g", length=ncol(model_sum)))), 
      type="latex", file=paste0("Tables/summary_model_", isos[1], ".tex"), 
      caption.placement="top", include.rownames=FALSE, 
      bookstab=TRUE, table.placement="H")
