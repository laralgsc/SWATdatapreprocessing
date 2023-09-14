library(sf)
library(tidyverse)
library(terra)


library(dplyr)
library(raster)
library(tidyterra)
library(lubridate)
library(mgsub)

past <- function(...) {
  paste(..., sep = "")
} # Modifica o default de paste para sep = ""

## Variaveis
dir <- "/Users/laralagesanches/Documents/BolsaDTI"; setwd(dir); dir()
folder_from <- "farinha"
folder_to <- "DataSetFarinha";
  subfolders <- c("0DatabaseLara", "1DelineateWatershed", "2CreateHRUs", "3EditInputs", "4Calibration")

DEM_folder <- "SRTMoriginal"
CLIM_folder <- "climate"
basin_file <- "farinha-basin.gpkg"
lulc_file <- "mapbiomas-brazil-collection-71-cerrado-2021.tif"
soils_file <- "pedo_area.shp"
EPSG <- 31983 #Sirgas2000 UTM23S

## Pre processing
file.path(folder_to, subfolders) %>% sapply(dir.create, recursive = T) #Cria minhas pastas de input do modelo
basin <- st_read(basin_file) %>%
  st_transform(EPSG) %>% 
  st_buffer(1000)

######################################## DEM ########################################

DEM <- list.files(file.path(folder_from, DEM_folder), pattern = ".hgt", full.names = T) %>% 
  sprc() %>% 
  merge() %>% 
  project(past("epsg:", EPSG), method = "bilinear", res = 30) %>% 
  crop(basin, mask = T) %>% 
  writeRaster(file.path(folder_to, subfolders[2], "SRTM1arc_proj.tif"), overwrite = T)
#plot(DEM)

#DEM <- rast(file.path(folder_to, subfolders[2], "SRTM1arc_proj.tif"))

######################################## LULC ########################################

# Mapa
LULC <- rast(file.path(folder_from, lulc_file)) %>% 
  project(DEM, method = "near") %>% 
  crop(basin, mask = T) %>% 
  writeRaster(file.path(folder_to, subfolders[3], "LULC.tif"), overwrite = T)
#plot(LULC)
#freq(LULC)

# Lookup table
unique(LULC, na.rm = T) %>% 
  rename(LANDUSE_ID = 1) %>% 
  mutate(SWAT_CODE = NA) %>% 
  write_csv(file.path(folder_to, subfolders[3], "lulc_code.csv"), na = "")

#LULC[is.na(LULC)] <- -1
#LULC[LULC==0] <- NA
#while (freq(LULC, value = NA) != 0) {
#LULC = focal(LULC, matrix(1,3,3), fun = modal, na.rm = T, pad = T, NAonly = T)
#}
#LULC[LULC == -1] <- NA
#writeRaster(LULC, file.path(paste(dir, "LULC_0rm.tif", sep = "")), overwrite = T)

#LULC <- rast(file.path(folder_to, subfolders[3], "LULC.tif"))

######################################## SOILS ########################################

# Mapa
soils <- st_read(soils_file) %>% 
  st_transform(EPSG) %>% 
  st_intersection(basin) %>% 
  mutate(cod_unidad = as.numeric(factor(nom_unidad)))

# Lookup table
soils %>%
  st_drop_geometry() %>% 
  dplyr::select(SOIL_ID = cod_unidad, SNAM = nom_unidad) %>% 
  mutate(SNAM = str_remove_all(SNAM, "[:digit:]")) %>% 
  arrange(SOIL_ID) %>% 
  unique() %>%
  na.omit() %>% 
  write_csv(file.path(folder_to, subfolders[3], "soils_code.csv"))

soils <- soils %>%
  rasterize(DEM, "cod_unidad") %>% 
  writeRaster(file.path(folder_to, subfolders[3], "soils.tif"), overwrite = T)
#plot(soils)
#freq(soils)

#soils <- rast(file.path(folder_to, subfolders[3], "soils.tif"))

####################################### CLIMATE #######################################

dir_clim <- file.path(dir, folder_from, CLIM_folder); dir(dir_clim)
product <- "rad"  #Name of climate variable as written in dir(dir_clim)

raster_list <- list.files(path = file.path(dir_clim, past(product, c("_CFSR", "_CFSv2"))),
                          pattern = '.grb2',
                          full.names = T)
t <- tibble(raster_list, date = parse_date(sub(".*(\\d{6}).*", "\\1", raster_list), "%Y%m")) %>% 
  arrange(date)
#t$raster_list

### Grided meteorological data to csv
## Reprojection
raster_v1 <- rast(raster_list[1]) #plot(raster_v1[[2]]); plot(basin[1], add = T)
raster_v2 <- rast(tail(t$raster_list, 1)) #plot(raster_v2[[1]], add = T)
#st_crs(raster_v1) == st_crs(raster_v2)
#res(raster_v1) == res(raster_v2)

project(raster_v1[[2]], past("epsg:", EPSG), method = "bilinear") %>% plot()

subbasins <- st_read(past(dir, "SWATinputs/basin/srtmclipsubbasins.shp")) %>% 
  st_as_sf() %>% 
  st_transform(st_crs(raster_v1))

#crop_v1 <- crop(raster_v1[[1]], subbasins); #plot(crop_v1)
#mask_v2 <- mask(raster_v2[[1]], subbasins); #plot(mask_v2, add = T)
#x <- project(raster_v2, crop_v1, method = "bilinear", mask = F, align = F); #plot(x[[1]])
#y <- project(mask_v2, crop_v1, method = "bilinear", mask = F, align = F); #plot(y[[1]])

stack_v1 <- rast(filter(t, grepl('CFSR', raster_list))$raster_list) %>% 
  terra::crop(subbasins)
stack_v2 <- rast(filter(t, grepl('CFSv2', raster_list))$raster_list) %>% 
  project(stack_v1, method = "bilinear", mask = F, align = F)
#res(stack_v1) == res(stack_v2)

## Grid point extration
raster <- rast(raster_list[1])
gridpts <- as.points(raster, values = F, na.all = T) %>% 
  st_as_sf()
gridpts$yx <- past(gsub('-|\\.', '', format(yFromCell(raster, 1:ncell(raster)), nsmall = 1)), "S",
                   gsub('-|\\.', '', format(xFromCell(raster, 1:ncell(raster)), nsmall = 1)), "W")

# Apagar pontos inuteis
subbasins <- st_read(past(dir, "SWATinputs/basin/srtmclipsubbasins.shp")) %>% 
  st_as_sf() %>% 
  st_transform(st_crs(gridpts))

nearpts <- st_nearest_feature(subbasins, gridpts) %>% 
  unique()
gridpts_near <- (gridpts[nearpts,])

#plot(raster[[1]]); plot(subbasins$geometry, add = T); plot(gridpts_near, add = T, col = "red")

st_write(gridpts_near, dsn = past(dir_clim, product, "gridpts.gpkg"), layer = "points", driver = "GPKG", append = F)
#dir(dir_clim, pattern = "_gridpts.gpkg")
#gridpts_near <- st_read(past(dir_clim, "T", "_gridpts.gpkg"))

## Climate Raster Stack

#stack <- rast(t$raster_list)
#stack <- rast(list(stack_v1, stack_v2))

assign("element", NULL)
i = 1
for (i in i:length(t$raster_list)) {
  element <- grep("GRIB_ELEMENT=", describe(t$raster_list[i]), value = T) %>% 
    str_remove("    GRIB_ELEMENT=") %>% 
    append(element)
}

assign("ref_time", NULL)
i = 1
for (i in i:length(t$raster_list)) {
  ref_time <- grep("GRIB_REF_TIME=", describe(t$raster_list[i]), value = T) %>% 
    parse_number() %>% 
    append(ref_time, after = 0)
}

ref_time <- as.POSIXct(ref_time, origin = "1970-01-01", tz = "UTC")
#plot(ref_time)

names(stack) <- paste(element, ref_time, sep = "_") %>% 
  mgsub(c(" ", ":00:00"), c("T", ":00:00Z"))

data <- terra::extract(stack, gridpts_near, method = "simple", xy = T, ID = F)

data$y <- past(gsub('-|\\.', '', format(round(data$y, digits = 1), nsmall = 1)),"S")
data$x <- past(gsub('-|\\.', '', format(round(data$x, digits = 1), nsmall = 1)), "W")

data <- unite(data, "yx", y:x, sep = "") %>% 
  data.table::transpose(keep.names = "id", make.names = "yx")
#View(data)

write_csv(data, past(dir_clim, unique(element), "_CFS.csv"))
          
          RAD6h_sum %>% 
            group_by(date = year(date)) %>% 
            summarise(across(2:11, sum)) %>% 
            pivot_longer(-date) %>%
            ggplot(aes(x = date, y = value, color = name)) + 
            geom_line() +
            scale_x_continuous(limits = c(NA,2011))
          
          with_tz(as.POSIXct(284040000, origin = "1970-01-01", tz = "UTC"), Sys.timezone())
          
          as.POSIXct(284022000, origin = "1970-01-01", tz = "UTC")

### Raw data formatting

#Regarding processing, we simply took every hour's data and
#added it to get total rainfall,
#max/min temperature,
#mean RH, windspeed, and
#sum for SR.
          
# Funcao para salvar um arquivo txt para cada estacao no formato requerido pelo SWAT,
# para funcionar a tabela de entrada deve ter a primeira coluna como data de cada registro e as demais colunas os registros de cada estacao
# Definicao de cada argumento:
# data = tabela com data e dados das estacoes
# dir = diretorio onde salvar
# variable.initial = texto padrao do SWAT para a variavel em questao

wrt_txt <- function(data, dir, variable.initial) {
  i=2
  for (i in i:length(data)) {
    c(str_remove_all(pull(data[1,1]), "-"),
      format(round(pull(data[i]), digits = 3), nsmall = 3, trim = T)) %>% 
      as_tibble() %>% 
      write_delim(paste(dir, variable.initial, "-", names(data[i]), ".txt", sep = ""), delim = ",", col_names = F)
  }
  }

products <- dir(dir_clim, patter = ".csv") %>% 
  str_remove("_CFS.csv")

stations_file = c("rh", "solar", "tmp", "wind")
data_files = c("r-", "s-", "t-", "w-")

## Radiation
# (sum of 6h averages) * (seconds in a 6h period)/1000000
{
product <- products[1]
data <- read_csv(past(dir_clim, product, "_CFS.csv"))

RAD_day <- data %>% 
  mutate(date = as.Date(str_remove(id, past(product, "_"))), .after = 1) %>% 
  group_by(date) %>% 
  summarise(across(where(is.double), sum)) %>% 
  mutate_if(is.numeric, ~.x * 21600/1000000)

wrt_txt(RAD_day, "/Users/laralanches/Documents/BolsaDTI/DadosEspaciais/SWATinputs/climate/", "s")

#gridpts_near <- st_read(paste(dir_clim, "T", "_gridpts.gpkg", sep = ""))
tibble(ID = seq(1:length(select_if(RAD_day, is.numeric))),
       NAME = past(data_files[j], names(select_if(RAD_day, is.numeric))),
       LAT = format(round(st_coordinates(gridpts_near)[,2], digits = 3)),
       LONG = format(round(st_coordinates(gridpts_near)[,1], digits = 3)),
       ELEVATION = NA) %>% 
  write_delim(paste("/Users/laralanches/Documents/BolsaDTI/DadosEspaciais/SWATinputs/climate/", stations_file[j], ".txt", sep = ""), delim = ",", col_names = T)
}

## Wind speed from u and v components
{
product <- products[5]
data <- read_csv(past(dir_clim, product, "_CFS.csv"))

u_comp <- filter(data, grepl("UGRD", id))
v_comp <- filter(data, grepl("VGRD", id))

ws <- sqrt(as.matrix(u_comp[2:7])^2 + as.matrix(v_comp[2:7])^2)
rownames(ws) <- gsub("UGRD_", "", u_comp$id)

# Wind speed by date at 1.7m
WND_day <- as_tibble(ws, rownames = "id") %>% 
  mutate(date = as.Date(id)) %>% 
  group_by(date) %>% 
  summarise(across(2:7, mean)) %>% 
  mutate(across(2:7, ~ .x * (1.7/10)^0.2))

wrt_txt(WND_day, "/Users/laralanches/Documents/BolsaDTI/DadosEspaciais/SWATinputs/climate/", "w")

#gridpts_near <- st_read(past(dir_clim, "RH", "_gridpts.gpkg"))
tibble(ID = seq(1:length(select_if(ws, is.numeric))),
       NAME = past(data_files[j], names(select_if(ws, is.numeric))),
       LAT = format(round(st_coordinates(gridpts_near)[,2], digits = 3), nsmall = 3),
       LONG = format(round(st_coordinates(gridpts_near)[,1], digits = 3), nsmall = 3),
       ELEVATION = NA) %>% 
  write_delim(past("/Users/laralanches/Documents/BolsaDTI/DadosEspaciais/SWATinputs/climate/", stations_file[j], ".txt"), delim = ",", col_names = T)
}

## Relative humidity
{
product <- products[2]
data <- read_csv(past(dir_clim, product, "_CFS.csv"))

RH_day <- data %>% 
  mutate(date = as.Date(str_remove(id, past(product, "_"))), .after = 1) %>% 
  group_by(date) %>% 
  summarise(across(where(is.double), mean)) %>% 
  mutate_if(is.numeric, ~.x * 0.01)

wrt_txt(RH_day, "/Users/laralanches/Documents/BolsaDTI/DadosEspaciais/SWATinputs/climate/", "r")
}

## Temperature
{
#TMAX
product <- products[3]
data <- read_csv(past(dir_clim, product, "_CFS.csv"))

Tmax_day <- data %>% 
  mutate(date = as.Date(str_remove(id, past(product, "_"))), .after = 1) %>% 
  group_by(date) %>% 
  summarise(across(where(is.double), max))

#TMIN
product <- products[4]
data <- read_csv(past(dir_clim, product, "_CFS.csv"))

Tmin_day <- data %>% 
  mutate(date = as.Date(str_remove(id, past(product, "_"))), .after = 1) %>% 
  group_by(date) %>% 
  summarise(across(where(is.double), min))

names(Tmax_day) == names(Tmin_day)

i = 2
for (i in i:length(Tmax_day)) {
ma <-  c(str_remove_all(pull(Tmax_day[1,1]), "-"),
  format(round(pull(Tmax_day[i]), digits = 3), nsmall = 3, trim = T))
mi <-  c("",
  format(round(pull(Tmin_day[i]), digits = 3), nsmall = 3, trim = T))

tibble(ma, mi) %>% 
  write_delim(past("/Users/laralanches/Documents/BolsaDTI/DadosEspaciais/SWATinputs/climate/", "t", "-", names(Tmax_day[i]), ".txt"), delim = ",", col_names = F)
}

}

## Precipitation

txt <- list.files(path = "/Users/laralanches/Documents/BolsaDTI/DadosEspaciais/Original/Climate/PCP_ANA",
           pattern = "([0-9]).txt", full.names = T)

stations <- tibble(.rows = 16098)
i=1
for (i in i:length(txt)){
  pcp <- tibble(read.table(txt[i]))
  data <- pcp[[4]] %>% 
    str_replace_all(",", ".") %>% 
    tibble()
  colnames(data) <- parse_number(txt[i])
  stations <- add_column(stations, data)
}

#past(pcp[[1,3]], "0", pcp[[1,2]], "0", pcp[[1,1]])

stations <- add_column(stations, date = seq(as.Date("1979-01-01"), as.Date("2023-01-27"), by = 1), .before = 1) %>% 
  mutate(across(where(is.character), as.numeric))

stations %>% 
  reshape2::melt(id.vars = 'date', variable.name = 'station') %>% 
  ggplot(aes(date, y = value, color = station)) +
  geom_line() +
  scale_x_date(limits = c(as.Date("1980-01-01"),NA)) +
  facet_wrap(vars(station), nrow = 4) +
  theme(legend.position="none")

empty <- which(colSums(stations[-1]) == -16098) + 1

stations <- stations %>% 
  select(-empty)

stations %>% 
  reshape2::melt(id.vars = 'date', variable.name = 'station') %>% 
  ggplot(aes(date, y = value, color = station)) +
  geom_line() +
  scale_x_date(limits = c(as.Date("1980-01-01"),NA)) +
  facet_wrap(vars(station), nrow = 4) +
  theme(legend.position="none")


stations %>% 
  reshape2::melt(id.vars = 'date', variable.name = 'station') %>% 
  ggplot(aes(date, y = value, color = station)) +
  geom_line() +
  scale_x_date(limits = c(as.Date("1980-01-01"),as.Date("2005-01-01"))) +
  facet_wrap(vars(station), nrow = 4) +
  theme(legend.position="none")


good_stations <- names(stations[-c(1,16)])


c("1643007", "1643038", "1743000", "1743002", "1743015",
  "1743016", "1743018", "1843002", "1843003", "1843011")

stations %>% 
  select(c(1,2,5,7,9,10)) %>% 
  reshape2::melt(id.vars = 'date', variable.name = 'station') %>% 
  ggplot(aes(date, y = value, color = station)) +
  geom_line() +
  scale_x_date(limits = c(as.Date("1980-01-01"), as.Date("1981-01-01")))


                  group_by(ws, year(date)) %>% 
                    ggplot(aes(x = date)) +
                    geom_line()
                  
                  tb <- tibble_row(as.numeric(str_remove_all(ws$date[1], "-")))
                  names(tb) <- names(ws[2])
                  add_row(tb, ws[2])
                  
                  write_csv(ws[2], paste(dir_clim, product, names(ws[2]), ".csv", sep = ""))

                  data <- read_table(unz("/Users/laralanches/Documents/southAmerica.zip", "55081_2022-02-01-02-54-17/s-170-434.txt"))
                  dir("/Users/laralanches/Documents/")
                  
                  files <- c("55081_2022-02-01-02-54-17/s-158-481.txt", "55081_2022-02-01-02-54-17/s-170-434.txt")
                  data <- read_table(unzip("/Users/laralanches/Documents/southAmerica.zip", files[2], exdir = "/Users/laralanches/Documents"))

######################################## WGEN ########################################

dir <- "/Users/laralanches/Documents/WGEN/"
dir(dir)

# Adicionar Text, data inicial (01/01/1979) e data final (10/31/2022)

files <- list.files(path = dir, full.names = T)

i=1
for (i in 1:length(files)) {
  read_delim(files[i], delim = ",", col_types = "c") %>% 
    add_row(tibble('19790101' = c("Text", "01/01/1979", "10/31/2022")), .before = 1) %>% 
    write_delim(files[i], delim = ",", quote = "none", col_names = F)
}

# Criar arquivos vazios de componentes climaticos das estacoes

variables <- c("pcp","tmp","slr","wnd","dwp", "hhr")

stations <- list.files(path = dir, pattern = '_', full.names = F) %>% 
  str_extract(".*_") %>% 
  unique()

i=1
for (i in 1:length(stations)) {
  j=1
  for(j in 1:length(variables)) {
    write_delim(data_frame(NA), paste(dir, stations[i], variables[j], ".txt", sep = ""), delim = ",", append = T, na = "", eol = "")
  }
}
