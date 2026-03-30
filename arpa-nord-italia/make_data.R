library(dplyr)
library(ggplot2)

dati <- read.csv("misurazioni_nord_italia_PM10_aggregate.csv", header = TRUE, sep = ";")

summary(dati)

dati$AggType = as.factor(dati$AggType)
dati$values_numeric = as.numeric(dati$values_numeric)
levels(dati$AggType)
head(dati)

library(sf)

centraline = st_read("centraline_nord_italia_PM10_aggregate.gpkg")
nord_italia = st_read("nord_italia.gpkg")

# impongo lon-lat
centraline = st_transform(centraline, crs=4326)
nord_italia = st_transform(nord_italia, crs = 4326)

#plot(nord_italia)
#plot(centraline)

# integraimo i due dataset tramite il codice nat

length(unique(centraline$air_quality_station_nat_code)) # 266
length(unique(dati$air_quality_station_nat_code))       # 266

locs = st_coordinates(centraline)
range(locs)

centraline.df = cbind(as.integer(centraline$air_quality_station_nat_code), locs[,1], locs[,2])
dim(centraline.df)
head(centraline.df)
colnames(centraline.df) <- c("air_quality_station_nat_code", "x", "y")

centraline.df = as.data.frame(centraline.df) 
centraline.df$air_quality_station_nat_code = as.integer(centraline.df$air_quality_station_nat_code)

# --------------

range(dati$Date)
anni = c("2018", "2019", "2020", "2021", "2022")
dati_anno = list()

dati <- dati %>% mutate(Date = as.Date(Date))
for(i in 1:length(anni))
{
  dati_anno[[i]] = dati %>% filter(format(Date, "%Y") == anni[i])
}


december_2018 = dati_anno[[1]] %>% filter(format(Date, "%m") == "12")
dim(december_2018)

# remove NAs 
december_2018 = december_2018[complete.cases(december_2018), ]
dim(december_2018)

december_2018$day = format( december_2018$Date, "%d" )

# merge !
merged_data <- december_2018 %>%
  left_join(centraline.df, by = "air_quality_station_nat_code")

names(merged_data)

merged_data = merged_data %>% arrange(air_quality_station_nat_code, Date)

write.csv(cbind(merged_data$x, merged_data$y), file="space_locs.csv")
write.csv(merged_data$day, file="time_locs.csv")

PM10 = data.frame( y = merged_data$values_numeric )
write.csv(PM10, file="PM10_dicembre_2018.csv")
# nb in questo modo "y" è la risposta

data_days <- list()
days <- unique(merged_data$day)
days <- sort(days)
for(i in 1:31){
  data_days[[i]] <- merged_data %>% filter(format(Date, "%d") == days[i])
  print(paste0("Giorno ", i, "\tnum data ", nrow( data_days[[i]])))
}

for(i in 1:31){
  write.csv(cbind(data_days[[i]]$x, data_days[[i]]$y), 
          file=paste0("data/2018.12.",i,"_space_locs.csv"))
  write.csv(data.frame(y=data_days[[i]]$values_numeric), 
          file=paste0("data/2018.12.",i,"_PM10.csv"))
}

# uniform mesh  
library(fdaPDE)
nodes <- as.matrix( read.csv("mesh_1/points.csv")[,2:3] )
triangles <- as.matrix( read.csv("mesh_1/cells.csv")[,2:4] )

mesh <- create.mesh.2D(nodes=nodes, triangles = triangles)
plot(mesh, pch=".")
points(data_days[[11]]$x, data_days[[11]]$y, pch=16, col="red")


