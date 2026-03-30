library(sf)
rm(list=ls())
south_italy <- st_read(dsn = "south-italy-no-sic.shx") # regions
plot(st_geometry(south_italy))

south_regions <- st_read(dsn = "south-italy-regions-no-sic.shx") # regions
plot(st_geometry(south_regions))


fire_data <- st_read(dsn = "../fire-data.shx")
points(st_coordinates(fire_data)[,1:2], pch=16, col="red")

fire_data <- st_filter(fire_data, south_italy)
points(st_coordinates(fire_data)[,1:2], pch=16, col="blue")
st_write(fire_data, dsn="fire-data.shp", delete_dsn = TRUE)

rm(list=ls())
graphics.off()
fire_data <- st_read("fire-data.shx")
south_italy <- st_read(dsn = "south-italy-no-sic.shx") 
south_regions <- st_read(dsn = "south-italy-regions-no-sic.shx") 

plot(st_geometry(south_regions))
points(st_coordinates(south_italy)[,1:2], col="red", type="l")
points(st_coordinates(fire_data)[,1:2], col="red", pch=16)


rm(list=ls())
output_path <- c("1e2/", "1e3/", "1e4/")
coeff <- c(1e2, 1e3, 1e4)

for(i in 1:length(output_path)){
  if(!dir.exists(output_path[i])) dir.create(output_path[i])
  
  fire_data <- st_read(dsn="fire-data.shx")
  geom_data <- st_geometry(fire_data) * coeff[i]
  st_geometry(fire_data) <- geom_data
  st_write(fire_data, dsn = paste0(output_path[i], "fire-data.shp"), 
           delete_dsn = TRUE)
  
  south_italy <- st_read(dsn = "south-italy-no-sic.shx") # regions
  geom <- st_geometry(south_italy) * coeff[i]
  st_geometry(south_italy) <- geom
  st_write(south_italy, dsn = paste0(output_path[i], "south-italy-no-sic.shp"), 
           delete_dsn = TRUE)
  
  south_regions <- st_read(dsn = "south-italy-regions-no-sic.shp") # regions
  geom <- st_geometry(south_regions) * coeff[i]
  st_geometry(south_regions) <- geom
  st_write(south_regions, dsn = paste0(output_path[i], "south-italy-regions-no-sic.shp"), 
           delete_dsn = TRUE)
}

# prova
i = 3
fire_data <- st_read(dsn=paste0(output_path[i], "fire-data.shp"))
south_italy <- st_read(dsn =paste0(output_path[i],"south-italy-no-sic.shp")) 
south_regions <- st_read(dsn = paste0(output_path[i],"south-italy-regions-no-sic.shp"))

plot(st_geometry(south_regions))
points(st_coordinates(south_italy)[,1:2], col="red", type="l")
points(st_coordinates(fire_data)[,1:2], col="red", pch=16)
