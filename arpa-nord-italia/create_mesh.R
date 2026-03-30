# ------------------- MESH ----------------------------------------------------- 
library(sf)
centraline = st_read("centraline_nord_italia_PM10_aggregate.gpkg")
nord_italia = st_read("nord_italia.gpkg")

centraline = st_transform(centraline, crs=4326)
nord_italia = st_transform(nord_italia, crs=4326)

plot(st_geometry(nord_italia))
plot(st_union(st_geometry(nord_italia)))
plot(st_geometry(centraline),add = T)

library(rmapshaper)

bd_line = st_union(st_geometry(nord_italia)) 
bd_line = st_cast(bd_line, "POLYGON")

perimeters <- st_length(st_cast(bd_line, "MULTILINESTRING"))
idx <- which.max(perimeters)
bd_line <- bd_line[idx, ]

#help("ms_simplify")
bd_simp = rmapshaper::ms_simplify(bd_line, keep=0.0025, keep_shapes = TRUE)
bd_simp
plot(bd_simp)

head(st_coordinates(bd_simp))

bd_nodes = st_coordinates(bd_simp)[,1:2]
# droppo ultimo nodo che è uguale al primo
bd_nodes = bd_nodes[-1,]

# check antiorario
plot(bd_nodes, pch=16, col="red")
points(bd_nodes[1,1], bd_nodes[1,2], pch=16, col="green2")
points(bd_nodes[2,1], bd_nodes[2,2], pch=16, col="blue2")
points(bd_nodes[3,1], bd_nodes[3,2], pch=16, col="black")
points(bd_nodes[nrow(bd_nodes),1], bd_nodes[nrow(bd_nodes),2], pch=16, col="black")

plot(bd_simp)
plot(st_geometry(centraline), add=T)

is_inside = st_within(st_geometry(centraline), bd_simp, sparse = F)
sum(is_inside)

plot(bd_simp)
plot(st_geometry(centraline)[is_inside==FALSE, ], add=T, col = c("red","blue","green"), pch=16)

centraline = centraline[is_inside,]

library(RTriangle)
bd_edges = cbind( 1:(nrow(bd_nodes)-1), 2:nrow(bd_nodes))
bd_edges = rbind(bd_edges, c(nrow(bd_nodes), 1))
pslg = pslg(P=bd_nodes, S=bd_edges)
plot(pslg)
mesh = RTriangle::triangulate(p = pslg, a = 0.025/2, q=20) # diminuire il coeff a -> "max area" per mesh più fini
# 1040
plot(mesh, pch=".")
dim(mesh$P)


bd_nodes = st_coordinates(bd_simp)[,1:2]
bd_nodes = bd_nodes[-1,]  # togli ultimo = primo

write.csv(bd_nodes, "mesh_1/border_nodes.csv", row.names = FALSE)



nodes = mesh$P
elements = mesh$T
boundary = mesh$PB
mesh = fdaPDE::create.mesh.2D(nodes = nodes, triangles = elements)

if(!dir.exists("mesh_1/")) dir.create("mesh_1/")
write.csv(format(mesh$nodes, digits=16), file = "mesh_1/points.csv")
write.csv(format(mesh$triangles, digits=16), file = "mesh_1/cells.csv")
write.csv(format(mesh$nodesmarkers, digits=16), file = "mesh_1/neigh.csv")


plot(mesh, pch=".")
points(st_coordinates(centraline), col="red3", pch=16)

if(!dir.exists("data/")) dir.create("data/")
write.csv(format(st_coordinates(centraline), digits=16), file = "data/locs.csv")

mesh = RTriangle::triangulate(p = pslg, a = 0.025/3, q=20)

nodes = mesh$P
elements = mesh$T
boundary = mesh$PB
mesh = fdaPDE::create.mesh.2D(nodes = nodes, triangles = elements)
plot(mesh, pch=".")
if(!dir.exists("mesh_2/")) dir.create("mesh_2/")
write.csv(format(mesh$nodes, digits=16), file = "mesh_2/points.csv")
write.csv(format(mesh$triangles, digits=16), file = "mesh_2/cells.csv")
write.csv(format(mesh$nodesmarkers, digits=16), file = "mesh_2/neigh.csv")
