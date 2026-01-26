import pandas as pd
import matplotlib.pyplot as plt

# carica il file
df = pd.read_csv("Meshes/Delaunay/sequenze/path_0_to_47.csv")   # cambia 123 con l'id del nodo v0

x = df["x"].to_numpy()
y = df["y"].to_numpy()

plt.figure()
plt.axis([0, 1, 0, 1])  # [xmin, xmax, ymin, ymax]  
plt.plot(x, y, marker="o", linewidth=1)     # linea + marker in ordine
plt.scatter(x[0], y[0], s=80)               # evidenzia start
plt.title("Geodesic visit order (Dijkstra pop order)") 
plt.savefig("Meshes/Delaunay/sequenze/path_0_to_47.png", dpi=200)
plt.show()
