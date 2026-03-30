import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from scipy.io import mmread


def load_mesh(folder_path):
    # =========================
    # Load points
    # =========================
    points_df = pd.read_csv(os.path.join(folder_path, "points.csv"))
    x = points_df["V1"].astype(str).str.strip().astype(float).to_numpy()
    y = points_df["V2"].astype(str).str.strip().astype(float).to_numpy()

    # =========================
    # Load cells
    # =========================
    cells_df = pd.read_csv(os.path.join(folder_path, "cells.csv"))
    cells = (
        cells_df[["V1", "V2", "V3"]]
        .apply(lambda col: col.astype(str).str.strip())
        .astype(int)
        .to_numpy()
    )
    cells -= 1  # 0-based indexing

    # =========================
    # Load f
    # =========================
    f = mmread(os.path.join(folder_path, "f.mtx"))
    if hasattr(f, "toarray"):
        f = f.toarray()
    f = np.array(f).flatten()

    triang = mtri.Triangulation(x, y, cells)

    return x, y, triang, f


def plot_three_meshes(folder_list, titles, output_path):

    # Carichiamo tutto
    meshes = [load_mesh(folder) for folder in folder_list]

    # Calcoliamo aspect ratio dalla prima mesh
    x0, y0, _, _ = meshes[0]
    width  = x0.max() - x0.min()
    height = y0.max() - y0.min()
    ratio  = width / height

    # Figura con 3 righe, 1 colonna
    fig, axes = plt.subplots(
        nrows=3, ncols=1,
        figsize=(8, 3 * 8 / ratio),
        constrained_layout=True
    )

    vmin = 0
    vmax = 72



    for ax, (x, y, triang, f) in zip(axes, meshes):
          tpc = ax.tripcolor(
          triang, f,
          shading='gouraud',
          cmap='viridis',
          vmin=vmin,
          vmax=vmax
          )

          ax.set_aspect("equal")
          ax.axis("off")
    

    # Colorbar unica
    cbar = fig.colorbar(
        tpc,
        ax=axes,
        orientation='vertical',
        fraction=0.03,
        pad=0.02
    )
    cbar.set_label("PM₁₀ concentration [µg/m³]")

    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved in {output_path}")



if __name__ == "__main__":

    folders = [
        "mesh_4",
        "mesh_5",
        "mesh_6"
    ]

    titles = [
        "Uniform mesh",
        "Nodal density-based adapted mesh",
        "Residual-based adapted mesh"
    ]

    plot_three_meshes(
        folders,
        titles,
        "1040_comparison.png"
    )