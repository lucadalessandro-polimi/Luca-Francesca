import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from scipy.io import mmread


def load_numeric_csv(path):
    df = pd.read_csv(path)

    # Rimuove colonne completamente vuote
    df = df.dropna(axis=1, how='all')

    # Tiene solo colonne numeriche
    df = df.select_dtypes(include=["number"])

    return df.to_numpy()


def plot_f_on_mesh(folder_path):

    points_path = os.path.join(folder_path, "points.csv")
    cells_path  = os.path.join(folder_path, "cells.csv")
    f_path      = os.path.join(folder_path, "f.mtx")

    # =========================
    # Load points
    # =========================
    points_df = pd.read_csv(points_path)

    # prendiamo SOLO V1 e V2
    x = points_df["V1"].astype(str).str.strip().astype(float).to_numpy()
    y = points_df["V2"].astype(str).str.strip().astype(float).to_numpy()

    print("x range:", x.min(), x.max())
    print("y range:", y.min(), y.max())

    # =========================
    # Load cells
    # =========================
    cells_df = pd.read_csv(cells_path)
    # prendiamo solo V1, V2, V3
    cells = cells_df[["V1", "V2", "V3"]].apply(lambda col: col.astype(str).str.strip()).astype(int).to_numpy()
    # Converti a 0-based indexing
    cells -= 1

    # =========================
    # Load f
    # =========================
    f = mmread(f_path)

    if hasattr(f, "toarray"):
        f = f.toarray()

    f = np.array(f).flatten()

    if len(f) != len(x):
        raise ValueError("f must be node-based and match number of points")

    # =========================
    # Triangulation
    # =========================
    triang = mtri.Triangulation(x, y, cells)

    # =========================
    # Plot
    # =========================
    '''
    width  = x.max() - x.min()
    height = y.max() - y.min()
    ratio  = width / height
    plt.figure(figsize=(8, 8/ratio))

    print(f.max())
    tpc = plt.tripcolor(triang, f, shading='gouraud', cmap='viridis', vmin=0, vmax=72)
    plt.colorbar(tpc, label="PM₁₀ concentration [µg/m³]")

    plt.gca().set_aspect("equal")
    #plt.title("Distribution of f on mesh")
    plt.axis("off")

    output_path = os.path.join(folder_path, "f_distribution.png")
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved in {output_path}")
    '''
        # =========================
    # Plot
    # =========================
    width  = x.max() - x.min()
    height = y.max() - y.min()
    ratio  = width / height

    fig, ax = plt.subplots(figsize=(8, 8/ratio))

    tpc = ax.tripcolor(triang, f, shading='gouraud',
                       cmap='viridis', vmin=0, vmax=72)

    ax.set_aspect("equal")
    ax.axis("off")

    # Colorbar ORIZZONTALE sotto
    cbar = fig.colorbar(
        tpc,
        ax=ax,
        orientation='horizontal',
        fraction=0.05,   # altezza della barra
        pad=0.08         # distanza dall'immagine
    )
    cbar.set_label("PM₁₀ concentration [µg/m³]")

    output_path = os.path.join(folder_path, "f_distribution.png")
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python plot_f.py <mesh_folder>")
        sys.exit(1)

    plot_f_on_mesh(sys.argv[1])