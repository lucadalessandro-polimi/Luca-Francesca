import json
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patheffects as pe


'''
def plot_dcel(filename):
    # Carica il file JSON esportato dalla DCEL
    with open(filename, 'r') as file:
        data = json.load(file)

    # Estrai i nodi
    nodes = {node["id"]: np.array(node["coords"]) for node in data["nodes"]}

    # Crea figura senza assi e sfondo
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_aspect('equal')
    ax.axis('off')  # Nasconde assi e numeri

    # Disegna ogni arco
    for edge in data["edges"]:
        from_id = edge["from"]
        to_id = edge["to"]
        is_sub = edge.get("segment", False)

        p1 = nodes[from_id]
        p2 = nodes[to_id]

        color = "black" if is_sub else "green"
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color, linewidth=1.5)

    # Disegna triangoli (celle)
    halfedge_to_node = {edge["id"]: edge["from"] for edge in data["edges"]}
    for cell in data["cells"]:
        try:
            cell_nodes = [nodes[halfedge_to_node[edge]] for edge in cell["edges"]]
            x_values = [pt[0] for pt in cell_nodes] + [cell_nodes[0][0]]
            y_values = [pt[1] for pt in cell_nodes] + [cell_nodes[0][1]]
            ax.plot(x_values, y_values, 'g-', linewidth=0.5)
        except KeyError:
            print(f"Error in drawing cell with edges {cell['edges']}")

    # Salvataggio immagine pulita
    plt.savefig("Meshes/Delaunay/delaunay_plot.png", dpi=300, bbox_inches='tight', pad_inches=0, transparent=True)
    print("Plot saved as 'Meshes/Delaunay/delaunay_plot.png'")

if __name__ == "__main__":
    plot_dcel("Meshes/Delaunay/delaunay_output.json")
'''

def plot_dcel(filename, out_png="Meshes/Delaunay/delaunay_plot.png", label_edges=True):
    # Carica il file JSON esportato dalla DCEL
    with open(filename, 'r') as file:
        data = json.load(file)

    # Estrai i nodi
    nodes = {node["id"]: np.array(node["coords"]) for node in data["nodes"]}

    # Pre-calcolo scala per l’offset delle etichette (dipende dalle dimensioni del dominio)
    all_pts = np.vstack(list(nodes.values()))
   
    span = float(np.ptp(all_pts, axis=0).max())

    # offset ~1% della dimensione massima (regolabile)
    label_offset = 0.0125 * (span if span > 0 else 1.0)

    # Crea figura senza assi e sfondo
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_aspect('equal')
    ax.axis('off')  # Nasconde assi e numeri

    # Disegna ogni arco + etichetta id
    for edge in data["edges"]:
        from_id = edge["from"]
        to_id = edge["to"]
        is_sub = edge.get("segment", False)

        p1 = nodes[from_id]
        p2 = nodes[to_id]

        color = "black" if is_sub else "green"
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color, linewidth=1.5, zorder=1)

        if label_edges:
            # punto medio e normale per posizionare l’etichetta "sopra" l’arco
            mid = 0.5 * (p1 + p2)
            dirv = p2 - p1
            n = np.array([-dirv[1], dirv[0]], dtype=float)  # normale 2D
            n_norm = np.linalg.norm(n)
            if n_norm > 0:
                n /= n_norm
            text_xy = mid + n * label_offset

            #ax.text(text_xy[0], text_xy[1], str(edge.get("id", "")),ha="center", va="center", fontsize=7, color=color, zorder=2,path_effects=[pe.withStroke(linewidth=2, foreground="white")])

    # Disegna triangoli (celle)
    halfedge_to_node = {edge["id"]: edge["from"] for edge in data["edges"]}
    for cell in data["cells"]:
        try:
            cell_nodes = [nodes[halfedge_to_node[edge]] for edge in cell["edges"]]
            x_values = [pt[0] for pt in cell_nodes] + [cell_nodes[0][0]]
            y_values = [pt[1] for pt in cell_nodes] + [cell_nodes[0][1]]
            ax.plot(x_values, y_values, 'g-', linewidth=0.5, zorder=0)
        except KeyError:
            print(f"Error in drawing cell with edges {cell['edges']}")

    # Salvataggio immagine pulita
    plt.savefig(out_png, dpi=300, bbox_inches='tight', pad_inches=0, transparent=True)
    print(f"Plot salvato come '{out_png}'")

if __name__ == "__main__":
    plot_dcel("Meshes/Delaunay/delaunay_output.json")

