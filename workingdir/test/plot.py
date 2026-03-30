import json
import matplotlib.pyplot as plt
import numpy as np

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
    #ax.grid(True, linestyle='--', alpha=0.5)

    # Disegna ogni arco
    for edge in data["edges"]:
        from_id = edge["from"]
        to_id = edge["to"]
        is_sub = edge.get("segment", False)

        p1 = nodes[from_id]
        p2 = nodes[to_id]

        if is_sub:
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color="black", linewidth=0.8)  # bordo nero più spesso
        else:
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color="green", linewidth=0.5)  # bordo verde sottile


    # Disegna triangoli (celle)
    halfedge_to_node = {edge["id"]: edge["from"] for edge in data["edges"]}
    for cell in data["cells"]:
        try:
            cell_nodes = [nodes[halfedge_to_node[edge]] for edge in cell["edges"]]
            x_values = [pt[0] for pt in cell_nodes] + [cell_nodes[0][0]]
            y_values = [pt[1] for pt in cell_nodes] + [cell_nodes[0][1]]
            #ax.plot(x_values, y_values, 'g-', linewidth=0.01)
        except KeyError:
            print(f"Error in drawing cell with edges {cell['edges']}")

    #highlight_point = (10315.03958, 43589.37505)
    #ax.plot(highlight_point[0], highlight_point[1], 'ro', markersize=6)
    black_points = np.array([
        [14.0433, 41.3938],
        [14.5049, 41.3851],
        [14.4562, 41.4288],
        [14.0433, 41.3938],
        [13.9788, 41.4635]
    ])
    #ax.scatter(black_points[:, 0], black_points[:, 1], color='black', s=20, label='Special Points')
    all_coords = np.array(list(nodes.values()))
    xmin, ymin = np.min(all_coords, axis=0)
    xmax, ymax = np.max(all_coords, axis=0)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    # Salvataggio immagine pulita
    fig.patch.set_alpha(0.0)
    ax.patch.set_alpha(0.0)

    plt.savefig("delaunay_plot.png", dpi=300, bbox_inches='tight', pad_inches=0, transparent=True)
    print("Plot saved as 'test/delaunay_plot.png'")

if __name__ == "__main__":
    plot_dcel("delaunay_output.json")
