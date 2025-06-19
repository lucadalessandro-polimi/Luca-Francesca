import json
import matplotlib.pyplot as plt
import numpy as np

def plot_dcel(filename):
    # Carica il file JSON esportato dalla DCEL
    with open(filename, 'r') as file:
        data = json.load(file)

    # Estrai i nodi

    nodes = {node["id"]: np.array(node["coords"]) for node in data["nodes"]}

    # Disegna ogni arco con colore diverso se è un subsegment
    plt.figure(figsize=(8, 6))

    for edge in data["edges"]:
        from_id = edge["from"]
        to_id = edge["to"]
        is_sub = edge.get("subsegment", False)

        p1 = nodes[from_id]
        p2 = nodes[to_id]

        color = "black" if is_sub else "green"
        plt.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color, linewidth=1.5)

    # Disegna anche i triangoli (celle)
    halfedge_to_node = {edge["id"]: edge["from"] for edge in data["edges"]}
    for cell in data["cells"]:
        try:
            cell_nodes = [nodes[halfedge_to_node[edge]] for edge in cell["edges"]]
            x_values = [pt[0] for pt in cell_nodes] + [cell_nodes[0][0]]
            y_values = [pt[1] for pt in cell_nodes] + [cell_nodes[0][1]]
            plt.plot(x_values, y_values, 'g-', linewidth=0.5)
        except KeyError:
            print(f"Errore nel disegnare la cella con edges {cell['edges']}")

    plt.title("DCEL Mesh Visualization")
    plt.axis("equal")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.grid(True)
    plt.savefig("Tests_mesh_efficiency/dcel_plot.png", dpi=300)
    print("Plot salvato come 'dcel_plot.png'")

if __name__ == "__main__":
    plot_dcel("Tests_mesh_efficiency/dcel_output.json")

