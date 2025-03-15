import json
import networkx as nx
import matplotlib.pyplot as plt

def plot_dcel(filename):
    # Carica il file JSON esportato dalla DCEL
    with open(filename, 'r') as file:
        data = json.load(file)

    # Estrarre i nodi
    nodes = {node["id"]: node["coords"] for node in data["nodes"]}
    boundary_nodes = {node["id"] for node in data["nodes"] if node["boundary"]}

    # Estrarre gli half-edges e creare una mappa halfedge_id -> nodo di partenza
    halfedge_to_node = {}
    edges = []
    edge_labels = {}

    for edge in data["edges"]:
        from_id = edge["from"]
        to_id = edge["to"]
        halfedge_id = edge["id"]
        twin_id = edge["twin"]

        edges.append((from_id, to_id))
        edge_labels[(from_id, to_id)] = f"{halfedge_id}/{twin_id}"  # Mostra id/twin_id
        halfedge_to_node[halfedge_id] = from_id  # Mappa: halfedge_id → nodo di partenza

    # Creare il grafo con NetworkX per una migliore visualizzazione
    G = nx.Graph()
    for node_id, (x, y) in nodes.items():
        G.add_node(node_id, pos=(x, y))

    for from_id, to_id in edges:
        G.add_edge(from_id, to_id)

    # Creare il plot
    plt.figure(figsize=(8, 6))

    # Disegna gli archi
    pos = {node_id: (x, y) for node_id, (x, y) in nodes.items()}
    nx.draw(G, pos, with_labels=True, node_size=300, node_color='black', edge_color='gray', font_color='red')

    # Disegna gli ID degli half-edges sugli archi
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=10, font_color='blue')

    # Disegna i nodi di bordo in blu
    for node_id in boundary_nodes:
        x, y = nodes[node_id]
        plt.scatter(x, y, color='blue', s=100, edgecolors='black', linewidths=1.5)

    print(data['cells'])
    # Disegna le celle
    for cell in data["cells"]:
        cell_nodes = [nodes[halfedge_to_node[edge]] for edge in cell["edges"]]
        x_values = [coord[0] for coord in cell_nodes] + [cell_nodes[0][0]]
        y_values = [coord[1] for coord in cell_nodes] + [cell_nodes[0][1]]
        plt.plot(x_values, y_values, 'g--', linewidth=1)

    # Mostrare la mesh
    plt.title("DCEL Mesh Visualization")
    plt.axis('equal')
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.grid(True)
    plt.savefig("dcel_plot.png", dpi=300)  # Salva l'immagine
    print("Plot salvato come 'dcel_plot.png'")



if __name__ == "__main__":
    plot_dcel("dcel_output.json")

