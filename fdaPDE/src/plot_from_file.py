import json
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

def circumcenter(A, B, C):
    """Calcola il circocentro del triangolo con vertici A, B, C"""
    D = 2 * (A[0] * (B[1] - C[1]) + B[0] * (C[1] - A[1]) + C[0] * (A[1] - B[1]))
    Ux = ((A[0]**2 + A[1]**2) * (B[1] - C[1]) +
          (B[0]**2 + B[1]**2) * (C[1] - A[1]) +
          (C[0]**2 + C[1]**2) * (A[1] - B[1])) / D
    Uy = ((A[0]**2 + A[1]**2) * (C[0] - B[0]) +
          (B[0]**2 + B[1]**2) * (A[0] - C[0]) +
          (C[0]**2 + C[1]**2) * (B[0] - A[0])) / D
    return np.array([Ux, Uy])

def circumradius(A, B, C, U):
    """Calcola il raggio della circonferenza circoscritta"""
    return np.linalg.norm(A - U)

def plot_dcel(filename):
    # Carica il file JSON esportato dalla DCEL
    with open(filename, 'r') as file:
        data = json.load(file)
    print(data['cells'])
    # Estrarre i nodi
    nodes = {node["id"]: np.array(node["coords"]) for node in data["nodes"]}
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
    for node_id, coords in nodes.items():
        G.add_node(node_id, pos=(coords[0], coords[1]))

    for from_id, to_id in edges:
        G.add_edge(from_id, to_id)

    # Creare il plot
    plt.figure(figsize=(8, 6))

    # Disegna gli archi
    pos = {node_id: (coords[0], coords[1]) for node_id, coords in nodes.items()}
    #nx.draw(G, pos, with_labels=True, node_size=100, node_color='black', edge_color='gray', font_color='white', font_size=8)


    # Disegna gli ID degli half-edges sugli archi
    #nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=8, font_color='red')

    # Disegna i nodi di bordo in blu
    for node_id in boundary_nodes:
        x, y = nodes[node_id]
        plt.scatter(x, y, color='blue', s=10, edgecolors='black', linewidths=1.5)

    # Disegna le celle (triangoli)
    for cell in data["cells"]:
        cell_nodes = [nodes[halfedge_to_node[edge]] for edge in cell["edges"]]
        x_values = [coord[0] for coord in cell_nodes] + [cell_nodes[0][0]]
        y_values = [coord[1] for coord in cell_nodes] + [cell_nodes[0][1]]
        plt.plot(x_values, y_values, 'g-', linewidth=1)

        # Calcolare e disegnare la circonferenza circoscritta
        '''if len(cell_nodes) == 3:  # Assicuriamoci che sia un triangolo
            A, B, C = cell_nodes
            U = circumcenter(A, B, C)
            R = circumradius(A, B, C, U)

            circle = plt.Circle(U, R, color='blue', fill=False, linestyle="dotted", linewidth=1.2)
            plt.gca().add_patch(circle)
            #plt.scatter(*U, color='purple', s=50, label="Circocentro")'''
    # Mostrare la mesh
    plt.title("")
    plt.axis('equal')
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.grid(True)
    plt.savefig("dcel_plot.png", dpi=300)  # Salva l'immagine
    print("Plot salvato come 'dcel_plot.png'")

if __name__ == "__main__":
    plot_dcel("dcel_output.json")

