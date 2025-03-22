import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle


def plot_mesh(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Separare nodi e celle
    nodes_section = lines[1:lines.index('\n')]  # Prende tutte le righe fino alla linea vuota
    cells_section = lines[lines.index('\n') + 2:]  # Salta la linea vuota e 'Cells:\n'

    nodes = []
    cells = []

    # Parsing dei nodi
    for line in nodes_section:
        parts = line.strip().split()
        if len(parts) == 4:  # id, x, y, boundary_marker
            nodes.append((int(parts[0]), float(parts[1]), float(parts[2]), int(parts[3])))

    # Parsing delle celle
    for line in cells_section:
        parts = line.strip().split()
        if len(parts) == 4:  # cell_id, node1, node2, node3
            cells.append((int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3])))

    # Plotting
    plt.figure(figsize=(10, 10))

    # Disegna i triangoli
    for cell_id, n1, n2, n3 in cells:
        x = [nodes[n1][1], nodes[n2][1], nodes[n3][1], nodes[n1][1]]
        y = [nodes[n1][2], nodes[n2][2], nodes[n3][2], nodes[n1][2]]
        plt.plot(x, y, color='green', lw=0.8)

        # Disegna l'ID del triangolo
        #centroid_x = np.mean([nodes[n1][1], nodes[n2][1], nodes[n3][1]])
        #centroid_y = np.mean([nodes[n1][2], nodes[n2][2], nodes[n3][2]])
        #plt.text(centroid_x, centroid_y, str(cell_id), fontsize=8, color='orange')

        # Calcola e disegna la circonferenza circoscritta
        A = np.array([nodes[n1][1], nodes[n1][2]])
        B = np.array([nodes[n2][1], nodes[n2][2]])
        C = np.array([nodes[n3][1], nodes[n3][2]])

        # Calcola il centro e il raggio della circonferenza circoscritta
        D = 2 * (A[0] * (B[1] - C[1]) + B[0] * (C[1] - A[1]) + C[0] * (A[1] - B[1]))
        if D != 0:
            Ux = ((np.linalg.norm(A)**2) * (B[1] - C[1]) + (np.linalg.norm(B)**2) * (C[1] - A[1]) + (np.linalg.norm(C)**2) * (A[1] - B[1])) / D
            Uy = ((np.linalg.norm(A)**2) * (C[0] - B[0]) + (np.linalg.norm(B)**2) * (A[0] - C[0]) + (np.linalg.norm(C)**2) * (B[0] - A[0])) / D
            center = np.array([Ux, Uy])
            radius = np.linalg.norm(A - center)

            circle = Circle(center, radius, color='blue', fill=False, linestyle='dotted', lw=0.8)
            #plt.gca().add_patch(circle)

    # Disegna i nodi
    for node_id, x, y, marker in nodes:
        plt.scatter(x, y, color='red' if marker == 1 else 'blue', s=20)
        #plt.text(x, y, str(node_id), fontsize=8, ha='right')

    plt.title('Triangulation Plot')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('equal')
    plt.grid(True)

    # Salva il grafico come immagine
    plt.savefig('mesh_plot.png')

if __name__ == "__main__":
    plot_mesh("mesh_output.txt")