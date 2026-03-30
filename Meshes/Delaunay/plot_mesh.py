import json
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patheffects as pe
from matplotlib.patches import Circle


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


def _load_xy_txt(filepath):
    """Carica un file 2 colonne (x y), separatore spazio/tab o virgola. Ritorna array (n,2) o None."""

    # individua il separatore dalla prima riga non vuota/non commentata
    delimiter = None
    with open(filepath, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            if "," in s:
                delimiter = ","
            break

    try:
        arr = np.loadtxt(filepath, delimiter=delimiter, comments="#")
    except Exception:
        arr = np.genfromtxt(filepath, delimiter=delimiter, comments="#")

    if arr is None or np.size(arr) == 0:
        return None
    arr = np.atleast_2d(arr)
    if arr.shape[1] < 2:
        return None
    return arr[:, :].astype(float)

def circumcircle(p1, p2, p3, eps=1e-14):
    """
    Ritorna (centro_x, centro_y), r del circumcerchio del triangolo p1,p2,p3.
    Se i punti sono quasi collineari, ritorna None.
    """
    ax, ay = float(p1[0]), float(p1[1])
    bx, by = float(p2[0]), float(p2[1])
    cx, cy = float(p3[0]), float(p3[1])

    d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    if abs(d) < eps:
        return None  # degenerato/quasi collineare

    a2 = ax*ax + ay*ay
    b2 = bx*bx + by*by
    c2 = cx*cx + cy*cy

    ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d
    uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d
    r = np.hypot(ux - ax, uy - ay)
    return (ux, uy), r


def plot_dcel(filenames, out_png="Meshes/Delaunay/delaunay_output.png", colors=None, label_edges=True, draw_circumcircles=False, plot_pts= False):

    if isinstance(filenames, str):
        filenames = [filenames]
    if colors is None:
        palette = ["green", "crimson", "royalblue", "orange", "purple"]
        colors = palette[:len(filenames)]

    # plotta su stesso asse
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_aspect('equal')
    ax.axis('off')
    for file_idx, filename in enumerate(filenames): 
        color = colors[file_idx % len(colors)]
          
        # Carica il file JSON esportato dalla DCEL
        with open(filename, 'r') as file:
            data = json.load(file)

        # Estrai i nodi
        nodes = {node["id"]: np.array(node["coords"]) for node in data["nodes"]}
        #for nid, p in nodes.items():ax.text(p[0], p[1],str(nid),fontsize=7,color="blue",ha="center",va="center",zorder=5,path_effects=[pe.withStroke(linewidth=2, foreground="white")])

        # Pre-calcolo scala per l’offset delle etichette (dipende dalle dimensioni del dominio)
        all_pts = np.vstack(list(nodes.values()))
    
        span = float(np.ptp(all_pts, axis=0).max())

        # offset ~1% della dimensione massima (regolabile)
        label_offset = 0.0125 * (span if span > 0 else 1.0)

        # Disegna ogni arco + etichetta id
        for edge in data["edges"]:
            from_id = edge["from"]
            to_id = edge["to"]
            is_sub = edge.get("segment", False)

            p1 = nodes[from_id]
            p2 = nodes[to_id]

            #else "blue" if( str(edge.get("id", ""))=="135" )  else "deepskyblue" if (str(edge.get("id", ""))=="529" or str(edge.get("id", ""))=="525" or str(edge.get("id", ""))=="527" or str(edge.get("id", ""))=="155" or str(edge.get("id", ""))=="195" or str(edge.get("id", ""))=="131" or str(edge.get("id", ""))=="95" or str(edge.get("id", ""))=="115")
            color = color if len(filenames) > 1 else "red" if is_sub  else "red"
            linewidth = 0.5 if color!="blue" and color!="deepskyblue" else 1.5
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color, linewidth=linewidth, zorder=1)

            if label_edges:
                # punto medio e normale per posizionare l’etichetta "sopra" l’arco
                mid = 0.5 * (p1 + p2)
                dirv = p2 - p1
                n = np.array([-dirv[1], dirv[0]], dtype=float)  # normale 2D
                n_norm = np.linalg.norm(n)
                if n_norm > 0:
                    n /= n_norm
                text_xy = mid + n * label_offset
                if str(edge.get("id", "")) == "309" or str(edge.get("id", "")) == "236":
                 ax.text(text_xy[0], text_xy[1], str(edge.get("id", "")),ha="center", va="center", fontsize=7, color=color, zorder=2,path_effects=[pe.withStroke(linewidth=2, foreground="white")])

        # Disegna triangoli (celle)
        halfedge_to_node = {edge["id"]: edge["from"] for edge in data["edges"]}
        for cell in data["cells"]:
            try:
                cell_nodes = [nodes[halfedge_to_node[edge]] for edge in cell["edges"]]
                x_values = [pt[0] for pt in cell_nodes] + [cell_nodes[0][0]]
                y_values = [pt[1] for pt in cell_nodes] + [cell_nodes[0][1]]
                ax.plot(x_values, y_values, 'g-', linewidth=0.5, zorder=0)
                # Calcola il centroide della cella (media dei vertici)
                centroid = np.mean(cell_nodes, axis=0)

                # Disegna l'etichetta con l'ID della cella in blu
                #if str(cell.get("id", "")) == "184" or str(cell.get("id", "")) == "182" or str(cell.get("id", "")) == "178" or str(cell.get("id", "")) == "185" or str(cell.get("id", "")) == "163" or str(cell.get("id", "")) == "16" or str(cell.get("id", "")) == "183" or str(cell.get("id", "")) == "177" or str(cell.get("id", "")) == "13" or str(cell.get("id", "")) == "39":
                #ax.text(centroid[0], centroid[1], str(cell.get("id", "")),ha="center", va="center", fontsize=8, color="blue", zorder=3,path_effects=[pe.withStroke(linewidth=2, foreground="white")])

                verts = []
                hids = cell.get("edges", [])
                try:
                    for hid in hids:
                        nid = halfedge_to_node[hid]
                        verts.append(np.asarray(nodes[nid], dtype=float))
                except KeyError as e:
                    print(f"Error in drawing cell {cell.get('id','?')} — missing key: {e}")
                    continue
                # circumcerchio (solo se triangolo)
                if draw_circumcircles and len(verts) == 3:
                    cc = circumcircle(verts[0], verts[1], verts[2])
                    if cc is not None:
                        (ux, uy), r = cc
                        circ = Circle((ux, uy), r, fill=False, linewidth=1,
                                    alpha=0.5, zorder=0, linestyle='-')
                        ax.add_patch(circ)

            except KeyError:
                print(f"Error in drawing cell with edges {cell['edges']}")

        pts = _load_xy_txt("Meshes/Delaunay/data_points.txt")
        if pts is not None and pts.size > 0 and plot_pts:
            ax.scatter(pts[:,0], pts[:,1], s=10, c="crimson", marker="o", edgecolors="white", linewidths=0.3, zorder=4, label="data")

        #from scipy.io import mmread
        #M = mmread("Meshes/Delaunay/locs.mtx")
        #if hasattr(M, "toarray"): M = M.toarray()
        #ax.scatter(M[:, 0], M[:, 1],s=15,c="black",marker="o",edgecolors="white",linewidths=0.3,zorder=3)

        '''
        pts = _load_xy_txt("Meshes/Delaunay/data_points2.txt")
        if pts is not None and pts.size > 0 and plot_pts:
            X = pts[:, 0]
            Y = pts[:, 1]
            DX = pts[:, 2]
            DY = pts = _load_xy_txt("Meshes/Delaunay/data_points2.txt")[:, 3] 
            ax.scatter(X, Y, s=18, c="blue", marker="o", edgecolors="white", linewidths=0.3, zorder=4, label="data")
            scale_factor = 10 # Esempio: prova diversi valori (10, 50, 100)
            ax.quiver(X, Y, DX, DY, color='red', scale=scale_factor, scale_units='xy', angles='xy',zorder=1,  width=0.005)
        '''
        '''
        pts=np.array([
            [0, 884.658],
            [332.985, 609.193],
            [539.077, 795.816],
            [681.723, 1065.16],
            [340.862, 1123.36],
            [130.327, 1250.53],
            [0, 1181.56]  
        ])
    if pts is not None and pts.size > 0 :
        ax.scatter(pts[:,0], pts[:,1], s=18, c="blue", marker="o", edgecolors="white", linewidths=0.3, zorder=4, label="data")
        for i, (x, y) in enumerate(pts, start=0):
            ax.text(x, y, str(i),fontsize=8, color="black",ha="center", va="bottom", zorder=5)
    '''

    # Salvataggio immagine pulita
    plt.savefig(out_png, dpi=300, bbox_inches='tight', pad_inches=0, transparent=True)
    print(f"Plot salvato come '{out_png}'")







if __name__ == "__main__":
    #plot_dcel("Meshes/Delaunay/delaunay_output.json",out_png="Meshes/Delaunay/Meshes new/50adapt_mix4_75_dati.png", label_edges=False,draw_circumcircles=False, plot_pts=True)
    #plot_dcel("Meshes/Delaunay/delaunay_output.json",out_png="Meshes/Delaunay/Meshes new/50adapt_mix4_75.png", label_edges=False,draw_circumcircles=False, plot_pts=False)
    plot_dcel("Meshes/Delaunay/delaunay_output.json",out_png="delaunay_output.png", label_edges=False,draw_circumcircles=False, plot_pts=False)

