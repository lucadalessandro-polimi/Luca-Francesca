import matplotlib.pyplot as plt

def read_node_file(filename):
    with open(filename) as f:
        lines = f.readlines()

    header = lines[0].strip().split()
    num_points = int(header[0])

    points = {}
    for line in lines[1:num_points+1]:
        parts = line.strip().split()
        idx = int(parts[0])
        x, y = float(parts[1]), float(parts[2])
        points[idx] = (x, y)

    return points

def read_ele_file(filename):
    with open(filename) as f:
        lines = f.readlines()

    header = lines[0].strip().split()
    num_triangles = int(header[0])

    triangles = []
    for line in lines[1:num_triangles+1]:
        parts = line.strip().split()
        v1, v2, v3 = int(parts[1]), int(parts[2]), int(parts[3])
        triangles.append((v1, v2, v3))

    return triangles

def plot_mesh(points, triangles):
    fig, ax = plt.subplots()
    for tri in triangles:
        x = [points[i][0] for i in tri] + [points[tri[0]][0]]
        y = [points[i][1] for i in tri] + [points[tri[0]][1]]
        ax.plot(x, y, 'k-')

    ax.set_aspect('equal')
    ax.set_title("Triangolazione da Triangle")
    plt.grid(True)
    plt.savefig("triangle_mesh.png", dpi=300)
    print("Mesh salvata come triangle_mesh.png")

# === MAIN ===
node_file = "regioni.1.node"
ele_file = "regioni.1.ele"

points = read_node_file(node_file)
triangles = read_ele_file(ele_file)
plot_mesh(points, triangles)
