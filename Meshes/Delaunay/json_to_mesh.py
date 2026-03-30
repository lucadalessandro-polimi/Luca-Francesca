import json
from typing import Dict, List, Tuple

def json_to_medit_mesh(json_path: str, mesh_path: str, dim: int = 2):
    with open(json_path, "r") as f:
        J = json.load(f)

    # ---- Nodes ----
    nodes = sorted(J["nodes"], key=lambda n: n["id"])
    node_by_id: Dict[int, dict] = {n["id"]: n for n in nodes}

    # Map JSON node id -> Medit 1-based index (compact 1..N)
    id_to_medit = {n["id"]: i + 1 for i, n in enumerate(nodes)}
    N = len(nodes)

    # ---- Edges ----
    edges_by_id = {e["id"]: e for e in J["edges"]}

    # ---- Triangles from cells (each cell lists 3 edge ids) ----
    # ---- Triangles from cells (use halfedge order!) ----
    cells = sorted(J["cells"], key=lambda c: c["id"])
    triangles: List[Tuple[int, int, int]] = []

    for c in cells:
        eids = c["edges"]   # already in h, h->next, h->next->next order
        if len(eids) != 3:
            raise ValueError(
                f"Cell {c['id']} is not a triangle: found {len(eids)} edges."
            )
        # vertices in halfedge order: from(h0), from(h1), from(h2)
        verts = [edges_by_id[eid]["from"] for eid in eids]
        if len(set(verts)) != 3:
            raise ValueError(
                f"Cell {c['id']} not a proper triangle loop: {verts}"
            )
        a, b, ccc = verts[0], verts[1], verts[2]
        # ---- OPTIONAL BUT STRONGLY RECOMMENDED: enforce CCW ----
        ax, ay = node_by_id[a]["coords"]
        bx, by = node_by_id[b]["coords"]
        cx, cy = node_by_id[ccc]["coords"]
        area2 = (bx-ax)*(cy-ay) - (by-ay)*(cx-ax)
        if area2 < 0:
            b, ccc = ccc, b   # flip orientation
        triangles.append((a, b, ccc))


    # ---- Constrained segments (Edges in Medit) ----
    # We write ONLY edges with segment=true (your fixed segments: boundary or internal)
    seg_set = set()
    segments: List[Tuple[int, int, int]] = []

    for e in J["edges"]:
        a, b = e["from"], e["to"]
        # edge is boundary if explicitly marked OR twin == -1
        is_boundary = e.get("boundary", False) or e.get("twin", -1) == -1
        # write edge if boundary OR segment
        if not (is_boundary or e.get("segment", False)):
            continue


        a, b = e["from"], e["to"]
        key = (min(a, b), max(a, b))
        if key in seg_set:
            continue
        seg_set.add(key)

        # classification boundary/internal
        # Prefer edge-level "boundary" if present; else fallback to node flags.
        if "boundary" in e:
            is_boundary = bool(e["boundary"])
        else:
            na = node_by_id[a]
            nb = node_by_id[b]
            is_boundary = bool(na.get("boundary", False) and nb.get("boundary", False))

        ref = 1 if is_boundary else 2
        segments.append((key[0], key[1], ref))

    # Deterministic ordering
    segments.sort(key=lambda s: (s[2], s[0], s[1]))
    triangles.sort()

    # ---- Write Medit .mesh ----
    out: List[str] = []
    out.append("MeshVersionFormatted 1")
    #out.append("1")
    out.append("")
    out.append(f"Dimension")
    out.append(str(dim))
    out.append("")

    # Vertices
    out.append(f"Vertices")
    out.append(str(N))
    for n in nodes:
        x, y = n["coords"]
        # Medit wants an integer "ref" per vertex. Use 1 for boundary nodes, 0 otherwise (or always 0).
        vref = 1 if n.get("boundary", False) else 0
        out.append(f"{x} {y} {vref}")
    out.append("")

    # Edges (segments)
    out.append(f"Edges")
    out.append(str(len(segments)))
    for (a, b, ref) in segments:
        out.append(f"{id_to_medit[a]} {id_to_medit[b]} {ref}")
    out.append("")

    # Triangles
    out.append(f"Triangles")
    out.append(str(len(triangles)))
    for (a, b, ccc) in triangles:
        # triangle ref: 0 is fine
        out.append(f"{id_to_medit[a]} {id_to_medit[b]} {id_to_medit[ccc]} 0")
    out.append("")
    out.append("End")
    out.append("")

    with open(mesh_path, "w", newline="\n") as f:
        f.write("\n".join(out))




#if __name__ == "__main__":
#    json_to_medit_mesh("Meshes/Delaunay/delaunay_output.json","Meshes/Delaunay/florida_base_400.mesh")

