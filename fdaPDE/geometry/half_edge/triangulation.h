#include <map> // Per mappare gli edge alle loro twin

template <int LocalDim, int EmbedDim>
class Triangulation {
private:
    std::list<Node<LocalDim, EmbedDim>> nodes; // Lista dei nodi (ex-vertici)
    std::list<std::array<int, 3>> faces; // Lista dei triangoli, rappresentati da tre indici

public:
    // Costruttore
    Triangulation(const std::list<Node<LocalDim, EmbedDim>>& verts,
                  const std::list<std::array<int, 3>>& fac)
        : nodes(verts), faces(fac) {}

    // Metodo per convertire la mesh in una rappresentazione half-edge
    std::list<HalfEdge<LocalDim, EmbedDim>> convert_to_half_edge();
};

// Implementazione del metodo convert_to_half_edge
template <int LocalDim, int EmbedDim>
std::list<HalfEdge<LocalDim, EmbedDim>> Triangulation<LocalDim, EmbedDim>::convert_to_half_edge() {
    std::list<HalfEdge<LocalDim, EmbedDim>> half_edges;
    std::map<std::pair<int, int>, HalfEdge<LocalDim, EmbedDim>*> edge_map; // Per memorizzare gli half-edge gemelli

    // Iterare su ogni triangolo della mesh (face-based)
    for (const auto& face : faces) {
        int v1_idx = face[0]; // Primo vertice del triangolo
        int v2_idx = face[1]; // Secondo vertice del triangolo
        int v3_idx = face[2]; // Terzo vertice del triangolo

        // Prendere i nodi corrispondenti ai vertici
        auto node1 = std::next(nodes.begin(), v1_idx);
        auto node2 = std::next(nodes.begin(), v2_idx);
        auto node3 = std::next(nodes.begin(), v3_idx);

        // Creare i tre half-edge del triangolo
        HalfEdge<LocalDim, EmbedDim>* he1 = new HalfEdge<LocalDim, EmbedDim>(*node1);
        HalfEdge<LocalDim, EmbedDim>* he2 = new HalfEdge<LocalDim, EmbedDim>(*node2);
        HalfEdge<LocalDim, EmbedDim>* he3 = new HalfEdge<LocalDim, EmbedDim>(*node3);

        // Impostare il ciclo degli half-edge per il triangolo (next e prev)
        he1->set_next(he2);
        he2->set_next(he3);
        he3->set_next(he1);

        he1->set_prev(he3);
        he2->set_prev(he1);
        he3->set_prev(he2);

        // Aggiungere gli half-edge alla lista
        half_edges.push_back(*he1);
        half_edges.push_back(*he2);
        half_edges.push_back(*he3);

        // Memorizzare gli half-edge per poter creare i collegamenti gemelli (twin)
        std::pair<int, int> edge1_key = std::minmax(v1_idx, v2_idx); // Spigolo tra v1 e v2
        std::pair<int, int> edge2_key = std::minmax(v2_idx, v3_idx); // Spigolo tra v2 e v3
        std::pair<int, int> edge3_key = std::minmax(v3_idx, v1_idx); // Spigolo tra v3 e v1

        // Collegare half-edge gemelli (twin) se lo spigolo è già stato visto
        if (edge_map.find(edge1_key) != edge_map.end()) {
            HalfEdge<LocalDim, EmbedDim>* twin_edge = edge_map[edge1_key];
            he1->set_twin(twin_edge);
            twin_edge->set_twin(he1);
        } else {
            edge_map[edge1_key] = he1;
        }

        if (edge_map.find(edge2_key) != edge_map.end()) {
            HalfEdge<LocalDim, EmbedDim>* twin_edge = edge_map[edge2_key];
            he2->set_twin(twin_edge);
            twin_edge->set_twin(he2);
        } else {
            edge_map[edge2_key] = he2;
        }

        if (edge_map.find(edge3_key) != edge_map.end()) {
            HalfEdge<LocalDim, EmbedDim>* twin_edge = edge_map[edge3_key];
            he3->set_twin(twin_edge);
            twin_edge->set_twin(he3);
        } else {
            edge_map[edge3_key] = he3;
        }
    }

    return half_edges; // Restituiamo la lista degli half-edge creati
}
