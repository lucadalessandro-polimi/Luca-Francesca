#ifndef DATA_STRUCTURE_H
#define DATA_STRUCTURE_H

#include <iostream>
#include <Eigen/Dense>

template <int LocalDim, int EmbedDim>
class HalfEdge;

template <int LocalDim, int EmbedDim>
class Node {
private:
    Eigen::Matrix<double, EmbedDim, 1> p;
    HalfEdge<LocalDim, EmbedDim>* e;

public:
    // Costruttore che accetta le coordinate
    template <typename... Coords>
    Node(Coords... coords) : p(coords...), e(nullptr) {
        static_assert(sizeof...(coords) == EmbedDim, "Numero di coordinate non corretto per lo spazio immerso");
    }

    // Getter per le coordinate
    Eigen::Matrix<double, EmbedDim, 1> point() const { return p; }

    // Getter e setter per l'HalfEdge associato
    void set_half_edge(HalfEdge<LocalDim, EmbedDim>* edge) { e = edge; }
    HalfEdge<LocalDim, EmbedDim>* half_edge() const { return e; }
};

template <int LocalDim, int EmbedDim>
class HalfEdge {
private:
    Node<LocalDim, EmbedDim>* v;
    HalfEdge<LocalDim, EmbedDim>* twin;
    HalfEdge<LocalDim, EmbedDim>* next;
    HalfEdge<LocalDim, EmbedDim>* prev;

public:
    HalfEdge() : v(nullptr), twin(nullptr), next(nullptr), prev(nullptr) {}

    HalfEdge(Node<LocalDim, EmbedDim>& vert, HalfEdge<LocalDim, EmbedDim>* tw = nullptr,
             HalfEdge<LocalDim, EmbedDim>* ne = nullptr, HalfEdge<LocalDim, EmbedDim>* pr = nullptr)
        : v(&vert), twin(tw), next(ne), prev(pr) {}

    // Setter e getter per il nodo (ora usando una reference)
    void set_node(Node<LocalDim, EmbedDim>& vert) { v = &vert; }
    Node<LocalDim, EmbedDim>* node() const { return v; }

    // Setter e getter per gli altri half-edge
    void set_twin(HalfEdge<LocalDim, EmbedDim>* tw) { twin = tw; }
    HalfEdge<LocalDim, EmbedDim>* twin() const { return twin; }

    void set_next(HalfEdge<LocalDim, EmbedDim>* ne) { next = ne; }
    HalfEdge<LocalDim, EmbedDim>* next() const { return next; }

    void set_prev(HalfEdge<LocalDim, EmbedDim>* pr) { prev = pr; }
    HalfEdge<LocalDim, EmbedDim>* prev() const { return prev; }
};

template <int LocalDim, int EmbedDim>
class Edge {
private:
    HalfEdge<LocalDim, EmbedDim>* e1;
    HalfEdge<LocalDim, EmbedDim>* e2;

public:
    Edge(Node<LocalDim, EmbedDim>& v1, Node<LocalDim, EmbedDim>& v2) {
        e1 = new HalfEdge<LocalDim, EmbedDim>(v1);
        e2 = new HalfEdge<LocalDim, EmbedDim>(v2);
        e1->set_twin(e2);
        e2->set_twin(e1);
    }

    ~Edge() {
        delete e1;
        delete e2;
    }

    HalfEdge<LocalDim, EmbedDim>* edge1() const { return e1; }
    HalfEdge<LocalDim, EmbedDim>* edge2() const { return e2; }
};

#endif // DATA_STRUCTURE_H
