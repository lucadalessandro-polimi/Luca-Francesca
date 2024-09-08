#include <iostream>
#include "data_structure.h"

// Funzioni di debug per controllare le connessioni degli HalfEdge
template <int LocalDim, int EmbedDim>
void check_twin(const HalfEdge<LocalDim, EmbedDim>& edge) {
    if (edge.twin() == nullptr) {
        std::cout << "HalfEdge doesn't have a twin." << std::endl;
    } else if (edge.twin()->twin() != &edge) {
        std::cout << "HalfEdge's twin is not correct." << std::endl;
    } else {
        std::cout << "Twin: Correct" << std::endl;
    }
}

template <int LocalDim, int EmbedDim>
void check_next(const HalfEdge<LocalDim, EmbedDim>& edge) {
    if (edge.next() == nullptr) {
        std::cout << "HalfEdge doesn't have a next." << std::endl;
    } else if (edge.next()->prev() != &edge) {
        std::cout << "The next HalfEdge's prev is not correct." << std::endl;
    } else {
        std::cout << "Next: Correct" << std::endl;
    }
}

template <int LocalDim, int EmbedDim>
void check_prev(const HalfEdge<LocalDim, EmbedDim>& edge) {
    if (edge.prev() == nullptr) {
        std::cout << "HalfEdge doesn't have a prev." << std::endl;
    } else if (edge.prev()->next() != &edge) {
        std::cout << "The previous HalfEdge's next is not correct." << std::endl;
    } else {
        std::cout << "Prev: Correct" << std::endl;
    }
}

int main() {
    // Creazione dei nodi (ex-vertici) e degli half-edge per la mesh
    Node<2, 2> node1(1.0, 2.0);
    Node<2, 2> node2(3.0, 4.0);
    Node<2, 2> node3(5.0, 6.0);

    // Creazione degli HalfEdge
    HalfEdge<2, 2> he1(node1);
    HalfEdge<2, 2> he2(node2);
    HalfEdge<2, 2> he3(node3);

    // Impostazione delle connessioni twin
    he1.set_twin(&he2);
    he2.set_twin(&he1);

    // Impostazione delle connessioni next e prev
    he1.set_next(&he2);
    he2.set_next(&he3);
    he3.set_next(&he1);

    he1.set_prev(&he3);
    he2.set_prev(&he1);
    he3.set_prev(&he2);

    // Esecuzione delle funzioni di debug per controllare le connessioni
    check_twin(he1);
    check_next(he1);
    check_prev(he1);

    check_twin(he2);
    check_next(he2);
    check_prev(he2);

    return 0;
}





