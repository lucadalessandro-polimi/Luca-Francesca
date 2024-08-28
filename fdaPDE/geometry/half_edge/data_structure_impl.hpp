#ifndef DATA_STRUCTURE_IMPL_HPP
#define DATA_STRUCTURE_IMPL_HPP

#include "data_structure.hpp"
#include <iostream>

// Funzione per controllare il twin di un HalfEdge
template <typename T>
bool checkTwin(const HalfEdge<T>* edge) {
    if (edge->getTwin() == nullptr) {
        std::cout << "HalfEdge non ha un twin." << std::endl;
        return false;
    }
    if (edge->getTwin()->getTwin() != edge) {
        std::cout << "Il twin dell'HalfEdge non è corretto." << std::endl;
        return false;
    }
    return true;
}

// Funzione per controllare il next di un HalfEdge
template <typename T>
bool checkNext(const HalfEdge<T>* edge) {
    if (edge->getNext() == nullptr) {
        std::cout << "HalfEdge non ha un next." << std::endl;
        return false;
    }
    if (edge->getNext()->getPrev() != edge) {
        std::cout << "Il prev dell'HalfEdge successivo non è corretto." << std::endl;
        return false;
    }
    return true;
}

// Funzione per controllare il prev di un HalfEdge
template <typename T>
bool checkPrev(const HalfEdge<T>* edge) {
    if (edge->getPrev() == nullptr) {
        std::cout << "HalfEdge non ha un prev." << std::endl;
        return false;
    }
    if (edge->getPrev()->getNext() != edge) {
        std::cout << "Il next dell'HalfEdge precedente non è corretto." << std::endl;
        return false;
    }
    return true;
}

// Funzione per controllare i collegamenti dell'HalfEdge
template <typename T>
void checkHalfEdgeConnections(const HalfEdge<T>* edge) {
    std::cout << "Verifica dei collegamenti del HalfEdge:" << std::endl;

    bool twin_ok = checkTwin(edge);
    bool next_ok = checkNext(edge);
    bool prev_ok = checkPrev(edge);

    if (twin_ok && next_ok && prev_ok) {
        std::cout << "Tutti i collegamenti del HalfEdge sono corretti." << std::endl;
    } else {
        std::cout << "Alcuni collegamenti del HalfEdge sono incorretti." << std::endl;
    }
}

#endif // DATA_STRUCTURE_IMPL_HPP


