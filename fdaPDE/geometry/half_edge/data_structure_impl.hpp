#ifndef DATA_STRUCTURE_IMPL_HPP
#define DATA_STRUCTURE_IMPL_HPP

#include "data_structure.hpp"
#include <iostream>

// Function to check the twin of a HalfEdge
template <typename T>
bool checkTwin(const HalfEdge<T>* edge) {
    if (edge->getTwin() == nullptr) {
        std::cout << "HalfEdge doesn't have a twin." << std::endl;
        return false;
    }
    if (edge->getTwin()->getTwin() != edge) {
        std::cout << "HalfEdge's twin is not correct." << std::endl;
        return false;
    }
    return true;
}

// Function to check the next of a HalfEdge
template <typename T>
bool checkNext(const HalfEdge<T>* edge) {
    if (edge->getNext() == nullptr) {
        std::cout << "HalfEdge doesn't have a next." << std::endl;
        return false;
    }
    if (edge->getNext()->getPrev() != edge) {
        std::cout << "The next HalfEdge's prev is not correct." << std::endl;
        return false;
    }
    return true;
}

// Function to check the prev of a HalfEdge
template <typename T>
bool checkPrev(const HalfEdge<T>* edge) {
    if (edge->getPrev() == nullptr) {
        std::cout << "HalfEdge doesn't have a prev." << std::endl;
        return false;
    }
    if (edge->getPrev()->getNext() != edge) {
        std::cout << "the precedent HalfEdge's prev is not correct." << std::endl;
        return false;
    }
    return true;
}

// Function to check the connections of a HalfEdge
template <typename T>
void checkHalfEdgeConnections(const HalfEdge<T>* edge) {
    std::cout << "Check HalfEdge's connections:" << std::endl;

    bool twin_ok = checkTwin(edge);
    bool next_ok = checkNext(edge);
    bool prev_ok = checkPrev(edge);

    if (twin_ok && next_ok && prev_ok) {
        std::cout << "All HalfEdge's connections are correct." << std::endl;
    } else {
        std::cout << "Some of the HalfEdge's connections are incorrect." << std::endl;
    }
}

#endif // DATA_STRUCTURE_IMPL_HPP


