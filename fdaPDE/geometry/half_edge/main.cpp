#include <iostream>
#include <list>
#include <Eigen/Dense>
#include "data_structure.hpp"

using namespace Eigen;

int main() {
    // Definizione dei punti
    Vector3d v1(1.0, 2.0, 3.0);
    Vector3d v2(4.0, 5.0, 6.0);
    Vector3d v3(7.0, 8.0, 9.0);

    // Creazione dei Vertex
    Vertex<Vector3d> vertex1(v1);
    Vertex<Vector3d> vertex2(v2);
    Vertex<Vector3d> vertex3(v3);

    // Creazione degli HalfEdge
    HalfEdge<Vector3d> he1, he2, he3;

    // Impostare i Vertex per gli HalfEdge
    he1.setVertex(&vertex1);
    he2.setVertex(&vertex2);
    he3.setVertex(&vertex3);

    // Collegamenti twin
    he1.setTwin(&he2);
    he2.setTwin(&he1);
    he2.setTwin(&he3);
    he3.setTwin(&he2);

    // Collegamenti next e prev per formare un ciclo
    he1.setNext(&he2);
    he2.setNext(&he3);
    he3.setNext(&he1);

    he1.setPrev(&he3);
    he2.setPrev(&he1);
    he3.setPrev(&he2);

    // Creazione delle liste
    std::list<Vertex<Vector3d>> vertexList;
    std::list<HalfEdge<Vector3d>> halfEdgeList;

    // Aggiunta di Vertex alla lista
    vertexList.push_back(vertex1);
    vertexList.push_back(vertex2);
    vertexList.push_back(vertex3);

    // Aggiunta di HalfEdge alla lista
    halfEdgeList.push_back(he1);
    halfEdgeList.push_back(he2);
    halfEdgeList.push_back(he3);

    // Stampa dei Vertex
    std::cout << "Vertex List:" << std::endl;
    for (const auto& vertex : vertexList) {
        std::cout << "Vertex Point: " << vertex.getPoint().transpose() << std::endl;
    }

    // Funzione per stampare i collegamenti degli HalfEdge
    auto printHalfEdgeConnections = [](const HalfEdge<Vector3d>& he) {
        std::cout << "HalfEdge Vertex: " << he.getVertex()->getPoint().transpose() << std::endl;
        if (he.getTwin()) {
            std::cout << "HalfEdge Twin Vertex: " << he.getTwin()->getVertex()->getPoint().transpose() << std::endl;
        }
        if (he.getNext()) {
            std::cout << "HalfEdge Next Vertex: " << he.getNext()->getVertex()->getPoint().transpose() << std::endl;
        }
        if (he.getPrev()) {
            std::cout << "HalfEdge Prev Vertex: " << he.getPrev()->getVertex()->getPoint().transpose() << std::endl;
        }
    };

    // Stampa degli HalfEdge e le loro connessioni
    std::cout << "HalfEdge List:" << std::endl;
    for (const auto& halfEdge : halfEdgeList) {
        printHalfEdgeConnections(halfEdge);
    }

    return 0;
}





