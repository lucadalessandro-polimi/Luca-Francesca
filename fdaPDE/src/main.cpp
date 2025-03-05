#include <iostream>
#include <fstream>  // Per scrivere su file
#include "src/geometry/dcel.h"

int main() {
    using namespace fdapde;

    // Creiamo i vertici di un esagono regolare di raggio 1
    Eigen::Matrix<double, 6, 2> hexagon_nodes;
    double angle_step = M_PI / 3.0; // 60 gradi tra ogni vertice
    for (int i = 0; i < 6; ++i) {
        double angle = i * angle_step;
        hexagon_nodes(i, 0) = cos(angle);  // X
        hexagon_nodes(i, 1) = sin(angle);  // Y
    }

    // Costruiamo il DCEL dell'esagono
    DCEL<2, 2> hexagon = DCEL<2, 2>::make_polygon(hexagon_nodes);

    // Stampiamo i nodi per verificare
    std::cout << "Nodi dell'esagono:\n" << hexagon_nodes << "\n";

    // Scriviamo i nodi in un file
    std::ofstream file("hexagon.txt");
    if (file.is_open()) {
        for (int i = 0; i < 6; ++i) {
            file << hexagon_nodes(i, 0) << " " << hexagon_nodes(i, 1) << "\n";
        }
        file.close();
        std::cout << "Coordinate salvate in hexagon.txt\n";
    } else {
        std::cerr << "Errore: impossibile scrivere il file.\n";
    }

    return 0;
}
