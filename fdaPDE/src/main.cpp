#include <iostream>
#include <fstream>  // Per scrivere su file
#include "src/geometry/dcel.h"


int main() {
    using namespace fdapde;
/*
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
    std::cout << "successo!" << std::endl;
*/

    Eigen::Matrix<double, 5, 2> points;
    points << 0, 0,  2, 0,  2, 2,  1, 1,  0, 2;
    DCEL<2, 2> dcel = DCEL<2, 2>::make_polygon(points);

    // Esportiamo la DCEL in un file JSON
    dcel.export_to_json("dcel_output.json");


    return 0;
}
