#include <iostream>
#include <fstream>  // Per scrivere su file
#include "src/geometry/dcel.h"


int main() {
    using namespace fdapde;


    Eigen::Matrix<double, 5, 2> points;
    points << 0, 0,  2, 0,  2, 2,  1, 1,  0, 2;
    DCEL<2, 2> dcel = DCEL<2, 2>::make_polygon(points);

    // Esportiamo la DCEL in un file JSON
    dcel.export_to_json("dcel_output.json");


    return 0;
}
