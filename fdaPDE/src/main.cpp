#include <iostream>
#include <fstream>  // Per scrivere su file
#include "src/geometry/dcel.h"
#include "src/geometry/delaunay.h"
#include <Eigen/Dense>  // Assicuriamoci che Eigen sia incluso

using coords_t = Eigen::Matrix<double, 2, 1>;  // Vettore 2D


int main() {
    using namespace fdapde;

/*
    Eigen::Matrix<double, 5, 2> points;
    points << 0, 0,  2, 0,  2, 2,  1, 1,  0, 2;
    DCEL<2, 2> dcel = DCEL<2, 2>::make_polygon(points);

    // Esportiamo la DCEL in un file JSON
    dcel.export_to_json("dcel_output.json");
*/

    
    // Creiamo un oggetto Delaunay
    Delaunay<> triangulation;

    Delaunay<>::node_t A(0, false, coords_t(0.0, 0.0));
    Delaunay<>::node_t B(1, false, coords_t(1.0, 0.0));
    Delaunay<>::node_t C(2, false, coords_t(0.0, 1.0));
    Delaunay<>::node_t D_inside(3, false, coords_t(0.2, 0.2));
    Delaunay<>::node_t D_outside(4, false, coords_t(2.0, 2.0));
    

    // Test `in_circle`
    std::cout << "D_inside è dentro il circumcerchio? " 
              << (triangulation.in_circle(A, B, C, D_inside) ? "Sì" : "No") << std::endl;

    std::cout << "D_outside è dentro il circumcerchio? " 
              << (triangulation.in_circle(A, B, C, D_outside) ? "Sì" : "No") << std::endl;

    // Test `is_point_inside_triangle`
    std::cout << "D_inside è dentro il triangolo ABC? " 
              << (triangulation.is_point_inside_triangle(D_inside, A, B, C) ? "Sì" : "No") << std::endl;

    std::cout << "D_outside è dentro il triangolo ABC? " 
              << (triangulation.is_point_inside_triangle(D_outside, A, B, C) ? "Sì" : "No") << std::endl;

    return 0;
}
