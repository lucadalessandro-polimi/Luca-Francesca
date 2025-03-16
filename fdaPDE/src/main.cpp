#include <iostream>
#include <fstream>  // Per scrivere su file
#include "src/geometry/dcel.h"
#include "src/geometry/delaunay.h"
#include <Eigen/Dense>  // Assicuriamoci che Eigen sia incluso

using coords_t = Eigen::Matrix<double, 2, 1>;  // Vettore 2D


int main() {
    using namespace fdapde;

    // Definiamo il bordo del dominio
    Eigen::Matrix<double,5 , 2> boundary;
    boundary << 0.0, 0.0,
                2.0, 0.0,
                2.0, 2.0,
                1.0, 1.0,  //casella delle lettere
                0.0, 2.0;




    // Definiamo alcuni punti interni strategici
    Eigen::Matrix<double, 1, 2> internal;
    internal << 0.5, 0.75;
                //1.5, 0.4,
                //1.0, 1.2,
                //1.0, 1.8;
                //1.7, 1.3,
                //1.0, 0.4;

    // Creiamo l'oggetto Delaunay con il dominio
    Delaunay<2, 2> delaunay(boundary, internal);

    auto find_halfedge_from_id = [&delaunay](int id) -> DCEL<2, 2>::halfedge_t* {
        for (auto it = delaunay.dcel().halfedges_begin(); it != delaunay.dcel().halfedges_end(); ++it) 
            if (it->id() == id) 
                return std::addressof(*it); 
        return nullptr; 
    };


    delaunay.initialize_triangulation();

    delaunay.build_triangulation();
    // Esportiamo la DCEL per visualizzazione

    delaunay.print_dcel();
    delaunay.dcel().export_to_json("dcel_output.json");


 
    return 0;
}
