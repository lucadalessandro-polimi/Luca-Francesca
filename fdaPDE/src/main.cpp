#include <iostream>
#include <fstream>  // Per scrivere su file
#include "src/geometry/dcel.h"
#include "src/geometry/delaunay.h"
#include <Eigen/Dense>  // Assicuriamoci che Eigen sia incluso

using coords_t = Eigen::Matrix<double, 2, 1>;  // Vettore 2D


int main() {
    using namespace fdapde;
/*
    // Definiamo il bordo del dominio
    Eigen::Matrix<double,5 , 2> boundary;
    boundary << 0.0, 0.0,
                2.0, 0.0,
                2.0, 2.0,
                1.0, 1.0,  //casella delle lettere
                0.0, 2.0;
*/

// Definizione del bordo della stella (10 vertici)
Eigen::Matrix<double, 10, 2> boundary;
boundary << 1.0, 2.0,  // Punto superiore
            2.0, 2.0,  // Incavo in alto a destra
            2.0, 1.3,  // Punta destra
            2.0, 0.3,  // Incavo in basso a destra
            1.7, 0.3,  // Punta inferiore destra
            1.0, 0.7,  // Incavo inferiore
            0.3, 0.3,  // Punta inferiore sinistra
            0.5, 1.0,  // Incavo in basso a sinistra
            0.0, 1.3,  // Punta sinistra
            0.0, 2.0;  // Incavo in alto a sinistra




    // Definiamo alcuni punti interni strategici
    Eigen::Matrix<double, 2, 2> internal;
    internal << 0.6, 0.6,
                1.5, 1.5;
             //   1.0, 0.2;
               // 1.0, 0.8;
            //    0.2, 1.5;
         //       1.0, 0.4;

    // Creiamo l'oggetto Delaunay con il dominio
    Delaunay<2, 2> delaunay(boundary, internal);

    auto find_halfedge_from_id = [&delaunay](int id) -> DCEL<2, 2>::halfedge_t* {
        for (auto it = delaunay.dcel().halfedges_begin(); it != delaunay.dcel().halfedges_end(); ++it) 
            if (it->id() == id) 
                return std::addressof(*it); 
        return nullptr; 
    };


    delaunay.initialize_triangulation();
    delaunay.dcel().remove_edge(find_halfedge_from_id(20));
    //delaunay.build_triangulation();
    // Esportiamo la DCEL per visualizzazione

    delaunay.print_dcel();
    delaunay.dcel().export_to_json("dcel_output.json");


 
    return 0;
}
