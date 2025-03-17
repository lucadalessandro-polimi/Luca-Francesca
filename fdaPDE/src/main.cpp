#include <iostream>
#include <fstream>  // Per scrivere su file
#include "src/geometry/dcel.h"
#include "src/geometry/delaunay.h"
#include <Eigen/Dense>  
#include <cmath>

using coords_t = Eigen::Matrix<double, 2, 1>;  // Vettore 2D
constexpr int N = 10; // Numero di punti del bordo esterno
constexpr double R = 1.0; // Raggio del dominio principale
constexpr double cx = 1.0, cy = 1.0; // Centro del dominio

constexpr int M = 8; // Numero di punti del buco
constexpr double r_hole = 0.4; // Raggio del buco
constexpr double cx_hole = 1.0, cy_hole = 1.0; // Centro del buco (coincide con il dominio)


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
/*
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
            0.0, 2.0;  // Incavo in alto a sinistra*/

// **Creazione del bordo esterno (antiorario)**
Eigen::Matrix<double, N, 2> boundary;
for (int i = 0; i < N; ++i) {
    double theta = 2.0 * M_PI * i / N;
    boundary(i, 0) = cx + R * cos(theta);
    boundary(i, 1) = cy + R * sin(theta);
}

// **Creazione del buco (orario)**
Eigen::Matrix<double, M, 2> hole;
for (int i = 0; i < M; ++i) {
    double theta = -2.0 * M_PI * i / M; // Cambio segno per ottenere senso orario
    hole(i, 0) = cx_hole + r_hole * cos(theta);
    hole(i, 1) = cy_hole + r_hole * sin(theta);
}

    // **Mettiamo il buco in un vettore di buchi**
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2>> holes = {hole};


    // Definiamo alcuni punti interni strategici
    Eigen::Matrix<double, 1, 2> internal;
    internal << 1.0, 1.5;
              //  1.5, 1.5;
             //   1.0, 0.2;
               // 1.0, 0.8;
            //    0.2, 1.5;
         //       1.0, 0.4;

    // Creiamo l'oggetto Delaunay con il dominio
    //CONVENZIONE : ANTIORARIA PER BORDO E ORARIA PER BUCHI, GESTITA MANUALMENTE ALL'INGRESSO CON LE MATRICI DALLO USER 
    Delaunay<2, 2> delaunay(boundary, internal, holes);
    //Delaunay<2, 2> delaunay(boundary, internal);

    auto find_halfedge_from_id = [&delaunay](int id) -> DCEL<2, 2>::halfedge_t* {
        for (auto it = delaunay.dcel().halfedges_begin(); it != delaunay.dcel().halfedges_end(); ++it) 
            if (it->id() == id) 
                return std::addressof(*it); 
        return nullptr; 
    };

   // delaunay.dcel().add_polygon(find_halfedge_from_id(0),internal);

    //delaunay.initialize_triangulation();
   // delaunay.dcel().remove_edge(find_halfedge_from_id(20));
    //delaunay.build_triangulation();
    // Esportiamo la DCEL per visualizzazione

    delaunay.print_dcel();
    delaunay.dcel().export_to_json("dcel_output.json");


 
    return 0;
}
