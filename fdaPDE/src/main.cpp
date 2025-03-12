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
    points << 0.0, 0.0,
              1.0, 0.0,
              2.0, 0.0,
              2.0, 2.0,
              0.0, 2.0;
    DCEL<2, 2> dcel = DCEL<2, 2>::make_polygon(points);
    Eigen::Matrix<double, 3, 2> points_polyg;
    points_polyg << 1.5, 0.5,
                    1.0, 1.0,
                    0.5, 1.5;

    auto find_halfedge_from_node = [&dcel](int node_id) -> DCEL<2, 2>::halfedge_t* {
        for (auto it = dcel.nodes_begin(); it != dcel.nodes_end(); ++it) 
            if (it->id() == node_id) 
                return it->halfedge(); 
        return nullptr; 
    };
    auto find_halfedge_from_id = [&dcel](int id) -> DCEL<2, 2>::halfedge_t* {
        for (auto it = dcel.halfedges_begin(); it != dcel.halfedges_end(); ++it) 
            if (it->id() == id) 
                return std::addressof(*it); 
        return nullptr; 
    };
    dcel.add_polygon(find_halfedge_from_id(4),points_polyg);

    Eigen::Matrix<double, 1, 2> points_triangle;
    points_triangle << 1.5, 1.5;
    dcel.add_polygon(find_halfedge_from_id(15),points_triangle);
   
    Eigen::Matrix<double, 1, 2> points_triangle_2;
    points_triangle_2 << 0.5, 0.5;
    dcel.add_polygon(find_halfedge_from_id(14),points_triangle_2);

    Eigen::Matrix<double, 1, 2> points_triangle_3;
    points_triangle_3 << 0.8, 0.5;  
    dcel.add_polygon(find_halfedge_from_id(25), points_triangle_3);

    Eigen::Matrix<double, 1, 2> points_triangle_4;
    points_triangle_4 << 0.1, 0.1;  
    dcel.add_polygon(find_halfedge_from_id(23), points_triangle_4);

    Eigen::Matrix<double, 1, 2> points_triangle_5;
    points_triangle_5 << 1.4, 0.5;  
    dcel.add_polygon(find_halfedge_from_id(29), points_triangle_5);

    Eigen::Matrix<double, 1, 2> points_triangle_6;
    points_triangle_6 << 0.3, 0.01;  
    dcel.add_polygon(find_halfedge_from_id(27), points_triangle_6);

    dcel.remove_edge(find_halfedge_from_id(16));
    dcel.remove_edge(find_halfedge_from_id(10));
    dcel.remove_edge(find_halfedge_from_id(12));
   
    Delaunay<2,2> triangulation(dcel);
   

    // Punto da inserire
    coords_t new_point;
    new_point << 0.75, 1.0;

    // Troviamo il triangolo contenente il punto
    const DCEL<2,2>::cell_t* triangle = triangulation.find_triangle(new_point);
  
    //iangulation.remove_triangle(triangle);
    

    if (triangle) {
        triangulation.insert_vertex(new_point, triangle);
        std::cout << "Punto inserito con successo!\n";
    } else {
        std::cout << "Errore: il punto non è contenuto in nessun triangolo!\n";
    }

    // Esportiamo per visualizzare la nuova mesh
    triangulation.dcel().export_to_json("dcel_output.json");

  
  
    //bool inside = triangulation.in_circle(Eigen::Vector2d(0.75, 1.0), Eigen::Vector2d(0.5, 0.5), Eigen::Vector2d(0.8, 0.5), Eigen::Vector2d(0.3, 0.01));
    //std::cout << "Punto D è dentro il circumcerchio? " << (inside ? "SI" : "NO") << std::endl;
*/



    
    // Definiamo il bordo del dominio
    Eigen::Matrix<double, 5, 2> boundary;
    boundary << 0.0, 0.0,
                2.0, 0.0,
                2.0, 2.0,
             //   1.0, 1.0, posta delle lettere
                0.0, 2.0,
                0.0, 0.0;  // Chiudiamo il poligono

    // Definiamo alcuni punti interni strategici
    Eigen::Matrix<double, 4, 2> internal;
    internal << 0.5, 0.5,
                1.5, 0.5,
                1.0, 1.2,
                1.0, 1.8;

    // Creiamo l'oggetto Delaunay con il dominio
    Delaunay<2, 2> delaunay(boundary, internal);

    auto find_halfedge_from_id = [&delaunay](int id) -> DCEL<2, 2>::halfedge_t* {
        for (auto it = delaunay.dcel().halfedges_begin(); it != delaunay.dcel().halfedges_end(); ++it) 
            if (it->id() == id) 
                return std::addressof(*it); 
        return nullptr; 
    };

    Eigen::Matrix<double, 1, 2> points_triangle_;
    points_triangle_ << 2.0, 2.0; 
    
    delaunay.dcel().add_polygon(find_halfedge_from_id(0),points_triangle_);
    //delaunay.dcel().add_polygon(find_halfedge_from_id(0),delaunay.dcel().adjacent(find_halfedge_from_id(0))->coords());
    //delaunay.initialize_triangulation();
    //delaunay.build_triangulation();
    // Esportiamo la DCEL per visualizzazione
    delaunay.dcel().export_to_json("dcel_output.json");


 
    return 0;
}
