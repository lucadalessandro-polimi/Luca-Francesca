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

    // TESTA INSERT_EDGE
    // Troviamo gli half-edges associati ai nodi 2 e 4
    DCEL<2, 2>::halfedge_t* h2 = find_halfedge_from_node(2);
    DCEL<2, 2>::halfedge_t* h4 = find_halfedge_from_node(4);
    if (h2 && h4) {
        auto ok= dcel.insert_edge(h2, h4);
        if (ok!=nullptr && ok!=h2 && ok!=h4) std::cout << "Edge inserted between node 2 and node 4.\n";
    } else {
        std::cerr << "Error: Could not find half-edges for nodes 0 and 3.\n";
    }

    //TESTA REMOVE_EDGE 
    DCEL<2, 2>::halfedge_t* h3 = find_halfedge_from_node(3); 
    auto hr= dcel.remove_edge(h3);
    std::cout << hr->id() << std::endl;
    */

    Eigen::Matrix<double, 5, 2> points;
    points << 0.0, 0.0,
              1.0, 0.0,
              2.0, 0.0,
              2.0, 2.0,
              0.0, 2.0;
    DCEL<2, 2> dcel = DCEL<2, 2>::make_polygon(points);
    // Funzione per trovare un half-edge associato a un nodo specifico
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
    

    //TESTA INSERT_NODE
    DCEL<2, 2>::node_t node1(dcel.n_nodes(),nullptr,false,0.5,1.), node2(dcel.n_nodes()+1,nullptr,false,1.5,1.);
    dcel.insert_node(node1);
    dcel.insert_node(node2);

    DCEL<2, 2>::halfedge_t h_i1(dcel.n_halfedges()+1000,&node1);
    DCEL<2, 2>::halfedge_t h_i2(dcel.n_halfedges()+1001,&node2);
    h_i1.set_cell(nullptr);
    h_i2.set_cell(nullptr);  //SE NO BISOGNA MODIFICARE IL COSTRUTTORE
    node1.set_halfedge(&h_i1);
    node2.set_halfedge(&h_i2);
    DCEL<2, 2>::halfedge_t* he1 = find_halfedge_from_node(1);
    DCEL<2, 2>::halfedge_t* hi1=&h_i1;
    DCEL<2, 2>::halfedge_t* hi2=&h_i2;

    DCEL<2, 2>::halfedge_t* h10= dcel.insert_edge(he1, hi1);
    DCEL<2, 2>::halfedge_t* h12= dcel.insert_edge(hi2, he1);
    DCEL<2, 2>::halfedge_t* h3=dcel.insert_edge(h10->twin(), h12);
    dcel.insert_edge(find_halfedge_from_id(3), find_halfedge_from_id(12));
    dcel.insert_edge(find_halfedge_from_id(14), find_halfedge_from_id(4));
    dcel.insert_edge(find_halfedge_from_id(14), find_halfedge_from_id(3));
    dcel.insert_edge(find_halfedge_from_id(10), find_halfedge_from_id(11));   

    
    dcel.remove_edge(find_halfedge_from_id(10));
    dcel.remove_edge(find_halfedge_from_id(13));
    dcel.remove_edge(find_halfedge_from_id(18));
    dcel.remove_edge(find_halfedge_from_id(21));
    dcel.remove_edge(find_halfedge_from_id(16));
    dcel.remove_edge(find_halfedge_from_id(14));
    
    // Esportiamo la DCEL in un file JSON
    dcel.export_to_json("dcel_output.json");

    
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
