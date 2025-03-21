#include <iostream>
#include <cmath>
#include <list>
#include <Eigen/Dense>
#include <cstdlib> 
#include <fstream>
#include <random> 
#include <nlohmann/json.hpp> 
using json = nlohmann::json;
#include "../geometry.h"
#include <chrono>
using namespace std::chrono;
using namespace fdapde;
using coords_t = Eigen::Matrix<double, 2, 1>;  // Vettore 2D

constexpr int N = 10; // Numero di punti del bordo esterno
constexpr double R = 1.0; // Raggio del dominio principale
constexpr double cx = 1.0, cy = 1.0; // Centro del dominio

constexpr int M = 8; // Numero di punti del buco
constexpr double r_hole = 0.4; // Raggio del buco
constexpr double cx_hole = 1.0, cy_hole = 1.0; // Centro del buco (coincide con il dominio)

int main() {

/*
    Eigen::Matrix<double, 16, 2> boundary;
    boundary << 0.0, 0.0,
                10.0, 0.0,
                20.0, 0.0,
                30.0, 0.0,
                40.0, 0.0,  // Lato inferiore con 4 punti intermedi
    
                40.0, 5.0,
                40.0, 10.0,
                40.0, 15.0,
                40.0, 20.0,  // Lato destro con 3 punti intermedi
    
                30.0, 20.0,
                20.0, 20.0,
                10.0, 20.0,
                0.0, 20.0,  // Lato superiore con 4 punti intermedi
    
                0.0, 15.0,
                0.0, 10.0,
                0.0, 5.0;  // Lato sinistro con 3 punti intermedi
                
*/

Eigen::Matrix<double, 4, 2> boundary;
boundary << 0.0, 0.0,
            10.0, 0.0,
            10.0, 10.0,
            0.0,10.0;

    
    /*
    // Definizione del bordo della stella (10 vertici)
    Eigen::Matrix<double, 10, 2> boundary;
    boundary << 1.0, 2.0,  // Punto superiore
                2.0, 2.0,  // Incavo in alto a destra
                2.0, 1.3,  // Punta destra
                2.0, 0.3,  // Incavo in basso a destra
                1.7, 0.0,  // Punta inferiore destra
                1.0, 0.7,  // Incavo inferiore
                0.3, 0.3,  // Punta inferiore sinistra
                0.5, 1.0,  // Incavo in basso a sinistra
                0.0, 1.3,  // Punta sinistra
                0.0, 2.0;  // Incavo in alto a sinistra
    */
  
/*    
    Eigen::Matrix<double, N, 2> boundary;
    for (int i = 0; i < N; ++i) {
        double theta = 2.0 * M_PI * i / N;
        boundary(i, 0) = cx + R * cos(theta);
        boundary(i, 1) = cy + R * sin(theta);
    }

    Eigen::Matrix<double, M, 2> hole;
    for (int i = 0; i < M; ++i) {
        double theta = -2.0 * M_PI * i / M; 
        hole(i, 0) = cx_hole + r_hole * cos(theta);
        hole(i, 1) = cy_hole + r_hole * sin(theta);
    }

    std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2>> holes = {hole};
 */   

Eigen::Matrix<double, 38, 2> internal;
internal << 11.5, 19.7,
            3.6, 4.9,
            9.4, 17.4,
            35.8, 7.6,
            0.3, 5.9,
            24.4, 5.0,
            13.1, 8.6,
            27.1, 1.2,
            17.8, 19.1,
            0.8, 13.1,
            8.6, 3.0,
            22.5, 3.0,
            7.2, 1.6,
            31.0, 10.4,
            19.7, 19.3,
            29.7, 17.7,
            17.7, 10.2,
            10.7, 11.9,
            2.4, 0.3,
            8.7, 5.4,
            32.5, 11.4,
            36.8, 3.2,
            0.3, 9.8,
            19.3, 9.2,
            12.7, 8.9,
            4.4, 10.1,
            29.9, 0.8,
            18.1, 4.3,
            27.6, 5.6,
            9.5, 5.1,
            35.7, 11.7,
            4.0, 18.9,
            2.6, 16.7,
            27.4, 6.5,
            23.8, 1.6,
            15.2, 5.3,
            1.7, 12.6,
            5.4, 19.2;
        //    23.5, 14.3;
    //        5.9, 2.1;


 

    Delaunay<2, 2> delaunay(boundary);
   // delaunay.set_internal_points(internal);
   // delaunay.build_triangulation();
   // std::cout << "PUNTO PROVA TROVATO IN CELLA ID : "<<delaunay.find_triangle(prova)->id() << std::endl;
   // delaunay.insert_vertex(prova, delaunay.find_triangle(prova));


    auto find_halfedge_from_id = [&delaunay](int id) -> DCEL<2, 2>::halfedge_t* {
        for (auto it = delaunay.dcel().halfedges_begin(); it != delaunay.dcel().halfedges_end(); ++it) 
            if (it->id() == id) 
                return std::addressof(*it); 
        return nullptr; 
    };


  //  delaunay.build_triangulation(1000);
    delaunay.dcel().remove_edge(find_halfedge_from_id(7));
    delaunay.print_dcel();
    delaunay.dcel().export_to_json("dcel_output.json");

///////////TEST OF COMPUTATIONAL EFFICIENCY////////////////
/*
    Delaunay<2, 2> delaunay_100(boundary); 
    Delaunay<2, 2> delaunay_1000(boundary); 
    Delaunay<2, 2> delaunay_10000(boundary); 

    std::ofstream outfile("fdaPDE/src/timing_results.csv");
    outfile << "NumPoints,TimeElapsed(ms)\n";



        int num_points = 100;

        auto start = high_resolution_clock::now();  // starting measuring time 
        
        delaunay_100.build_triangulation(num_points);   // function to test
        
        auto end = high_resolution_clock::now();    // ending measuring time 
        auto duration = duration_cast<milliseconds>(end - start).count();

        std::cout << "elapsed time for " << num_points << " points: " << duration << " ms" << std::endl;
        outfile << num_points << "," << duration << "\n";  

        num_points = 1000;

        start = high_resolution_clock::now();  // starting measuring time 
        
        delaunay_1000.build_triangulation(num_points);   // function to test
        
        end = high_resolution_clock::now();    // ending measuring time 
        duration = duration_cast<milliseconds>(end - start).count();

        std::cout << "elapsed time for " << num_points << " points: " << duration << " ms" << std::endl;
        outfile << num_points << "," << duration << "\n";  

        num_points = 10000;

        start = high_resolution_clock::now();  // starting measuring time 
        
        delaunay_10000.build_triangulation(num_points);   // function to test
        
        end = high_resolution_clock::now();    // ending measuring time 
        duration = duration_cast<milliseconds>(end - start).count();

        std::cout << "elapsed time for " << num_points << " points: " << duration << " ms" << std::endl;
        outfile << num_points << "," << duration << "\n"; 

    outfile.close();
    std::cout << "completed test in 'timing_results.txt'." << std::endl;
*/
    return 0;
}
