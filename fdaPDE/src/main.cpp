#include <iostream>
#include <cmath>
#include <list>
#include <unordered_map>
#include <vector>
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

    /*
    Eigen::Matrix<double, 13, 2> boundary;  //concave
    boundary << 0.5, 1.0, 
                0.3, 0.3,  
                1.0, 0.7,  
                1.7, 0.0, 
                1.85, 0.15, 
                2.0, 0.3,  
                2.0, 1.3,  
                2.0, 2.0,
                1.0, 2.0,
                0.0, 2.0,  
                0.0, 1.7,
                0.0, 1.5,
                0.0, 1.3;*/
    /*
    Eigen::Matrix<double, 18, 2> boundary; //scala
    boundary << 0.0, 0.0, 
                3.0, 0.0,  
                6.0, 0.0,  
                9.0, 0.0, 
                12.0, 0.0,
                12.0, 2.0, 
                9.0, 2.0,  
                9.0, 4.0,
                8.0, 4.0,
                7.0, 4.0,
                6.0, 4.0,  
                6.0, 6.0,
                3.0, 6.0,
                3.0, 8.0,
                0.0, 8.0,
                0.0, 6.0,
                0.0, 4.0,
                0.0, 2.0;*/
    

    Eigen::Matrix<double, 4, 2> internal;
    internal << 5.0, 5.0,
                5.0, 15.0,
                27.0, 11.0,
                34.0, 8.5;

    Eigen::Matrix<double, 1, 2> prova;
    prova << 3.8, 5.0;
                
          
    /*    
    Eigen::Matrix<double, 4, 2> boundary;
    boundary <<  0.0, 0.0,  
                10.0, 0.0,  
                10.0, 10.0,
            //   0.5, 0.5,  
                0.0, 10.0;*/

    Delaunay<2, 2> delaunay(boundary);
  //  Delaunay<2, 2> delaunay(boundary,holes);

    auto find_halfedge_from_id = [&delaunay](int id) -> DCEL<2, 2>::halfedge_t* {
        for (auto it = delaunay.dcel().halfedges_begin(); it != delaunay.dcel().halfedges_end(); ++it) 
            if (it->id() == id) 
                return std::addressof(*it); 
        return nullptr; 
    };

    delaunay.build_triangulation(internal, boundary);
    //delaunay.build_triangulation(1000,boundary);
    //std::cout<<"-------------------------------FLIP DI CONTROLLO-------------------------------"<<std::endl;
    //delaunay.flip();
    delaunay.split_triangle(find_halfedge_from_id(78)->cell());
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
        
        delaunay_100.build_triangulation(num_points,boundary);   // function to test
        
        auto end = high_resolution_clock::now();    // ending measuring time 
        auto duration = duration_cast<milliseconds>(end - start).count();

        std::cout << "elapsed time for " << num_points << " points: " << duration << " ms" << std::endl;
        outfile << num_points << "," << duration << "\n";  

        num_points = 1000;

        start = high_resolution_clock::now();  // starting measuring time 
        
        delaunay_1000.build_triangulation(num_points,boundary);   // function to test
        
        end = high_resolution_clock::now();    // ending measuring time 
        duration = duration_cast<milliseconds>(end - start).count();

        std::cout << "elapsed time for " << num_points << " points: " << duration << " ms" << std::endl;
        outfile << num_points << "," << duration << "\n";  

        num_points = 10000;

        start = high_resolution_clock::now();  // starting measuring time 
        
        delaunay_10000.build_triangulation(num_points,boundary);   // function to test
        
        end = high_resolution_clock::now();    // ending measuring time 
        duration = duration_cast<milliseconds>(end - start).count();

        std::cout << "elapsed time for " << num_points << " points: " << duration << " ms" << std::endl;
        outfile << num_points << "," << duration << "\n"; 

    outfile.close();
*/  
    return 0;
}
