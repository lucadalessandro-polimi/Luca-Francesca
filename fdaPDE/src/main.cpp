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
    boundary << 0.0, 1.3,
                0.5, 1.0, 
                0.3, 0.3,  
                1.0, 0.7,  
                1.7, 0.0, 
                1.85, 0.15, 
                2.0, 0.3,  
                2.0, 1.3,  
                1.0, 1.3,
                1.0, 2.0,
                0.0, 2.0,  
                0.0, 1.7,
                0.0, 1.5;
    Eigen::Matrix<double, 4, 2> hole;
    hole << 0.5, 1.5,
            0.5, 1.7,
            0.3, 1.7,
            0.3, 1.5;
    */           
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
    

    Eigen::Matrix<double, 1, 2> internal;
    internal << 0.0, 0.0;
                
    
    Eigen::Matrix<double, 108, 2> C_boundary;
    C_boundary << -0.910947171536292, -0.160624564341911,
                -0.869215674226965, -0.316368632576244,
                -0.801073498500606, -0.4625,
                -0.708591109885055, -0.594578538960049,
                -0.594578538960049, -0.708591109885055,
                -0.4625, -0.801073498500606,
                -0.316368632576244, -0.869215674226965,
                -0.160624564341911, -0.910947171536292,
                5.6638043864285e-17, -0.925,
                0.166666666666667, -0.925,
                0.333333333333333, -0.925,
                0.5, -0.925,
                0.666666666666667, -0.925,
                0.833333333333333, -0.925,
                1, -0.925,
                1.16666666666667, -0.925,
                1.33333333333333, -0.925,
                1.5, -0.925,
                1.66666666666667, -0.925,
                1.83333333333333, -0.925,
                2, -0.925,
                2.16666666666667, -0.925,
                2.33333333333333, -0.925,
                2.5, -0.925,
                2.66666666666667, -0.925,
                2.83333333333333, -0.925,
                3, -0.925,
                3.16072704159334, -0.89302940365474,
                3.29698484809835, -0.80198484809835,
                3.38802940365474, -0.665727041593338,
                3.42, -0.505,
                3.38802940365474, -0.344272958406662,
                3.29698484809835, -0.20801515190165,
                3.16072704159334, -0.11697059634526,
                3, -0.085,
                2.83333333333333, -0.085,
                2.66666666666667, -0.085,
                2.5, -0.085,
                2.33333333333333, -0.085,
                2.16666666666667, -0.085,
                2, -0.085,
                1.83333333333333, -0.085,
                1.66666666666667, -0.085,
                1.5, -0.085,
                1.33333333333333, -0.085,
                1.16666666666667, -0.085,
                1, -0.085,
                0.833333333333333, -0.085,
                0.666666666666667, -0.085,
                0.5, -0.085,
                0.333333333333333, -0.085,
                0.166666666666667, -0.085,
                5.2045770037451e-18, -0.085,
                -0.085, 1.04091540074902e-17,
                5.2045770037451e-18, 0.085,
                0.166666666666667, 0.085,
                0.333333333333333, 0.085,
                0.5, 0.085,
                0.666666666666667, 0.085,
                0.833333333333333, 0.085,
                1, 0.085,
                1.16666666666667, 0.085,
                1.33333333333333, 0.085,
                1.5, 0.085,
                1.66666666666667, 0.085,
                1.83333333333333, 0.085,
                2, 0.085,
                2.16666666666667, 0.085,
                2.33333333333333, 0.085,
                2.5, 0.085,
                2.66666666666667, 0.085,
                2.83333333333333, 0.085,
                3, 0.085,
                3.16072704159334, 0.11697059634526,
                3.29698484809835, 0.20801515190165,
                3.38802940365474, 0.344272958406662,
                3.42, 0.505,
                3.38802940365474, 0.665727041593338,
                3.29698484809835, 0.80198484809835,
                3.16072704159334, 0.89302940365474,
                3, 0.925,
                2.83333333333333, 0.925,
                2.66666666666667, 0.925,
                2.5, 0.925,
                2.33333333333333, 0.925,
                2.16666666666667, 0.925,
                2, 0.925,
                1.83333333333333, 0.925,
                1.66666666666667, 0.925,
                1.5, 0.925,
                1.33333333333333, 0.925,
                1.16666666666667, 0.925,
                1, 0.925,
                0.833333333333333, 0.925,
                0.666666666666667, 0.925,
                0.5, 0.925,
                0.333333333333333, 0.925,
                0.166666666666667, 0.925,
                5.6638043864285e-17, 0.925,
                -0.160624564341911, 0.910947171536292,
                -0.316368632576244, 0.869215674226965,
                -0.4625, 0.801073498500606,
                -0.594578538960049, 0.708591109885055,
                -0.708591109885055, 0.594578538960049,
                -0.801073498500606, 0.4625,
                -0.869215674226965, 0.316368632576244,
                -0.910947171536292, 0.160624564341911,
                -0.925, 1.1327608772857e-16;
    
   

    //Delaunay<2, 2> delaunay(C_boundary);
    //HOLES
    //std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2>> holes;
    //holes.push_back(hole);
    //Delaunay<2, 2> delaunay(boundary,holes);

    Eigen::Matrix<double, 4, 2> boundary1;
    boundary1 << 0.0, 0.0,
                40.0, 0.0,  
                40.0, 20.0,  
                0.0, 20.0;

    Eigen::Matrix<double, 4, 2> boundary_clockwise;
    boundary_clockwise << 0.0, 0.0,
                        0.0, 20.0,
                        40.0, 20.0,
                        40.0, 0.0;

    Eigen::Matrix<double, 10, 2> star;
    star <<0.0,    1.0,        // punta in alto
          -0.2245, 0.3090,     // interno
          -0.9511, 0.3090,     // punta sinistra
          -0.3633,-0.1180,     // interno
          -0.5878,-0.8090,     // punta in basso a sinistra
          0.0,   -0.3820,     // interno
          0.5878,-0.8090,     // punta in basso a destra
          0.3633,-0.1180,     // interno
          0.9511, 0.3090,     // punta destra
          0.2245, 0.3090;     // interno

                        
    Delaunay<2, 2> mesh(star,internal);
    /*
    auto dcel = Delaunay<2, 2>::Triangulation_to_DCEL(mesh); //if TRiangulation_to_DCEL is public
    auto find_halfedge_from_id = [&dcel](int id) -> DCEL<2, 2>::halfedge_t* {
        for (auto it = dcel.halfedges_begin(); it != dcel.halfedges_end(); ++it)
            if (it->id() == id)
                return std::addressof(*it);
        return nullptr;
    };
    auto* he = find_halfedge_from_id(0);
    Delaunay<2,2>::split_subsegment(dcel, he);
    */

    //mesh.Ruppert_refinement(2.0);

    



///////////TEST OF COMPUTATIONAL EFFICIENCY////////////////
    /*   
    Delaunay<2, 2> delaunay_100(boundary); 
    Delaunay<2, 2> delaunay_1000(boundary); 
    Delaunay<2, 2> delaunay_10000(boundary); 
    Delaunay<2, 2> delaunay_100000(boundary); 

    std::ofstream outfile("fdaPDE/src/timing_results.csv");
    outfile << "NumPoints,TimeElapsed(ms)\n";
        
        int num_points = 100;

        auto start = high_resolution_clock::now();  // starting measuring time 
        
        delaunay_100.build_triangulation(num_points,boundary);   // function to test
        
        auto end = high_resolution_clock::now();    // ending measuring time 
        auto duration = duration_cast<milliseconds>(end - start).count();

        std::cout << "elapsed time for " << num_points << " points: " << duration << " ms" << std::endl;
        outfile << num_points << "," << duration << "\n";  
        
        
        int num_points = 1000;

        auto start = high_resolution_clock::now();  // starting measuring time 
        
        delaunay_1000.build_triangulation(num_points,boundary);   // function to test
        
        auto end = high_resolution_clock::now();    // ending measuring time 
        auto duration = duration_cast<milliseconds>(end - start).count();

        std::cout << "elapsed time for " << num_points << " points: " << duration << " ms" << std::endl;
        outfile << num_points << "," << duration << "\n";  

        num_points = 10000;

        start = high_resolution_clock::now();  // starting measuring time 
        
        delaunay_10000.build_triangulation(num_points,boundary);   // function to test
        
        end = high_resolution_clock::now();    // ending measuring time 
        duration = duration_cast<milliseconds>(end - start).count();

        std::cout << "elapsed time for " << num_points << " points: " << duration << " ms" << std::endl;
        outfile << num_points << "," << duration << "\n"; 

        num_points = 100000;

        start = high_resolution_clock::now();  // starting measuring time 
        
        delaunay_100000.build_triangulation(num_points,boundary);   // function to test
        
        end = high_resolution_clock::now();    // ending measuring time 
        duration = duration_cast<milliseconds>(end - start).count();

        std::cout << "elapsed time for " << num_points << " points: " << duration << " ms" << std::endl;
        outfile << num_points << "," << duration << "\n"; 

    outfile.close();
    */
    return 0;
}
