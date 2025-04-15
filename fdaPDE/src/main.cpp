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

double rho_from_min_angle(double theta_min_deg) {
    if (theta_min_deg <= 0.0 || theta_min_deg >= 180.0) {
        return std::numeric_limits<double>::infinity();
    }

    double theta_rad = theta_min_deg * M_PI / 180.0;
    return 1.0 / (2.0 * std::sin(theta_rad));
}

int main() {
    
    Eigen::Matrix<double, 4, 2> boundary;
    boundary << 0.0, 0.0,
                //10.0, 0.0,
                //20.0, 0.0,
                //30.0, 0.0,
                40.0, 0.0,  // Lato inferiore con 4 punti intermedi
    
                //40.0, 5.0,
                //40.0, 10.0,
                //40.0, 15.0,
                40.0, 20.0,  // Lato destro con 3 punti intermedi
    
                //30.0, 20.0,
                //20.0, 20.0,
                //10.0, 20.0,
                0.0, 20.0; // Lato superiore con 4 punti intermedi
    
                //0.0, 15.0,
                //0.0, 10.0,
                //0.0, 5.0;  // Lato sinistro con 3 punti intermedi

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
    

    Eigen::Matrix<double, 14, 2> internal;
    internal << //1.0, 4.0,
                9.0, 4.0,
                11.0, 4.0,
                //19.0, 1.0,
                //21.0, 1.0,
                29.0, 4.0,
                31.0, 4.0,
                39.0, 4.0,
                //1.0, 16.0,
                9.0, 16.0,
                11.0, 16.0,
                //19.0, 19.0,
                //21.0, 19.0,
                29.0, 16.0,
                31.0, 16.0,
                //39.0, 16.0,
                1.0, 4.0,
                //1.0, 6.0,
                4.0, 9.0,
                4.0, 11.0,
                //1.0, 14.0,
                //1.0, 16.0,
                //1.0, 19.0,
                //39.0, 4.0,
                //39.0, 6.0,
                36.0, 9.0,
                36.0, 11.0;
                //39.0, 14.0,
                //39.0, 16.0,
                //39.0, 19.0;
                
                
          
/*        
    Eigen::Matrix<double, 4, 2> boundary;
    boundary <<  0.0, 0.0,  
                10.0, 0.0,  
                10.0, 10.0,
            //   0.5, 0.5,  
                0.0, 10.0;*/

    Eigen::Matrix<double, 8, 2> boundary_concave_ref;
    boundary_concave_ref <<  0, 0, 
                 4, 0, 
                 4, 2, 
                 2.5, 1.5, 
                 4, 4, 
                 0, 4, 
                 0, 2, 
                 1.5, 2.5;
    Eigen::Matrix<double, 1, 2> int_pt;
    int_pt << 0., 0.;

    Eigen::Matrix<double,10,2> star_points_ref;
    star_points_ref << 
    0.0, 1.0,               // Punto 1: punta superiore
    -0.2245, 0.3090,        // Punto 2
    -0.9511, 0.3090,        // Punto 3: punta sinistra
    -0.3633, -0.1180,       // Punto 4
    -0.5878, -0.8090,       // Punto 5: punta inferiore sinistra
    0.0, -0.3820,           // Punto 6
    0.5878, -0.8090,        // Punto 7: punta inferiore destra
    0.3633, -0.1180,        // Punto 8
    0.9511, 0.3090,         // Punto 9: punta destra
    0.2245, 0.3090;  
    

    
    //Delaunay<2, 2> mesh(star_points_ref,int_pt);  //non aggiunge nodi in più perchè sto usando initialize_internal con polygon
    //Delaunay<2,2>::Ruppert_refinement(mesh, rho_from_min_angle(20.0));

    /*
    fdapde::DCEL<2, 2> dcel = fdapde::DCEL<2, 2>::make_polygon(boundary_concave_ref);
    auto find_node_from_id = [&dcel](int id) -> DCEL<2, 2>::node_t* {
        for (auto it = dcel.nodes_begin(); it != dcel.nodes_end(); ++it) 
            if (it->id() == id) 
                return std::addressof(*it); 
        return nullptr; 
    };
    std::cout << (fdapde::internals::is_angle_acute(find_node_from_id(0)->coords(),find_node_from_id(1)->coords(),find_node_from_id(2)->coords(),90.0) ? "yes" : "no") << std::endl;
    */

    //------------------------------- REFINEMENT POTENZIATO ESEMPI----------------------------------------------------------------
    
    Eigen::Matrix<double, 6, 2> triang_isoscele;
    triang_isoscele << 
    -10, 0.0,      // B
    //-0.5, 1.,      // M_AB
    -0.7, std::sqrt(3)/2,
     0.0, 2.1,      // A
     0.5, 1.,      // M_AC
     1.1, 0.0,      // C
     -0.55, 0.0;

    Eigen::Matrix<double, 1, 2> int_is;
    int_is << 1.0, 1.0;
             //0.1, 1.5,
             //0.7, 0.5;
            //-0.6, std::sqrt(3)/20;
    Delaunay<2, 2> mesh(triang_isoscele,int_is);

    /*auto find_halfedge_from_id = [&mesh](int id) -> DCEL<2, 2>::halfedge_t* {
        for (auto it = delaunay.dcel().halfedges_begin(); it != delaunay.dcel().halfedges_end(); ++it) 
            if (it->id() == id) 
                return std::addressof(*it); 
        return nullptr; 
    };
    auto find_halfedge_from_id = [&dcel](int id) -> DCEL<2, 2>::halfedge_t* {
            for (auto it = dcel.halfedges_begin(); it != dcel.halfedges_end(); ++it) 
                if (it->id() == id) 
                    return std::addressof(*it); 
            return nullptr; 
        };
    std::cout << (is_edge_seditious(find_halfedge_from_id(10)) ? "yes" : "no") << std::endl;*/


    //std::cout<<"-------------------------------FLIP DI CONTROLLO-------------------------------"<<std::endl;
    //delaunay.flip();
    
    

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
