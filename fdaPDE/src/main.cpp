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
    
    Eigen::Matrix<double, 24, 2> boundary; 
    boundary << 0.0, 0.0,
                5.0, 0.0,  //24
                10.0, 0.0,
                15.0, 0.0,
                20.0, 0.0,
                25.0, 0.0,
                30.0, 0.0,
                35.0, 0.0,
                40.0, 0.0, 
                40.0, 5.0, 
                40.0, 10.0,
                40.0, 15.0,
                40.0, 20.0,
                35.0, 20.0,
                30.0, 20.0,
                25.0, 20.0,
                20.0, 20.0,
                15.0, 20.0,
                10.0, 20.0,
                5.0, 20.0,
                0.0, 20.0,
                0.0, 15.0,
                0.0, 10.0,
                0.0, 5.0;
   

    
    Eigen::Matrix<double, 13, 2> boundary_concave;  
    boundary_concave << 0.0, 1.3,
                0.5, 1.0, 
                0.3, 0.3,  
                1.0, 0.7,  
                1.0, 0.0, 
                1.85, 0.15, 
                2.0, 0.3,  
                2.0, 1.3,  
                1.0, 1.3,
                1.0, 2.0,
                0.0, 2.0,  
                0.0, 1.7,
                0.0, 1.5;
    Eigen::Matrix<double, 4, 2> hole_concave1;
    hole_concave1 << 0.5, 1.5,
            0.5, 1.7,
            0.3, 1.7,
            0.3, 1.5;
            

    Eigen::Matrix<double, 4, 2> hole_concave2;
    hole_concave2 << 1.3, 1.2,
            1.6, 1.2,
            1.6, 0.5,
            1.3, 0.5;    
    
    Eigen::Matrix<double, 10, 2> boundary_stairs; //scala
    boundary_stairs << 0.0, 0.0, 
                //3.0, 0.0,  
                //6.0, 0.0,  
                //9.0, 0.0, 
                12.0, 0.0,
                12.0, 2.0, 
                9.0, 2.0,  
                9.0, 4.0,
                //8.0, 4.0,
                //7.0, 4.0,
                6.0, 4.0,  
                6.0, 6.0,
                3.0, 6.0,
                3.0, 8.0,
                0.0, 8.0;
                //0.0, 6.0,
                //0.0, 4.0,
                //0.0, 3.0;
    

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

    Eigen::Matrix<double,108,2> C;
    C << 
    -0.910947171536292,-0.160624564341911,
    -0.869215674226965,-0.316368632576244,
    -0.801073498500606,-0.4625,
    -0.708591109885055,-0.594578538960049,
    -0.594578538960049,-0.708591109885055,
    -0.4625,-0.801073498500606,
    -0.316368632576244,-0.869215674226965,
    -0.160624564341911,-0.910947171536292,
    5.6638043864285e-17,-0.925,
    0.166666666666667,-0.925,
    0.333333333333333,-0.925,
    0.5,-0.925,
    0.666666666666667,-0.925,
    0.833333333333333,-0.925,
    1,-0.925,
    1.16666666666667,-0.925,
    1.33333333333333,-0.925,
    1.5,-0.925,
    1.66666666666667,-0.925,
    1.83333333333333,-0.925,
    2,-0.925,
    2.16666666666667,-0.925,
    2.33333333333333,-0.925,
    2.5,-0.925,
    2.66666666666667,-0.925,
    2.83333333333333,-0.925,
    3,-0.925,
    3.16072704159334,-0.89302940365474,
    3.29698484809835,-0.80198484809835,
    3.38802940365474,-0.665727041593338,
    3.42,-0.505,
    3.38802940365474,-0.344272958406662,
    3.29698484809835,-0.20801515190165,
    3.16072704159334,-0.11697059634526,
    3,-0.085,
    2.83333333333333,-0.085,
    2.66666666666667,-0.085,
    2.5,-0.085,
    2.33333333333333,-0.085,
    2.16666666666667,-0.085,
    2,-0.085,
    1.83333333333333,-0.085,
    1.66666666666667,-0.085,
    1.5,-0.085,
    1.33333333333333,-0.085,
    1.16666666666667,-0.085,
    1,-0.085,
    0.833333333333333,-0.085,
    0.666666666666667,-0.085,
    0.5,-0.085,
    0.333333333333333,-0.085,
    0.166666666666667,-0.085,
    5.2045770037451e-18,-0.085,
    -0.085,1.04091540074902e-17,
    5.2045770037451e-18,0.085,
    0.166666666666667,0.085,
    0.333333333333333,0.085,
    0.5,0.085,
    0.666666666666667,0.085,
    0.833333333333333,0.085,
    1,0.085,
    1.16666666666667,0.085,
    1.33333333333333,0.085,
    1.5,0.085,
    1.66666666666667,0.085,
    1.83333333333333,0.085,
    2,0.085,
    2.16666666666667,0.085,
    2.33333333333333,0.085,
    2.5,0.085,
    2.66666666666667,0.085,
    2.83333333333333,0.085,
    3,0.085,
    3.16072704159334,0.11697059634526,
    3.29698484809835,0.20801515190165,
    3.38802940365474,0.344272958406662,
    3.42,0.505,
    3.38802940365474,0.665727041593338,
    3.29698484809835,0.80198484809835,
    3.16072704159334,0.89302940365474,
    3,0.925,
    2.83333333333333,0.925,
    2.66666666666667,0.925,
    2.5,0.925,
    2.33333333333333,0.925,
    2.16666666666667,0.925,
    2,0.925,
    1.83333333333333,0.925,
    1.66666666666667,0.925,
    1.5,0.925,
    1.33333333333333,0.925,
    1.16666666666667,0.925,
    1,0.925,
    0.833333333333333,0.925,
    0.666666666666667,0.925,
    0.5,0.925,
    0.333333333333333,0.925,
    0.166666666666667,0.925,
    5.6638043864285e-17,0.925,
    -0.160624564341911,0.910947171536292,
    -0.316368632576244,0.869215674226965,
    -0.4625,0.801073498500606,
    -0.594578538960049,0.708591109885055,
    -0.708591109885055,0.594578538960049,
    -0.801073498500606,0.4625,
    -0.869215674226965,0.316368632576244,
    -0.910947171536292,0.160624564341911,
    -0.925,1.1327608772857e-16;

    


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
    int_is << 0.0, 1.0;
             //0.1, 1.5,
             //0.7, 0.5;
            //-0.6, std::sqrt(3)/20;
    //Delaunay<2, 2> mesh(triang_isoscele,int_is);
    

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
    
    Eigen::Matrix<double, 4, 2> hole;
    hole << 10., 5.,
    10., 15.,
    30., 15.,
    30., 5.;
    Eigen::Matrix<double, 8, 2> hole3;
    hole3 << 10., 5.,
    10., 10.,
    10., 15.,
    20., 15.,
    30., 15.,
    30., 10.,
    30., 5.,
    20., 5.;
    Eigen::Matrix<double, 4, 2> hole2;
    hole2 << 3.,4.,
    3., 11.,
    6., 11.,
    6., 4.;
    Eigen::Matrix<double, 5, 2> hole4;
    hole4 << 36.5, 15.,
    37.5, 16.5,
    36.5, 18.,
    38.5, 18.,
    38.5, 15.;
    
    Eigen::Matrix<double, 1, 2> internal1;
    internal1 << 18.9533, 16.7049;
    
    Eigen::Matrix<double, 6, 2> hole_stairs;
    hole_stairs << 4.0, 2.5,   
         5.0, 2.5,   
         5.3, 3.0,   
         5.0, 3.5,   
         4.0, 3.5,   
         3.7, 3.0; 
    Eigen::Matrix<double, 6, 2> heart_hole;
    heart_hole <<  
    1.5, 6.6,     // punta inferiore
    2.1, 7.1,     // lato curvo destro
    1.8, 7.6,     // top destro
    1.5, 7.4,     // vertice superiore centrale
    1.2, 7.6,     // top sinistro
    0.9, 7.1;     // lato curvo sinistro

            
    
    Delaunay<2, 2> del(boundary_stairs, 0, std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2>> {hole_stairs, heart_hole});
    //del.Ruppert_refinement(2.0);  //COSA NON VA  ?????

///////////TEST OF COMPUTATIONAL EFFICIENCY////////////////


    std::ofstream outfile("fdaPDE/src/timing_results.csv");
    outfile << "NumPoints,TimeElapsed(ms)\n";
    
    /*Delaunay<2, 2> del_10000(boundary,10000); 
    auto start = high_resolution_clock::now();  // starting measuring time 
    //Delaunay<2, 2> del_1000(boundary,1000);  
    del_10000.Ruppert_refinement(2.0);
    auto end = high_resolution_clock::now();    // ending measuring time 
    auto duration = duration_cast<milliseconds>(end - start).count();
    std::cout << "elapsed time for " << del_10000.dcel().n_nodes() << " points: " << duration << " ms" << std::endl;
    outfile << del_10000.dcel().n_nodes() << "," << duration << "\n"; 

    /*Delaunay<2, 2> del_100000(boundary,100000);
    start = high_resolution_clock::now();  // starting measuring time 
    //Delaunay<2, 2> del_10000(boundary,10000);  
    del_100000.Ruppert_refinement(2.0);
    end = high_resolution_clock::now();    // ending measuring time 
    duration = duration_cast<milliseconds>(end - start).count();
    std::cout << "elapsed time for " << del_100000.dcel().n_nodes() << " points: " << duration << " ms" << std::endl;
    outfile << del_100000.dcel().n_nodes() << "," << duration << "\n"; */
/*
    Delaunay<2, 2> del_100000(boundary,100000);
    start = high_resolution_clock::now();  // starting measuring time 
    //Delaunay<2, 2> del_10000(boundary,10000);  
    del_100000.Ruppert_refinement(2.0);
    end = high_resolution_clock::now();    // ending measuring time 
    duration = duration_cast<milliseconds>(end - start).count();
    std::cout << "elapsed time for " << del_100000.dcel().n_nodes() << " points: " << duration << " ms" << std::endl;
    outfile << del_100000.dcel().n_nodes() << "," << duration << "\n"; 

    start = high_resolution_clock::now();  // starting measuring time 
    Delaunay<2, 2> del_100000(boundary,100000);  
    //del_100000.Ruppert_refinement(2.0);
    end = high_resolution_clock::now();    // ending measuring time 
    duration = duration_cast<milliseconds>(end - start).count();
    std::cout << "elapsed time for " << 100000 << " points: " << duration << " ms" << std::endl;
    outfile << 100000 << "," << duration << "\n"; */

    /*start = high_resolution_clock::now();  // starting measuring time 
    Delaunay<2, 2> delaunay_1000000(boundary,1000000);  
    end = high_resolution_clock::now();    // ending measuring time 
    duration = duration_cast<milliseconds>(end - start).count();
    std::cout << "elapsed time for " << 1000000 << " points: " << duration << " ms" << std::endl;
    outfile << 1000000 << "," << duration << "\n"; */

    /*auto start = high_resolution_clock::now();  // starting measuring time 
    Delaunay<2, 2> delaunay_8000(boundary,1000, std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2>> { hole2, hole3});  
    auto end = high_resolution_clock::now();    // ending measuring time 
    auto duration = duration_cast<milliseconds>(end - start).count();
    std::cout << "elapsed time for " << 8000 << " points: " << duration << " ms" << std::endl;
    outfile << 8000 << "," << duration << "\n"; 

    start = high_resolution_clock::now();  // starting measuring time 
    Delaunay<2, 2> delaunay_10000(boundary,5000, std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2>> { hole2, hole3});  
    end = high_resolution_clock::now();    // ending measuring time 
    duration = duration_cast<milliseconds>(end - start).count();
    std::cout << "elapsed time for " << 10000 << " points: " << duration << " ms" << std::endl;
    outfile << 10000 << "," << duration << "\n"; 

    start = high_resolution_clock::now();  // starting measuring time 
    Delaunay<2, 2> delaunay_20000(boundary,10000, std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2>> { hole2, hole3});  
    end = high_resolution_clock::now();    // ending measuring time 
    duration = duration_cast<milliseconds>(end - start).count();
    std::cout << "elapsed time for " << 20000 << " points: " << duration << " ms" << std::endl;
    outfile << 20000 << "," << duration << "\n"; 
    
    outfile.close();*/
    
    return 0;
}
