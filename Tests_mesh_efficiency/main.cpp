
#include <fstream>
#include <nlohmann/json.hpp> 
using json = nlohmann::json;
#include <geometry.h>
using namespace fdapde;
#include "domains.h"


int main() {

//-----------------------CALL TO BUILD THE MESH AND PRINT ITS QUALITY MEASURES-------------------------------------
Eigen::Matrix<double, 2, 2> internal;
internal << 30.5,0.5,
            24.5,0.5;
Eigen::Matrix<double, 4, 2> region1;
region1 << 
    0.0, 0.0,
    2000.0, 0.0,
    2000.0, 1000.0,
    0.0, 1000.0;
Eigen::Matrix<double, 4, 2> region2;
region2 << 
    2000.0, 0.0,
    4000.0, 0.0,
    4000.0, 1000.0,
    2000.0, 1000.0;
Eigen::Matrix<double, 4, 2> region3;
region3 << 
    2000.0, 1000.0,
    4000.0, 1000.0,
    4000.0, 2000.0,
    2000.0, 2000.0;
Eigen::Matrix<double, 4, 2> region4;
region4 << 
    0.0, 1000.0,
    2000.0, 1000.0,
    2000.0, 2000.0,
    0.0, 2000.0;

std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2>>> holes(5); 

// === Regione 1 ===
Eigen::Matrix<double, 4, 2> hole1_r1;
hole1_r1 << 300.0, 200.0,  500.0, 200.0,  500.0, 400.0,  300.0, 400.0;

Eigen::Matrix<double, 4, 2> hole2_r1;
hole2_r1 << 1500.0, 600.0,  1700.0, 600.0,  1700.0, 800.0,  1500.0, 800.0;

holes[0].push_back(hole1_r1);
holes[0].push_back(hole2_r1);
holes[1].push_back(hole1_r1);
holes[1].push_back(hole2_r1);

// === Regione 2 ===
Eigen::Matrix<double, 4, 2> hole1_r2;
hole1_r2 << 2300.0, 300.0,  2500.0, 300.0,  2500.0, 500.0,  2300.0, 500.0;

Eigen::Matrix<double, 4, 2> hole2_r2;
hole2_r2 << 3500.0, 100.0,  3700.0, 100.0,  3700.0, 300.0,  3500.0, 300.0;

holes[0].push_back(hole1_r2);
holes[0].push_back(hole2_r2);
holes[3].push_back(hole1_r2);
holes[3].push_back(hole2_r2);

// === Regione 3 ===
Eigen::Matrix<double, 4, 2> hole1_r3;
hole1_r3 << 2200.0, 1200.0,  2400.0, 1200.0,  2400.0, 1400.0,  2200.0, 1400.0;

Eigen::Matrix<double, 4, 2> hole2_r3;
hole2_r3 << 3200.0, 1700.0,  3400.0, 1700.0,  3400.0, 1900.0,  3200.0, 1900.0;

holes[0].push_back(hole1_r3);
holes[0].push_back(hole2_r3);
holes[2].push_back(hole1_r3);
holes[2].push_back(hole2_r3);

// === Regione 4 ===
Eigen::Matrix<double, 4, 2> hole1_r4;
hole1_r4 << 200.0, 1400.0,  400.0, 1400.0,  400.0, 1600.0,  200.0, 1600.0;

Eigen::Matrix<double, 4, 2> hole2_r4;
hole2_r4 << 1400.0, 1100.0,  1600.0, 1100.0,  1600.0, 1300.0,  1400.0, 1300.0;

holes[0].push_back(hole1_r4);
holes[0].push_back(hole2_r4);
holes[4].push_back(hole1_r4);
holes[4].push_back(hole2_r4);



Delaunay<2, 2> del(std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2>> {rectangle, region1, region3, region4, region2},0);
std::cout << typeid(del).name() << std::endl;
del.dcel().export_to_json("Tests_mesh_efficiency/dcel_output.json");
del.print_statistics();

    
//-----------------------TEST OF COMPUTATIONAL EFFICIENCY OF CONFLICT GRAPH ALGORITHM------------------------------


    /*std::ofstream outfile("Tests_mesh_efficiency/timing.csv");
    outfile << "NumPoints,TimeElapsed(ms)\n";

    auto start = high_resolution_clock::now(); 
    Delaunay<2, 2> del_10000(std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2>> {rectangle}, 10000);
    auto end = high_resolution_clock::now();  

    auto duration = duration_cast<milliseconds>(end - start).count();
    std::cout << "elapsed time for 10000 internal points: " << duration << " ms" << std::endl;
    outfile << 10000 << "," << duration << "\n"; 

    start = high_resolution_clock::now(); 
    Delaunay<2, 2> del_100000(std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2>> {rectangle}, 100000);  
    end = high_resolution_clock::now(); 

    duration = duration_cast<milliseconds>(end - start).count();
    std::cout << "elapsed time for 100000 internal points: " << duration << " ms" << std::endl;
    outfile << 100000 << "," << duration << "\n"; 

    start = high_resolution_clock::now(); 
    Delaunay<2, 2> del_1000000(std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2>> {rectangle}, 1000000);  
    end = high_resolution_clock::now(); 

    duration = duration_cast<milliseconds>(end - start).count();
    std::cout << "elapsed time for 1000000 internal points: " << duration << " ms" << std::endl;
    outfile << 1000000 << "," << duration << "\n"; 
    
    outfile.close(); */

//-----------------------END TEST OF COMPUTATIONAL EFFICIENCY OF CONFLICT GRAPH ALGORITHM------------------------------


//-----------------------TEST OF COMPUTATIONAL EFFICIENCY OF REFINEMENT ALGORITHM------------------------------


    /*std::ofstream outfile("Tests_mesh_efficiency/timing.csv");
    outfile << "NumPoints,TimeElapsed(ms)\n";
    
    Delaunay<2, 2> del1(std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2>> {rectangle});
    auto start = high_resolution_clock::now();    
    del1.Ruppert_refinement(20, del1.domain_area()/1000.0);
    auto end = high_resolution_clock::now();  

    auto duration = duration_cast<milliseconds>(end - start).count();
    std::cout << "elapsed time for "<<del1.dcel().n_nodes() <<" points: " << duration << " ms" << std::endl; 
    outfile << del1.dcel().n_nodes() << "," << duration << "\n"; 

    Delaunay<2, 2> del2(std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2>> {rectangle});
    start = high_resolution_clock::now();    
    del2.Ruppert_refinement(20, del2.domain_area()/10000.0);
    end = high_resolution_clock::now();  

    duration = duration_cast<milliseconds>(end - start).count();
    std::cout << "elapsed time for "<<del2.dcel().n_nodes() <<" points: " << duration << " ms" << std::endl; 
    outfile << del2.dcel().n_nodes() << "," << duration << "\n"; 

    Delaunay<2, 2> del3(std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2>> {rectangle});
    start = high_resolution_clock::now();    
    del3.Ruppert_refinement(20, del3.domain_area()/50000.0);
    end = high_resolution_clock::now();  

    duration = duration_cast<milliseconds>(end - start).count();
    std::cout << "elapsed time for "<<del3.dcel().n_nodes() <<" points: " << duration << " ms" << std::endl; 
    outfile << del3.dcel().n_nodes() << "," << duration << "\n"; 

    Delaunay<2, 2> del4(std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2>> {rectangle});
    start = high_resolution_clock::now();    
    del4.Ruppert_refinement(20, del4.domain_area()/100000.0);
    end = high_resolution_clock::now();  

    duration = duration_cast<milliseconds>(end - start).count();
    std::cout << "elapsed time for "<<del4.dcel().n_nodes() <<" points: " << duration << " ms" << std::endl; 
    outfile << del4.dcel().n_nodes() << "," << duration << "\n"; 

    outfile.close();
    */
    //-----------------------END TEST OF COMPUTATIONAL EFFICIENCY OF REFINEMENT ALGORITHM------------------------------
    
    return 0;
}