
#include <fstream>
#include <nlohmann/json.hpp> 
using json = nlohmann::json;
#include <geometry.h>
using namespace fdapde;
#include "domains.h"


int main() {

//-----------------------CALL TO BUILD THE MESH AND PRINT ITS QUALITY MEASURES-------------------------------------
Delaunay<2, 2> del(std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2>> {skyline},20,1.0,0, holes_skyline);
del.dcel().export_to_json("Tests_mesh_efficiency/dcel_output.json");
//del.print_statistics();

    
//-----------------------TEST OF COMPUTATIONAL EFFICIENCY OF CONFLICT GRAPH ALGORITHM------------------------------


/*    std::ofstream outfile("Tests_mesh_efficiency/timing.csv");
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


    std::ofstream outfile("Tests_mesh_efficiency/timing.csv");
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

    //-----------------------END TEST OF COMPUTATIONAL EFFICIENCY OF REFINEMENT ALGORITHM------------------------------
    
    return 0;
}