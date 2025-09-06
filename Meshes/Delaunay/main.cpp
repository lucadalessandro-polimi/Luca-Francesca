#include <fstream>
#include <nlohmann/json.hpp> 
using json = nlohmann::json;
#include <geometry.h>
using namespace fdapde;
#include "domains.h"


int main() {

    //-----------------------CALL TO BUILD THE MESH AND PRINT ITS QUALITY MEASURES-------------------------------------

    std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2>> domain={skyline, building1, building2, building3};  //insert the domain you want to triangulate, with its eventual subregions
    
    Delaunay<2, 2> del({skyline}, 0); 
    std::cout << "Area senza offset: " << del.domain_area() << std::endl;

    auto dom = fdapde::offset_polygon(skyline, -0.008);
    Delaunay<2, 2> del2({dom}, 0);
    std::cout << "Area con offset: " << del2.domain_area() << std::endl;     
    del2.dcel().export_to_json("Meshes/Delaunay/delaunay_output.json");
    //del.print_statistics();

    return 0;
}