#include <fstream>
#include <nlohmann/json.hpp> 
using json = nlohmann::json;
#include <geometry.h>
using namespace fdapde;
#include "domains.h"


int main() {

    //-----------------------CALL TO BUILD THE MESH AND PRINT ITS QUALITY MEASURES-------------------------------------

    Delaunay<2, 2> del(std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2>> {skyline, building1, building2, building3},20,1,0, holes_skyline);
    del.dcel().export_to_json("Meshes/Delaunay/delaunay_output.json");
    del.print_statistics();

    return 0;
}