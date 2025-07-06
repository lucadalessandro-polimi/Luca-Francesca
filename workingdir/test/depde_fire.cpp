#include <fdaPDE/fdapde.h>
#include <unsupported/Eigen/SparseExtra>

using namespace fdapde;

int main(){
    
    std::string datadir = "../data/fire-data-south/";
   
    // read shape file -----------------------------------------------------------------------------
    std::string italy_str = datadir + "south-italy-regions-no-sic.shp";
    
    // solo bordo esterno senza vincoli su regioni
    std::string italy_ext = datadir + "south-italy-no-sic.shp";
    
    auto italy_shp = read_shp(italy_str);
    auto italy_ext_shp = read_shp(italy_ext);
    
    std::cout << italy_shp << std::endl;

    std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2>> coords(1+italy_shp.n_records());    
    
    coords[0] = italy_ext_shp.polygon(0).nodes()[0];

    std::cout << "\nREGIONS " << std::endl;
    for(int k = 1; k < italy_shp.n_records() +1 ; k++){
        
        coords[k] = italy_shp.polygon(k-1).nodes()[0];
        std::cout << "rows: " << coords[k].rows() << " cols: "<< coords[k].cols() << std::endl;
        std::cout << "first pt: " << coords[k](0,0) << " " << coords[k](0,1) << std::endl;
        std::cout << "last pt: " << coords[k](coords[k].rows()-1,0) << " " 
                                 << coords[k](coords[k].rows()-1,1) << "\n" << std::endl;
        std::cout << std::endl;
    }
    
    
    Delaunay<2, 2> del(coords, 0);
    del.refinement(20., del.domain_area()/700);
    Triangulation<2,2> D=del.triangulation();

    
    // read data -----------------------------------------------------------------------------------
    std::string points_str = datadir + "fire-data.shp";
    auto points_shp = read_shp(points_str);
    
    std::cout << "\n FIRE LOCATIONS" << std::endl;
    std::cout << points_shp <<std::endl;
    Eigen::MatrixXd locations = Eigen::MatrixXd::Zero(points_shp.n_records(), 2);
    
    for(int i=0; i < points_shp.n_records(); ++i){ 
        locations(i,0) = points_shp.point(i).x; 
        locations(i,1) = points_shp.point(i).y;
    }
    
  
    GeoFrame data(D);
    auto& l = data.insert_scalar_layer<POINT>("l", locations);
    std::cout << l << std::endl;
    // DEPDE ---------------------------------------------------------------------------------------
    
    // physics
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction  v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);
    // modeling
    DEPDE m(data, fe_de_elliptic(a, F));
    
    // Uniform init 
    Eigen::MatrixXd g_init = Eigen::MatrixXd::Ones(D.nodes().rows(),1);
    g_init /= D.measure();
    g_init.array() = g_init.array().log();
    
    std::cout <<"measure: " << D.measure() <<"\n" << std::endl;
    
    m.set_llik_tolerance(1e-4);
    m.fit(0.1, g_init, GradientDescent<Dynamic, BacktrackingLineSearch> {500, 1e-4, 1e-2});
    
    Eigen::saveMarket(D.nodes(), datadir + "dofs.mtx");
    Eigen::saveMarket(m.density(), datadir + "coeffs.mtx");
    Eigen::saveMarket(D.cells(), datadir + "cells.mtx");
   
    return 0;
}
