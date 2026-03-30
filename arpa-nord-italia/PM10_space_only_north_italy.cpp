#include<fdaPDE/fdapde.h>
#include<string>
#include<unsupported/Eigen/SparseExtra>
#include<utils.h>
using namespace fdapde;

int main(int argc, char *argv[]){

    using vector_t = Eigen::Matrix<double, Dynamic, 1>; 
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>; 
    
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    using sparsesolver_t = Eigen::SparseLU<sparse_matrix_t>;
    constexpr int local_dim = 2;
    using PointT = Eigen::Matrix<double, local_dim, 1>;
    
    std::string mesh_path = std::string( argv[1] ) + "/"; //"mesh_1/"; // PASSARE QUESTO COME ARGV
    std::string data_path = "data/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "cells.csv", mesh_path + "neigh.csv", true, true);

    std::string train_locs_path = argv[2];
    std::string train_obs_path  = argv[3];
    std::string test_locs_path  = argv[4];
    std::string pred_out_path   = argv[5];
        
    // data
    //matrix_t locs = read_csv<double>(data_path + "2022.10.789_space_locs.csv").as_matrix();
    matrix_t locs = read_csv<double>(train_locs_path).as_matrix();
    GeoFrame data(D);
    auto& layer = data.insert_scalar_layer<POINT>("layer", locs);
    //layer.load_csv<double>(data_path + "2022.10.789_PM10.csv");   
    layer.load_csv<double>(train_obs_path);
    std::cout << layer << std::endl;


    FeSpace Vh(D, P1<1>);   // linear finite element in space
    TrialFunction f(Vh);
    TestFunction  v(Vh);
    auto a_D = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u_D;
    auto F_D = integral(D)(u_D * v);

    double alpha = 0.95;
    QSRPDE model("y ~ f", data, alpha, fe_ls_elliptic(a_D, F_D));
    std::cout << "model created" << std::endl;
    vector_t exps = vector_t::LinSpaced(100,-8,0);
    int n_lambda = exps.size();
    vector_t lambdas = vector_t::Ones(n_lambda);
    for(int i=0; i<n_lambda;++i) lambdas[i] = std::pow(10, exps[i]);

    GridSearch<1> opt;

    // GCV
    auto start = std::chrono::steady_clock::now();
    opt.optimize(model.gcv(200, 1234), lambdas);
    //model.fit(alpha, 1e-3);  
    auto end = std::chrono::steady_clock::now();
    auto seconds =
        std::chrono::duration_cast<std::chrono::seconds>(end - start);

    std::cout << "elapsed time: " << seconds.count() << " s" << std::endl;
    
    // salvare risultati in mesh path ?
    std::cout << "\n lambda opt = argmin(gcv) = " << opt.optimum() << " , min(gcv) = "  << opt.value() << std::endl;

    std::cout << "gcv values: " << opt.values() << std::endl; // opt.values() è uno std::vector! operator<< in include/utils.h

    vector_t values = vector_t::Zero(opt.values().size());
    for(int i = 0; i < opt.values().size(); ++i) values[i] = opt.values()[i];

    //Eigen::saveMarket(values, mesh_path + "gcv_values.mtx");
    //Eigen::saveMarket(lambdas, mesh_path + "lambdas.mtx");
    //Eigen::saveMarket(opt.optimum(), mesh_path + "lambda_opt.mtx");
    
    model.fit(alpha, opt.optimum()[0]); // è un eigen vector di dim 1, prende double

    
    //Eigen::saveMarket(model.f(), mesh_path + "f.mtx");
    //Eigen::saveMarket(model.misfit(), mesh_path + "g.mtx");

    vector_t coeffs = model.f();
    FeFunction f_est(Vh, coeffs);
    matrix_t locs_test = read_csv<double>(test_locs_path).as_matrix();
    vector_t pred_test(locs_test.rows());

    for(int i = 0; i < locs_test.rows(); ++i){
        PointT p = locs_test.row(i);
        pred_test[i] = f_est(p);
    }
    Eigen::saveMarket(pred_test, pred_out_path);


    return 0;
}
