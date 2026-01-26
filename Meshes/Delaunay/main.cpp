#include <fstream>
#include <nlohmann/json.hpp> 
using json = nlohmann::json;
#include <geometry.h>
using namespace fdapde;
#include "domains.h"



using coords_t = typename Eigen::Matrix<double, 1, 2>;
using node_t = typename DCEL<2, 2>::node_t;


// function to obtain .mesh and .sol files to use in MMG adaptation
void write_mesh_and_metric(fdapde::DCEL<2,2>& dcel,const std::unordered_map<node_t*, Eigen::Matrix2d>& node_metrics,const std::string& mesh_filename,const std::string& sol_filename) {
    std::ofstream mesh_out(mesh_filename);
    std::ofstream sol_out(sol_filename);

    const std::size_t N = dcel.n_nodes();
    const std::size_t T = dcel.n_cells();

    std::unordered_map<node_t*, std::size_t> node_id;
    node_id.reserve(N);
    std::size_t id = 1;
    for (auto it = dcel.nodes_begin(); it != dcel.nodes_end(); ++it) {
        node_id[&(*it)] = id++;
    }

    mesh_out << "MeshVersionFormatted 2\n\n";
    mesh_out << "Dimension 2\n\n";
    mesh_out << "Vertices\n";
    mesh_out << N << "\n";

    for (auto it = dcel.nodes_begin(); it != dcel.nodes_end(); ++it) {
        const auto& coords = it->coords();
        mesh_out << coords(0) << " " << coords(1) << " 0\n";
    }

    mesh_out << "\nTriangles\n";
    mesh_out << T << "\n";

    for (auto cit = dcel.cells_begin(); cit != dcel.cells_end(); ++cit) {
        const auto* c = &(*cit);
        auto h = c->halfedge();
        node_t* n0 = h->node();
        node_t* n1 = h->next()->node();
        node_t* n2 = h->next()->next()->node();
        mesh_out << node_id.at(n0) << " " << node_id.at(n1) << " " << node_id.at(n2) << " 0\n";
    }

    mesh_out << "\nEnd\n";

    sol_out << "MeshVersionFormatted 2\n\n";
    sol_out << "Dimension 2\n\n";
    sol_out << "SolAtVertices\n";
    sol_out << N << "\n";
    sol_out << "1 3\n"; 

    for (auto it = dcel.nodes_begin(); it != dcel.nodes_end(); ++it) {
        node_t* node = &(*it);
        auto mit = node_metrics.find(node);
        const Eigen::Matrix2d& M = mit->second;
        sol_out << " " << M(0,0) << " " << M(0,1) << " " << M(1,1) << "\n";
    }

    sol_out << "\nEnd\n";
}


struct TensorDecomposition {
    Eigen::Vector2d ev;       // ev(0) = lambda_min, ev(1) = lambda_max
    Eigen::Matrix2d R;        // eigenvectors matrix
};
TensorDecomposition decomposeTensor(const Eigen::Matrix2d& M) {
    Eigen::Matrix2d Ms = 0.5 * (M + M.transpose()); 
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es(Ms);
    
    TensorDecomposition result;
    result.ev = es.eigenvalues();
    result.R = es.eigenvectors();
    
    // order ev(0) = min, ev(1) = max 
    if (result.ev(0) > result.ev(1)) {
        std::swap(result.ev(0), result.ev(1));
        result.R.col(0).swap(result.R.col(1));
    }
    return result;
}

// build the scaled metric at a node given its tensor decomposition
Eigen::Matrix2d build_metric_scaled(const TensorDecomposition& decomposition,double R_max, double Beta, double hmax_tilde_sq  ) {
    double ev0 = decomposition.ev(0);
    double ev1 = decomposition.ev(1);

    // alpha-->density, gamma-->shape
    double alpha_global = std::sqrt(ev0 * ev1);
    double gamma_global = ev0 / alpha_global;   
    double gammabeta = std::pow(gamma_global, Beta); 
    double sqlaK = alpha_global * gammabeta;   
    
    // inverse final eigenvectors (1/h^2)
    double invsqlambda1 = 1.0 / (hmax_tilde_sq * alpha_global * gammabeta);   
    double invsqlambda2 = gammabeta / (hmax_tilde_sq * alpha_global); 

    // anisotropy constraint
    double ratio = invsqlambda2 / invsqlambda1;
    if (ratio > R_max) { 
        invsqlambda1 = invsqlambda2 / (R_max * R_max);
    }  
    
    // final reconstruction of M = R * diag(invsqlambda1, invsqlambda2) * R^T
    Eigen::Matrix2d M_final = decomposition.R.col(0) * decomposition.R.col(0).transpose() * invsqlambda1 + decomposition.R.col(1) * decomposition.R.col(1).transpose() * invsqlambda2;

    return M_final;
}

// function that normalizes the metric on dcel vertices given a raw metric function
std::unordered_map<node_t*, Eigen::Matrix<double, 2, 2>> normalize_metric_fun(std::function<Eigen::Matrix2d(const coords_t&)> Mraw, fdapde::DCEL<2, 2>& dcel_, const double H_MAX_USER  = 0.4,const double BETA_POWER  = 2.0, const double R_MAX_ANISOTROPY = 2.5) {
 
    std::vector<TensorDecomposition> Decomposed_nodes(dcel_.n_nodes());
    std::unordered_map<node_t*, Eigen::Matrix<double, 2, 2>> node_metrics;  

    auto it = dcel_.nodes_begin();
    for (size_t i = 0; i < dcel_.n_nodes(); ++i, ++it) {
        Eigen::Matrix2d Mraw_i = Mraw(it->coords()); 
        Decomposed_nodes[i] = decomposeTensor(Mraw_i);           
    }

    double sqlaK_max = 0.0;
    for (const auto& decomp : Decomposed_nodes) {
        double ev0 = decomp.ev(0); // lambda_min
        double ev1 = decomp.ev(1); // lambda_max
        
        double alphaV_i = std::sqrt(ev0 * ev1);
        double gammaV_i = ev0 / alphaV_i;
        
        double sqlaK_i = alphaV_i * std::pow(gammaV_i, BETA_POWER);
        if (sqlaK_i > sqlaK_max) {
            sqlaK_max = sqlaK_i;
        }
    }

    double hmax_tilde_sq = (H_MAX_USER * H_MAX_USER) / sqlaK_max;

    it = dcel_.nodes_begin();
    for (size_t i = 0; i < dcel_.n_nodes(); ++i, ++it) {
        const auto& decomp = Decomposed_nodes[i];
        node_metrics[&(*std::next(dcel_.nodes_begin(), i))] = build_metric_scaled(decomp,R_MAX_ANISOTROPY,BETA_POWER,hmax_tilde_sq);
    }

    return node_metrics;
}

std::unordered_map<node_t*, Eigen::Matrix<double, 2, 2>> normalize_matrix(std::unordered_map<DCEL<2,2>::node_t*,Eigen::Matrix2d> map, fdapde::DCEL<2, 2>& dcel_, const double HMAX  = 0.4,const double C_DENSITY = 1.0,    const double BETA  = 2.0, const double R_MAX = 2.5) {
    
    const size_t N = dcel_.n_nodes();
    std::vector<TensorDecomposition> dec(N);

    int i = 0;
    for (auto it = dcel_.nodes_begin(); it != dcel_.nodes_end(); ++it, ++i)
        dec[i] = decomposeTensor(map.at(&(*it)));

    double sqlaK_max = 0.0;
    for (auto& d : dec) {
        double ev0 = d.ev(0), ev1 = d.ev(1);
        double alpha = std::sqrt(ev0*ev1);
        double gamma = ev0 / alpha;
        sqlaK_max = std::max(sqlaK_max, alpha * std::pow(gamma, BETA));
    }
    const double hmax_tilde_sq = (HMAX * HMAX) / sqlaK_max;

    i = 0;
    for (auto it = dcel_.nodes_begin(); it != dcel_.nodes_end(); ++it, ++i){
        map[&(*it)] = C_DENSITY* build_metric_scaled(dec[i], R_MAX, BETA, hmax_tilde_sq);
        //std::cout << "Node " << it->id() << ": Metric after normalization:\n" << map[&(*it)] << "\n";
    }

    return map;
}







int main() {

    //-----------------------CALL TO BUILD THE MESH AND PRINT ITS QUALITY MEASURES-------------------------------------

    std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2>> domain={skyline, building1, building2, building3};  //insert the domain you want to triangulate, with its eventual subregions
    

    auto random_points_in_rectangle = [](int n, double xmin, double xmax, double ymin, double ymax, unsigned long long seed)-> Eigen::Matrix<double, Eigen::Dynamic, 2>
    {
        Eigen::Matrix<double, Eigen::Dynamic, 2> P(n,2);
        std::mt19937_64 rng(seed);
        std::uniform_real_distribution<double> Ux(xmin, xmax);
        std::uniform_real_distribution<double> Uy(ymin, ymax);
        for (std::size_t i=0; i<n; ++i) {
            P(i,0) = Ux(rng);
            P(i,1) = Uy(rng);
        }
        return P;
    };

    const auto random_points_around_segment = [](int n,const Eigen::Vector2d& A,const Eigen::Vector2d& B,double sigma_perp,double sigma_tang = 0.0,unsigned long long seed = 42ULL)-> Eigen::Matrix<double, Eigen::Dynamic, 2>
    {
        Eigen::Matrix<double, Eigen::Dynamic, 2> P(n, 2);
        std::mt19937_64 rng(seed);
        std::uniform_real_distribution<double> Ut(0.0, 1.0);
        std::normal_distribution<double> Z(0.0, 1.0);

        Eigen::Vector2d d = B - A;
        double L = d.norm();
        Eigen::Vector2d t_hat, n_hat;
        if (L > 0.0) {
            t_hat = d / L;
            n_hat = Eigen::Vector2d(-t_hat.y(), t_hat.x()); // rotazione +90°
        } else {
            t_hat = Eigen::Vector2d::UnitX();
            n_hat = Eigen::Vector2d::UnitY();
        }

        for (int i = 0; i < n; ++i) {
            double t = Ut(rng);
            Eigen::Vector2d base = A + t * d;
            double z_perp = Z(rng);
            double z_tang = (sigma_tang > 0.0 ? Z(rng) : 0.0);
            Eigen::Vector2d p = base + sigma_perp * z_perp * n_hat + sigma_tang * z_tang * t_hat;
            P.row(i) = p;
        }
        return P;
    }; 

    auto export_data_to_txt = [](const Eigen::Matrix<double, Eigen::Dynamic, 2>& data, const std::string& filename) {
        std::ofstream f(filename);
        if (!f) throw std::runtime_error("Cannot open " + filename);
        for (int i = 0; i < data.rows(); ++i)
            f << data(i,0) << " " << data(i,1) << "\n";
    };
    


    using coords_t = Eigen::Matrix<double, 1, 2>;
    auto metric_fun = [&](coords_t coords, double alpha = 1) -> Eigen::Matrix<double, 2, 2> {
        /*const double gamma = std::cos(4.0 * M_PI * (coords(1) - coords(0))) + 2.0;
        const double nux = 1.0 / std::sqrt(2.0);
        const double nuy = 1.0 / std::sqrt(2.0);
        const double d11 = alpha * (nux*nux*gamma + nuy*nuy / gamma);
        const double d12 = alpha * (nux*nuy*(gamma - 1.0 / gamma));
        const double d22 = alpha * (nuy*nuy*gamma + nux*nux / gamma);*/   //LINEARE
        const double gamma = std::cos(4 * M_PI * (coords(1) - coords(0)*coords(0))) + 2;
        const double nux   = 1. / sqrt(4 * coords(0)*coords(0) + 1);
        const double nuy   = 2 * coords(0) / sqrt(4 * coords(0)*coords(0) + 1);
        const double d11 = alpha * (nux*nux * gamma + nuy*nuy / gamma);
        const double d12 = alpha * (nux * nuy * (gamma - 1. / gamma));
        const double d22 = alpha * (nuy*nuy * gamma + nux*nux / gamma);
        Eigen::Matrix<double, 2, 2> M;
        M << d11, d12, d12, d22;
        return M;
    };

    auto metric_fun_isotropa = [&](coords_t coords) -> Eigen::Matrix<double, 2, 2> {
        double x = coords(0);
        double y = coords(1);
        double eps = 1e-6;
        Eigen::Matrix<double, 2, 2> M;
        M << x*x+eps, 0.0, 0.0, x*x+eps;  
        return M;
    };

    auto metric_fun_anisotropa = [&](coords_t coords) -> Eigen::Matrix<double, 2, 2> {
        double r = coords.norm();
        double h_min = 0.0003 + 0.0002*r;
        double h_max = 4 + 2*r;
        double lambda_max = 1.0/(h_min*h_min);
        double lambda_min = 1.0/(h_max*h_max);
        coords_t n; n << 1/std::sqrt(2.0), 1/std::sqrt(2.0);
        coords_t nperp; nperp << -n(1), n(0);
        Eigen::Matrix2d M = lambda_max * (n.transpose() * n) + lambda_min * (nperp.transpose() * nperp);
        return M;
    };




    //auto data1 = random_points_around_segment(300, Eigen::Vector2d(0.1, 0.1), Eigen::Vector2d(0.9, 0.9), 0.05, 0.05, 12345);
    /*auto data1 = random_points_in_rectangle(250, 0, 0.6, 0, 1, 67890);
    auto data2 = random_points_in_rectangle(250, 0, 1, 0.5, 1, 67890);
    Eigen::MatrixXd data(data1.rows() + data2.rows(), data1.cols());
    data << data1, data2;*/
    auto data = random_points_in_rectangle(50, 0, 0.2, 0.7, 1, 67890);
    export_data_to_txt(data, "Meshes/Delaunay/data_points.txt");

    
    Delaunay<2, 2> del(std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2>>{U}, 0);
    del.refinement(30, del.domain_area()/100);
    del.dcel().export_to_json("Meshes/Delaunay/delaunay_output.json");

    //Adaptivity<2,2, fdapde::AdaptiveStrategy::GradientMagnitude> adapt(del,data, 20, del.domain_area()/100, 1e-1);

    Delaunay<2,2>& delaunay = del;
    double prev_diameter = delaunay.min_triangle_diameter();
    double tol = 5e-3;

    for(int i=0; i<100; ++i){
        std::cout << "Iter " << i+1 << " of adaptivity." << std::endl;
        std::cout << "Current number of cells: " << delaunay.dcel().n_cells() << std::endl;

        NodeMetric<2,2,AdaptiveStrategy::NodeDensity> node_metrics_obj(delaunay.dcel(), data, del.domain_area()/100);
        auto node_metrics = node_metrics_obj.node_density_knn_geodesic(); 
        node_metrics = normalize_matrix(node_metrics, delaunay.dcel(), 0.05, 0.1);  //hmax,c

        //auto node_metrics = normalize_metric_fun(metric_fun, delaunay.dcel(), 0.075, 2.5, 2.5);
        //write_mesh_and_metric(delaunay.dcel(), node_metrics, "Meshes/Delaunay/mesh_init100.mesh", "Meshes/Delaunay/metric100.sol");

        Adaptivity<2,2> adapt(delaunay, node_metrics);
        adapt.adaptivity_cycle(20, del.domain_area());
        adapt.dcel().export_to_json("Meshes/Delaunay/delaunay_output.json");  //()"Meshes/Delaunay/c_var/datapoints_" + std::to_string(c) + ".json");
        std::cout << adapt.dcel().n_cells() << " cells after adaptivity." << std::endl;

        // stopping criterion
        double curr_diameter = adapt.delaunay().min_triangle_diameter();
        std::cout << "curr: " << curr_diameter << "   prev: " << prev_diameter << std::endl;
        // errore tra ultima iter di adapt e iter i+1-esima è perchè i valori interpolati dentro al ciclo sono diversi dai valori reali della funzione
        if(std::abs(curr_diameter - prev_diameter)/prev_diameter < tol){
            std::cout << "Converged: relative change in min triangle diameter is below threshold. Iter " << i+1 << std::endl;
            break;
        }
        
        prev_diameter = curr_diameter;
        delaunay = adapt.delaunay();
        
    }

    




    /*int N = 360;             
    double r = 1.0;            
    double eps = 1;      // absolute tolerance     
    double eps_rel = 0.9;  // relative tolerance based on domain metric (choose the minimum value between eps, eps_rel*segm_len/2.0, segm_len/2.0)   
    Eigen::Matrix<double, Eigen::Dynamic, 2> boundary(N, 2);
    for (int i = 0; i < N; ++i) {
        double theta = 2.0 * M_PI * i / N;
        boundary(i, 0) = r * std::cos(theta);
        boundary(i, 1) = r * std::sin(theta);
    }
    std::cout << "Original boundary has " << boundary.rows() << " points." << std::endl;
    DCEL<2, 2> dcel = DCEL<2, 2>::make_polygon(boundary);
    dcel.export_to_json("Meshes/Delaunay/dcel.json");
    // Applica la semplificazione
    Eigen::Matrix<double, Eigen::Dynamic, 2> simplified = simplify_boundary(boundary, eps, eps_rel);
    std::cout << "Simplified boundary has " << simplified.rows() << " points." << std::endl;
    DCEL<2, 2> dcel_new = DCEL<2, 2>::make_polygon(simplified);
    dcel_new.export_to_json("Meshes/Delaunay/dcel_2.json");*/

    /*DCEL<2, 2> dcel1 = DCEL<2, 2>::make_polygon(skyline);
    dcel1.export_to_json("Meshes/Delaunay/dcel_2.json");
    auto boundary = convex_offset(skyline, 2.0);
    DCEL<2, 2> dcel = DCEL<2, 2>::make_polygon(boundary);
    dcel.export_to_json("Meshes/Delaunay/dcel.json");*/

    return 0;
}