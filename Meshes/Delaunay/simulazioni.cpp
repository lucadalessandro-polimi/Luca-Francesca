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


Eigen::Matrix<double, Eigen::Dynamic, 2> read_locs_mtx(const std::string& path){
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open file: " + path);

    std::string line;

    // 1) Read header line
    if (!std::getline(in, line))
        throw std::runtime_error("Empty file: " + path);

    // 2) Skip comment lines starting with '%'
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] != '%') break;
    }
    if (line.empty())
        throw std::runtime_error("Missing size line in: " + path);

    // 3) Read size line: rows cols nnz
    int rows = 0, cols = 0;
    long long nnz = 0;
    {
        std::istringstream iss(line);
        if (!(iss >> rows >> cols >> nnz))
            throw std::runtime_error("Invalid size line: " + line);
    }

    if (cols != 2)
        throw std::runtime_error("Expected 2 columns (x,y). Found cols=" + std::to_string(cols));
    if (rows <= 0)
        throw std::runtime_error("Invalid rows=" + std::to_string(rows));

    // We'll build a dense matrix rows x 2
    Eigen::Matrix<double, Eigen::Dynamic, 2> locs(rows, 2);
    locs.setZero();

    // 4) Read triplets: i j value  (MatrixMarket uses 1-based indices)
    for (long long k = 0; k < nnz; ++k) {
        int i = 0, j = 0;
        double v = 0.0;
        if (!(in >> i >> j >> v))
            throw std::runtime_error("Unexpected EOF while reading triplets at k=" + std::to_string(k));

        if (i < 1 || i > rows || j < 1 || j > cols)
            throw std::runtime_error("Out-of-bounds triplet: (" + std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(v) + ")");

        locs(i - 1, j - 1) = v; // convert to 0-based
    }

    return locs;
}


// 2D Gaussian density with mean (mx,my) and covariance
// Sigma = [s11 s12; s12 s22]
double gauss2(double x, double y,double mx, double my, double s11, double s12, double s22)
{
    const double detS = s11 * s22 - s12 * s12;
    if (detS <= 0.0) {
        throw std::runtime_error("Covariance matrix is not positive definite (det <= 0).");
    }

    // Inverse of 2x2 symmetric matrix
    const double inv11 =  s22 / detS;
    const double inv12 = -s12 / detS;
    const double inv22 =  s11 / detS;

    const double dx = x - mx;
    const double dy = y - my;

    // Quadratic form: [dx dy] * invS * [dx; dy]
    const double q = inv11 * dx * dx + 2.0 * inv12 * dx * dy + inv22 * dy * dy;

    const double norm = 1.0 / (2.0 * M_PI * std::sqrt(detS));
    return norm * std::exp(-0.5 * q);
}

// mixture of 4 bivariate Gaussians
double mix4_gauss(coords_t coords)
{
    double x = coords(0);
    double y = coords(1);

    // weights
    constexpr double w1 = 0.25, w2 = 0.25, w3 = 0.25, w4 = 0.25;

    // centers (means)
    constexpr double mx1 = 0.3, my1 = 0.3;
    constexpr double mx2 = 0.3, my2 = 0.7;
    constexpr double mx3 = 0.7, my3 = 0.3;
    constexpr double mx4 = 0.7, my4 = 0.7;

    // covariance entries as in your FreeFem script:
    // Sigma_k = [s11k s12k; s12k s22k]
    constexpr double s111 = 0.0064, s121 = 0.0,   s221 = 0.0064;
    constexpr double s112 = 0.0036, s122 = 0.0,   s222 = 0.0100;
    constexpr double s113 = 0.01,   s123 = 0.0,   s223 = 0.0036;
    constexpr double s114 = 0.01,   s124 = 0.009, s224 = 0.01;

    // weighted sum
    return w1 * gauss2(x, y, mx1, my1, s111, s121, s221) + w2 * gauss2(x, y, mx2, my2, s112, s122, s222)
         + w3 * gauss2(x, y, mx3, my3, s113, s123, s223) + w4 * gauss2(x, y, mx4, my4, s114, s124, s224);
}
// Hessian of 2D Gaussian
inline Eigen::Matrix2d hessian_gauss2(double x, double y,double mx, double my,double s11, double s12, double s22){
    const double detS = s11*s22 - s12*s12;

    const double inv11 =  s22 / detS;
    const double inv12 = -s12 / detS;
    const double inv22 =  s11 / detS;

    const double dx = x - mx;
    const double dy = y - my;

    const double q = inv11*dx*dx + 2.0*inv12*dx*dy + inv22*dy*dy;

    const double norm = 1.0 / (2.0 * M_PI * std::sqrt(detS));
    const double g = norm * std::exp(-0.5 * q);

    const double a = inv11*dx + inv12*dy;
    const double b = inv12*dx + inv22*dy;
    // Hessian
    Eigen::Matrix2d H;
    H(0,0) = g * (a*a - inv11);   
    H(1,1) = g * (b*b - inv22);   
    H(0,1) = g * (a*b - inv12);   
    H(1,0) = H(0,1);

    return H;
}




Eigen::MatrixXd sample_mix4_gauss_direct_matrix(int N, uint32_t seed, bool truncate_to_unit_square) {
    constexpr double w1 = 0.25, w2 = 0.25, w3 = 0.25, w4 = 0.25;

    constexpr double mx1 = 0.3, my1 = 0.3;
    constexpr double mx2 = 0.3, my2 = 0.7;
    constexpr double mx3 = 0.7, my3 = 0.3;
    constexpr double mx4 = 0.7, my4 = 0.7;

    constexpr double s111 = 0.0064, s121 = 0.0,   s221 = 0.0064;
    constexpr double s112 = 0.0036, s122 = 0.0,   s222 = 0.0100;
    constexpr double s113 = 0.01,   s123 = 0.0,   s223 = 0.0036;
    constexpr double s114 = 0.01,   s124 = 0.009, s224 = 0.01;

    auto make_L = [](double s11, double s12, double s22) {
        Eigen::Matrix2d S;
        S << s11, s12,
             s12, s22;
        Eigen::LLT<Eigen::Matrix2d> llt(S);
        if (llt.info() != Eigen::Success) {
            S(0,0) += 1e-12; S(1,1) += 1e-12;
            llt.compute(S);
        }
        return llt.matrixL();
    };

    struct Comp2Dcol {
        double w;
        Eigen::Vector2d mu;   // colonna 2x1
        Eigen::Matrix2d L;
    };

    std::array<Comp2Dcol,4> comps = {{
        {w1, Eigen::Vector2d(mx1,my1), make_L(s111,s121,s221)},
        {w2, Eigen::Vector2d(mx2,my2), make_L(s112,s122,s222)},
        {w3, Eigen::Vector2d(mx3,my3), make_L(s113,s123,s223)},
        {w4, Eigen::Vector2d(mx4,my4), make_L(s114,s124,s224)}
    }};

    std::array<double,4> cdf = {
        comps[0].w,
        comps[0].w + comps[1].w,
        comps[0].w + comps[1].w + comps[2].w,
        comps[0].w + comps[1].w + comps[2].w + comps[3].w
    };

    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    std::normal_distribution<double> nd(0.0, 1.0);  // <-- QUI era il pezzo mancante

    Eigen::MatrixXd X(N, 2);
    int i = 0;

    while (i < N) {
        double u = u01(gen);
        int k = (u < cdf[0]) ? 0 : (u < cdf[1]) ? 1 : (u < cdf[2]) ? 2 : 3;

        Eigen::Vector2d z(nd(gen), nd(gen));
        Eigen::Vector2d x = comps[k].mu + comps[k].L * z;

        if (!truncate_to_unit_square ||
            (0.0 <= x(0) && x(0) <= 1.0 && 0.0 <= x(1) && x(1) <= 1.0)) {
            X(i,0) = x(0);
            X(i,1) = x(1);
            ++i;
        }
    }

    return X;
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
    }

    return map;
}







int main() {


    auto export_data_to_txt = [](const Eigen::Matrix<double, Eigen::Dynamic, 2>& data, const std::string& filename) {
        std::ofstream f(filename);
        if (!f) throw std::runtime_error("Cannot open " + filename);
        for (int i = 0; i < data.rows(); ++i)
            f << data(i,0) << " " << data(i,1) << "\n";
    };
    


    using coords_t = Eigen::Matrix<double, 1, 2>;
    using Mat2 = Eigen::Matrix2d;
    auto hessian_mix4_gauss = [&](coords_t coords) -> Eigen::Matrix<double, 2, 2> {
        const double x = coords(0);
        const double y = coords(1);

        // weights
        constexpr double w1 = 0.25, w2 = 0.25, w3 = 0.25, w4 = 0.25;

        // centers (means)
        constexpr double mx1 = 0.3, my1 = 0.3;
        constexpr double mx2 = 0.3, my2 = 0.7;
        constexpr double mx3 = 0.7, my3 = 0.3;
        constexpr double mx4 = 0.7, my4 = 0.7;

        // covariance entries
        constexpr double s111 = 0.0064, s121 = 0.0,   s221 = 0.0064;
        constexpr double s112 = 0.0036, s122 = 0.0,   s222 = 0.0100;
        constexpr double s113 = 0.01,   s123 = 0.0,   s223 = 0.0036;
        constexpr double s114 = 0.01,   s124 = 0.009, s224 = 0.01;

        Mat2 H = Mat2::Zero();
        H += w1 * hessian_gauss2(x, y, mx1, my1, s111, s121, s221);
        H += w2 * hessian_gauss2(x, y, mx2, my2, s112, s122, s222);
        H += w3 * hessian_gauss2(x, y, mx3, my3, s113, s123, s223);
        H += w4 * hessian_gauss2(x, y, mx4, my4, s114, s124, s224);
        return H;
    };
    auto make_Mraw_from_hessian_abs = [=](auto hessian_fun, double eps = 1e-3){ return [=](const coords_t& coords) -> Mat2 {
            Mat2 H  = hessian_fun(coords);
            Mat2 Hs = 0.5 * (H + H.transpose());

            Eigen::SelfAdjointEigenSolver<Mat2> es(Hs);
            if (es.info() != Eigen::Success) {
                return Mat2::Identity();
            }

            Eigen::Vector2d lam = es.eigenvalues().cwiseAbs();
            double w_min = 1e-4;   // non coarsen troppo nelle zone di alta curvatura
            double w_max = 1e+4;   // non raffinare all’infinito nelle zone piatte

            Eigen::Vector2d w = (lam.array() + eps).pow(-0.25).matrix();
            w(0) = std::clamp(w(0), w_min, w_max);
            w(1) = std::clamp(w(1), w_min, w_max);

            Mat2 R = es.eigenvectors();
            return R * w.asDiagonal() * R.transpose();
        };
    };




    //auto data = sample_mix4_gauss_direct_matrix(1000,42,true);
    auto data = read_locs_mtx("Meshes/Delaunay/locs_florida.mtx");  //locs_nord //2018.12.11_space_locs
    export_data_to_txt(data, "Meshes/Delaunay/data_points.txt");
    auto norditalia = read_locs_mtx("Meshes/Delaunay/bd_norditalia.mtx");
    auto florida = read_locs_mtx("Meshes/Delaunay/bd_florida.mtx");

    double min_angle = 20.0; 
    double max_area = 0.112; //  0.025/2;  // 0.25/9;  //0.02

    Delaunay<2, 2> del(std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2>>{florida}, min_angle, max_area, 0);   //caso 2

    //Delaunay<2, 2> del(std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2>>{rectangle_reduced}, 0);
    //del.refinement(20., del.domain_area()/6000);
    del.dcel().export_to_json("Meshes/Delaunay/delaunay_output.json");

    Delaunay<2,2>& delaunay = del;
    double prev_diameter = delaunay.min_triangle_diameter();
    double tol = 5e-2;

    for(int i=0; i< 0; ++i){      
        std::cout << "Iter " << i+1 << " of adaptivity." << std::endl;

        NodeMetric<2,2,AdaptiveStrategy::NodeDensity> node_metrics_obj(delaunay.dcel(), data, max_area);
        auto node_metrics = node_metrics_obj.node_density_knn_geodesic(); 
        node_metrics = normalize_matrix(node_metrics, delaunay.dcel(), 1.27, 1.);  // obiettivo: 750

        Adaptivity<2,2> adapt(delaunay, node_metrics);
        adapt.adaptivity_cycle(min_angle, max_area*100);
        adapt.dcel().export_to_json("Meshes/Delaunay/delaunay_output.json");  
        std::cout << adapt.dcel().n_cells() << " cells after adaptivity." << std::endl;

        // stopping criterion
        double curr_diameter = adapt.delaunay().min_triangle_diameter();
        std::cout << "curr: " << curr_diameter << "   prev: " << prev_diameter << std::endl;
        if(std::abs(curr_diameter - prev_diameter)/prev_diameter < tol){
            std::cout << "Converged: relative change in min triangle diameter is below threshold. Iter " << i+1 << std::endl;
            break;
        }
        
        prev_diameter = curr_diameter;
        delaunay = adapt.delaunay();
    }

    std::cout << "Number of nodes: " << delaunay.dcel().n_nodes() << std::endl;


    return 0;
}


