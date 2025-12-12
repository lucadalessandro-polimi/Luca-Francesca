#include <fstream>
#include <nlohmann/json.hpp> 
using json = nlohmann::json;
#include <geometry.h>
using namespace fdapde;
#include "domains.h"

using namespace std::chrono;

using Mat2 = Eigen::Matrix<double, Eigen::Dynamic, 2>;

// rettangolo [0,1]x[0,1]
/*static inline Mat2 make_rectangle_unit()
{
    Mat2 rect(4,2);
    rect << 0.0, 0.0,
            1.0, 0.0,
            1.0, 1.0,
            0.0, 1.0;
    return rect;
}

// genera n punti uniformi nel rettangolo axis-aligned [xmin,xmax] x [ymin,ymax]
static inline Mat2 random_points_in_rectangle(std::size_t n, double xmin, double xmax, double ymin, double ymax, uint64_t seed = 42)
{
    Mat2 P(n,2);
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> Ux(xmin, xmax);
    std::uniform_real_distribution<double> Uy(ymin, ymax);
    for (std::size_t i=0; i<n; ++i) {
        P(i,0) = Ux(rng);
        P(i,1) = Uy(rng);
    }
    return P;
}

// costruisce mesh con diversa densità variando la soglia area per la refinement
// se area_thresh <= 0 -> NESSUN refinement (rimangono 4 nodi)
template<int LocalDim, int EmbedDim>
static inline Delaunay<LocalDim, EmbedDim>
build_mesh(const std::vector<Eigen::Matrix<double, Eigen::Dynamic, EmbedDim>>& boundary,double area_thresh_scale) // fraction
{
    Delaunay<LocalDim, EmbedDim> del(boundary, 0);
    if (area_thresh_scale > 0.0) {
        const double area = del.domain_area();
        del.refinement(10, area * area_thresh_scale);
    }
    return del;
}*/

int main() {

    /*using LD = std::integral_constant<int,2>;
    using ED = std::integral_constant<int,2>;
    // ----------------- setup rettangolo -----------------
    Eigen::Matrix<double, Eigen::Dynamic, 2> rectangle = make_rectangle_unit();
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2>> bd = { rectangle };
    const double xmin = 0.0, xmax = 1.0, ymin = 0.0, ymax = 1.0;
    const std::vector<std::size_t> data_sizes = { 500, 1000, 2000, 4000};
    // ----------------- livelli di mesh -----------------
    const std::vector<double> mesh_area_scales = {
        1.0/10.0,
        1.0/100.0,
        1.0/1000.0
    };*/

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

    // Punti attorno a una RETTA y = m x + q (x ∈ [xmin, xmax])
    const auto random_points_around_line = [](int n,double m, double q,double xmin, double xmax,double sigma_perp,double sigma_tang = 0.0,unsigned long long seed = 42ULL)-> Eigen::Matrix<double, Eigen::Dynamic, 2>
    {
        Eigen::Matrix<double, Eigen::Dynamic, 2> P(n, 2);
        std::mt19937_64 rng(seed);
        std::uniform_real_distribution<double> Ux(xmin, xmax);
        std::normal_distribution<double> Z(0.0, 1.0);

        double denom = std::sqrt(1.0 + m*m);
        Eigen::Vector2d t_hat(1.0/denom, m/denom);
        Eigen::Vector2d n_hat(-m/denom, 1.0/denom);

        for (int i = 0; i < n; ++i) {
            double x = Ux(rng);
            double y = m * x + q;
            Eigen::Vector2d base(x, y);
            double z_perp = Z(rng);
            double z_tang = (sigma_tang > 0.0 ? Z(rng) : 0.0);
            Eigen::Vector2d p = base + sigma_perp * z_perp * n_hat + sigma_tang * z_tang * t_hat;
            P.row(i) = p;
        }
        return P;
    };




    // ------------------ PARAMETRI GLOBALI -------------------
    int n_points = 200;           // numero di punti dati
    double xmin = 1500, xmax = 2500.0;
    double ymin = 500.0, ymax = 1500.0;
    unsigned long long seed = 12345; //42ULL;

    // Dominio rettangolare
    Eigen::MatrixXd rectangle(4, 2);
    rectangle << xmin, ymin,
                 xmax, ymin,
                 xmax, ymax,
                 xmin, ymax;

    // ------------------ CREAZIONE DATI -------------------

    // ------------------ FILE DI OUTPUT -------------------
    std::ofstream csv("Meshes/Delaunay/complexity_test.csv");
    csv << "Lmin,Lmax,tol, n_data, n_iter,n_nodes,time_ms\n";
    //csv << "tol, n_data, n_iter,n_nodes,time_ms\n";

    // ------------------ PARAMETRI DI TEST -------------------
    std::vector<double> Lmins = {0.7};  //{0.6, 0.7, 0.8};
    std::vector<double> Lmaxs = {1.3};  //{1.2, 1.3, 1.4};
    std::vector<double> tols  = {0.001}; //{0.05, 0.1, 0.5};
    std::vector<int> n_data = {10, 100, 1000, 10000};

    // ------------------ LOOP PRINCIPALE -------------------
    for (double Lmin : Lmins) {
        for (double Lmax : Lmaxs) {
            for (double tol : tols) {
                for (int n : n_data) {

                    std::cout << "\n=== Running test with Lmin=" << Lmin << ", Lmax=" << Lmax << ", tol=" << tol << ", n_data=" << n << " ===" << std::endl;
                    //std::cout << "\n=== Running test with " << "tol=" << tol << ", n_data=" << n << " ===" << std::endl;

                    // 1️⃣ Crea mesh iniziale
                    Delaunay<2, 2> del({rectangle}, 0);
                    del.refinement(20, del.domain_area() / 100.0);
                    Eigen::MatrixXd data = random_points_in_rectangle(n, xmin, xmax, ymin, ymax, seed);

                    // 2️⃣ Inizializza Adaptivity
                    auto start = high_resolution_clock::now();
                    Adaptivity<2, 2> adapt(del, data, 20, del.domain_area() / 100.0, tol);

                    // 3️⃣ Imposta parametri metrici
                    adapt.metric().set_Lmin(Lmin);
                    adapt.metric().set_Lmax(Lmax);

                    // 4️⃣ Esegui ciclo adattivo
                    int n_iter = adapt.adaptivity_cycle_two(20, del.domain_area() / 100.0, tol);   
                    
                    // 5️⃣ Misura tempo
                    auto end = high_resolution_clock::now();
                    double ms = duration_cast<milliseconds>(end - start).count();

                    // 6️⃣ Conta vertici finali
                    int n_nodes = adapt.dcel().n_nodes();

                    std::cout << "Completed in " << n_iter << " iterations, "<< n_nodes << " nodes, time: " << ms << " ms.\n";

                    // 7️⃣ Salva su file
                    csv << Lmin << "," << Lmax << "," << tol << "," << n << "," << n_iter << "," << n_nodes << "," << ms << "\n";
                    //csv << tol << "," << n << "," << n_iter << "," << n_nodes << "," << ms << "\n";

                }
            }
        }
    }

    csv.close();
    std::cout << "\nComplexity test finished. Results saved to complexity_test.csv\n";

    /*std::ofstream out("Meshes/Delaunay/timing_adaptivity.csv");
    out << "DataN,AreaScale,NumPoints,Time_ms\n";
    for (double scale : mesh_area_scales) {
        for (std::size_t N : data_sizes) {
            auto del = build_mesh<LD::value, ED::value>(bd, scale);
            const std::size_t mesh_nodes = del.dcel().n_nodes();

            // genera dati Nx2 nel rettangolo
            Mat2 data = random_points_in_rectangle(N, xmin, xmax, ymin, ymax, 123456 + N);

            // ------------ misura tempo di costruzione Adaptivity ------------
            auto t0 = high_resolution_clock::now();

            auto n_nodes = del.dcel().n_nodes();
            Adaptivity<LD::value, ED::value> adapt(del, data, n_nodes -1);  //4

            // se vuoi misurare anche una fase di esecuzione (es. adapt.run()):
            // adapt.run();  // <--- decommenta se esiste

            auto t1 = high_resolution_clock::now();
            using ms_f = std::chrono::duration<double, std::milli>;
            double ms = std::chrono::duration_cast<ms_f>(t1 - t0).count();

            std::cout << "DataN=" << N << ", NumPoints=" << mesh_nodes - adapt.dcel().n_nodes()  << " -> Adaptivity ctor time = " << ms << " ms\n";
            //out << N << "," << scale << "," << mesh_nodes - adapt.dcel().n_nodes()  << "," << ms << "\n";
            out << N << "," << scale << "," << adapt.dcel().n_nodes()  << "," << ms << "\n";
        }
    }
    out.close();
    std::cout << "Wrote CSV to Meshes/Delaunay/timing_adaptivity.csv\n";*/
    return 0;
}
