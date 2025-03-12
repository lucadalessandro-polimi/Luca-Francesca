#ifndef __FDAPDE_DELAUNAY_H__
#define __FDAPDE_DELAUNAY_H__

#include "src/geometry/dcel.h"
//#include "src/geometry/primitives.h" //finche non risolviamo questione include non usiamo le due fun in primitives.h

namespace fdapde {
  
template <int LocalDim = 2, int EmbedDim = 2>
class Delaunay {
   public:
    static constexpr int local_dim = LocalDim;
    static constexpr int embed_dim = EmbedDim;

    using coords_t = Eigen::Matrix<double, embed_dim, 1>;
    using node_t = typename DCEL<local_dim, embed_dim>::node_t;
    using halfedge_t = typename DCEL<local_dim, embed_dim>::halfedge_t;
    using cell_t = typename DCEL<local_dim, embed_dim>::cell_t;

    // Costruttori
    // Costruttore di default
    Delaunay() = default;
    
    // Costruttore che inizializza la DCEL con il bordo e memorizza i punti interni
    Delaunay(const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& boundary,
        const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& internal)
    : boundary_points_(boundary),
    internal_points_(internal),
    dcel_(DCEL<local_dim, embed_dim>::make_polygon(boundary)) { }



    // Getter
    const DCEL<local_dim, embed_dim>& dcel() const { return dcel_; }
    // OVERLOAD NON CONSTANTE PER LAVORARE NEL MAIN 
    DCEL<local_dim, embed_dim>& dcel() { return dcel_; }
    
    const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& boundary_points() const { return boundary_points_; }
    const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& internal_points() const { return internal_points_; }

    void set_boundary_points(const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& points) {
        boundary_points_ = points;
    }

    void set_internal_points(const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& points) {
        internal_points_ = points;
    }


    bool in_circle(const coords_t& A, const coords_t& B, const coords_t& C, const coords_t& D) const {
  
    
        // Costruzione della matrice 3x3 
        Eigen::Matrix3d M;
        M << (A.x() - D.x()), (A.y() - D.y()), (A.x() - D.x()) * (A.x() - D.x()) + (A.y() - D.y()) * (A.y() - D.y()),
             (B.x() - D.x()), (B.y() - D.y()), (B.x() - D.x()) * (B.x() - D.x()) + (B.y() - D.y()) * (B.y() - D.y()),
             (C.x() - D.x()), (C.y() - D.y()), (C.x() - D.x()) * (C.x() - D.x()) + (C.y() - D.y()) * (C.y() - D.y());
    
        double det = M.determinant();
        
        std::cout << "Determinante InCircle: " << det << " -> " << (det > 0 ? "Dentro" : "Fuori") << std::endl;
        
        return det > 0;
    }
    

    bool is_point_inside_triangle(const coords_t& P, const coords_t& A, const coords_t& B, const coords_t& C) const{
        // baricenter test

        Eigen::Matrix<double, LocalDim, LocalDim> M;
        M << (B - A), (C - A); 
        coords_t lambda = M.inverse() * (P - A);
        double lambda1 = lambda.x();
        double lambda2 = lambda.y();
        double lambda3 = 1 - lambda1 - lambda2;
        return (lambda1 >= 0 && lambda1<=1 && lambda2 >= 0 && lambda2<=1 && lambda3 >= 0 && lambda3<=1);  
    }

    const cell_t* find_triangle(const coords_t& P) {
        for (auto it = dcel_.cells_begin(); it != dcel_.cells_end(); ++it) {
            cell_t* cell = &(*it);  // Ora il tipo corrisponde correttamente
    
            const coords_t& A = cell->halfedge()->node()->coords();
            const coords_t& B = cell->halfedge()->next()->node()->coords();
            const coords_t& C = cell->halfedge()->prev()->node()->coords();
    
            if (is_point_inside_triangle(P, A, B, C)) {
                return cell;
            }
        }
        return nullptr;
    }

    void remove_triangle(const cell_t* triangle) {
        dcel_.remove_polygon(triangle);  // Riutilizziamo remove_polygon
    }


    void dig_cavity(const coords_t& u, halfedge_t* vw) { 
        node_t* x = dcel_.adjacent(vw);
        if (!x) {
            std::cout << "Half-edge " << vw->id() << " non ha un nodo adiacente.\n";
            return;
        }
        std::cout << "Nodo adiacente a half-edge " << vw->id() << ": " << x->id() << "\n";
        
        std::cout << "Controllo in_circle su nodo " << x->id() << " rispetto al punto inserito " << u.transpose() << std::endl;
        
        if (in_circle(u, vw->node()->coords(), vw->next()->node()->coords(), x->coords())) {  
            std::cout << "Punto " << x->id() << " è dentro il circumcerchio. Rimuovo triangolo e continuo.\n";
            
            halfedge_t* wv = vw->twin();
            halfedge_t* vx = vw->twin()->next();
            halfedge_t* xw = vw->twin()->prev(); 

            std::cout << "Half-edges del triangolo: " << std::endl;
            std::cout << " - wv: " << wv->id() << " | vx: " << vx->id() << " | xw: " << xw->id() << std::endl;
    
            dig_cavity(u, vx);
            dig_cavity(u, xw);

            remove_triangle(vw->twin()->cell());
        } else {
            std::cout << "Punto " << x->id() << " NON è dentro il circumcerchio. Aggiungo nuovo triangolo.\n";
            dcel_.add_triangle(vw->node()->coords(),vw->next()->node()->coords(),u);
            return;
        } 
    }
    
    
    void insert_vertex(const coords_t& u, const cell_t* triangle) {
        std::cout << "Inserimento vertice: " << u.transpose() << " nel triangolo con ID: " << triangle->id() << std::endl;
        
        halfedge_t* vw = triangle->halfedge();
        halfedge_t* wx = vw->next();
        halfedge_t* xv = vw->prev();
    
        std::cout << "Half-edges del triangolo: " << std::endl;
        std::cout << " - vw: " << vw->id() << " | wx: " << wx->id() << " | xv: " << xv->id() << std::endl;
        
        // Espandiamo la cavità e ricostruiamo la triangolazione
        dig_cavity(u, vw);
        dig_cavity(u, wx);
        dig_cavity(u, xv);

        // Rimuoviamo il triangolo attuale
        remove_triangle(triangle);
        std::cout << "Triangolo rimosso correttamente.\n";
    }

    void build_triangulation() {
        std::cout << "🔷 Inizio costruzione della triangolazione...\n";
        
        for (int i = 0; i < internal_points_.rows(); ++i) {
            coords_t u = internal_points_.row(i);
            std::cout << "🔍 Inserimento del punto interno: " << u.transpose() << std::endl;
            
            const cell_t* triangle = find_triangle(u);
            
            if (!triangle) {
                std::cerr << "❌ Errore: Nessun triangolo trovato per il punto " << u.transpose() << "!" << std::endl;
                continue;
            }
            
            insert_vertex(u, triangle);
        }
        
        std::cout << "✅ Triangolazione completata con successo!\n";
    }
    
    void initialize_triangulation() {
        std::cout << "🔷 Inizializzazione triangolazione...\n";
    
        int n = boundary_points_.rows();
        if (n < 3) {
            std::cerr << "❌ Errore: Non ci sono abbastanza punti per formare un triangolo!\n";
            return;
        }
    
        // Creiamo una prima triangolazione con i punti al bordo
        for (int i = 1; i < n - 1; ++i) {
            dcel_.add_triangle(boundary_points_.row(0), boundary_points_.row(i), boundary_points_.row(i + 1));
        }
    
        std::cout << "✅ Triangolazione iniziale completata!\n";
    }

   private:
    DCEL<local_dim, embed_dim> dcel_;  
    Eigen::Matrix<double, Eigen::Dynamic, embed_dim> boundary_points_;
    Eigen::Matrix<double, Eigen::Dynamic, embed_dim> internal_points_;
};
  
}   // namespace fdapde

#endif // __DELAUNAY_H__