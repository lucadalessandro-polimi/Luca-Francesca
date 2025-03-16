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
        double tol = 1e-6; // Tolleranza numerica per gestire il caso in cui il punto è su un lato
        
        for (auto it = dcel_.cells_begin(); it != dcel_.cells_end(); ++it) {
            cell_t* cell = &(*it);  // Ora il tipo corrisponde correttamente
    
            const coords_t& A = cell->halfedge()->node()->coords();
            const coords_t& B = cell->halfedge()->next()->node()->coords();
            const coords_t& C = cell->halfedge()->prev()->node()->coords();
    
            if (is_point_inside_triangle(P, A, B, C)) {
                return cell;
            }
        // Se il punto è su un lato, rimuoviamo il lato corrispondente
        halfedge_t* removed_edge = nullptr;
        if (is_point_on_edge(P, A, B, tol)) {
            removed_edge = cell->halfedge();
        } else if (is_point_on_edge(P, B, C, tol)) {
            removed_edge = cell->halfedge()->next();
        } else if (is_point_on_edge(P, C, A, tol)) {
            removed_edge = cell->halfedge()->prev();
        }

        if (removed_edge) {
            std::cout << "⚠️ Il punto è su un lato, rimuovo l'edge " << removed_edge->id() << std::endl;
            dcel_.remove_edge(removed_edge);

            // **Scorriamo sugli half-edges della cella e aggiungiamo P con add_polygon**
            halfedge_t* h = cell->halfedge();
            halfedge_t* h1 = cell->halfedge()->next();
            halfedge_t* h2 = cell->halfedge()->next()->next();
            halfedge_t* h3 = cell->halfedge()->next()->next()->next();

            dcel_.add_polygon(h, P.transpose());
            dcel_.add_polygon(h1, P.transpose());
            dcel_.add_polygon(h2, P.transpose());
            dcel_.add_polygon(h3, P.transpose());

            //dig_cavity(P.transpose(),h);
            //dig_cavity(P.transpose(),h1);
            //dig_cavity(P.transpose(),h2);
            //dig_cavity(P.transpose(),h3);
        /*    do {
                std::cout<<"CIAOOOOOOOOOO:"<<h->id()<<std::endl;
                dcel_.add_polygon(h, P.transpose());  // Collegamento con il nuovo punto
                h = h->next();
            } while (h != cell->halfedge());*/

            return nullptr; // Restituiamo la cella dopo l’aggiornamento
        }

       }
        return nullptr;
    }

    bool is_point_on_edge(const coords_t& P, const coords_t& A, const coords_t& B, double tol) {
        double cross = (P.y() - A.y()) * (B.x() - A.x()) - (P.x() - A.x()) * (B.y() - A.y());
        if (std::abs(cross) > tol) return false; // Non è collineare
    
        double dot = (P.x() - A.x()) * (B.x() - A.x()) + (P.y() - A.y()) * (B.y() - A.y());
        if (dot < 0) return false; // Punto fuori dal segmento
    
        double len_sq = (B.x() - A.x()) * (B.x() - A.x()) + (B.y() - A.y()) * (B.y() - A.y());
        if (dot > len_sq) return false; // Punto fuori dal segmento
    
        return true; // Il punto è sul segmento
    }
    


    void dig_cavity(const coords_t& u, halfedge_t* vw) { 
        if(vw->on_boundary()){
            std::cout<<"SONO AL BRODO CON : "<<vw->id()<<std::endl;
            dcel_.add_polygon(vw,u.transpose());
            return;
        }

        node_t* x = dcel_.adjacent(vw);
        if (!x) {
            std::cout << "Half-edge " << vw->id() << " non ha un nodo adiacente.\n";
            return;
        }
        std::cout << "Nodo adiacente a half-edge " << vw->id() << ": " << x->id() << "\n";
        
        std::cout << "Controllo in_circle su nodo " << x->id() << " rispetto al punto inserito " << u.transpose() << std::endl;
        std::cout<<"CAZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZO:   "<<vw->id()<<std::endl;
        std::cout<<"A: "<<u<<" B: "<<vw->node()->id()<<" C: "<<vw->twin()->node()->id()<<" X: "<<x->id()<<std::endl;
        if (in_circle(u, vw->node()->coords(), vw->twin()->node()->coords(), x->coords())) {  
       // if (in_circle(u, vw->twin()->node()->coords(), vw->node()->coords(), x->coords())) {  
            std::cout << "Punto " << x->id() << " è dentro il circumcerchio";
            
            halfedge_t* wv = vw->twin();
            halfedge_t* vx = vw->twin()->next();
            halfedge_t* xw = vw->twin()->prev(); 

            std::cout << "Half-edges del triangolo: " << std::endl;
            std::cout << " - wv: " << wv->id() << " | vx: " << vx->id() << " | xw: " << xw->id() << std::endl;
    
            dcel_.remove_edge(vw);
            dig_cavity(u, vx);
            dig_cavity(u, xw);
           
            std::cout<<"REMOVE: "<<vw->twin()->id()<<std::endl;
            //dcel_.remove_polygon(vw->twin()->cell());
            //dcel_.remove_edge(vw);
        } else {
            std::cout << "Punto " << x->id() << " NON è dentro il circumcerchio. Aggiungo nuovo triangolo.\n";
            dcel_.add_polygon(vw,u.transpose());
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
        //dcel_.remove_polygon(triangle);
    
    }

    void build_triangulation() {
        std::cout << "🔷 Inizio costruzione della triangolazione...\n";
        
        for (int i = 1; i < internal_points_.rows(); ++i) {
            coords_t u = internal_points_.row(i);
            std::cout << "🔍 Inserimento del punto interno: " << u.transpose() << std::endl;
            
            const cell_t* triangle = find_triangle(u);
            
            if (!triangle) {
                std::cerr << "❌ Errore: Nessun triangolo trovato per il punto " << u.transpose() << "!" << std::endl;
                continue;
            }
            
            insert_vertex(u, triangle);
        }
        int cont = 0;
        for(auto it = dcel_.cells_begin();it!=dcel_.cells_end();++it){
            it->set_id(cont);
            cont++;
        }
/*
        int cont_h = 0;
        for(auto it = dcel_.halfedges_begin();it!=dcel_.halfedges_end();++it){
            it->set_id(cont_h);
            cont_h++;
        }
 */      
        std::cout << "✅ Triangolazione completata con successo!\n";
    }

    void initialize_triangulation() {
        std::cout << "🔷 Inizializzazione della triangolazione...\n";
    
        if (internal_points_.rows() == 0) {
            std::cerr << "❌ Errore: Nessun punto interno disponibile per inizializzare la triangolazione!\n";
            return;
        }
        Eigen::Matrix<double, 1, embed_dim> first_internal = internal_points_.row(0);
        auto it = dcel_.halfedges_begin();
        int count = dcel_.n_nodes();

        for (int i = 0; i < boundary_points_.rows(); ++i, ++it) {
            halfedge_t* he = &(*it);
            dcel_.add_polygon(he, first_internal);
        }
    
        std::cout << "✅ Triangolazione iniziale completata!\n";
    }
    
    void print_dcel() {
        std::cout << "==============================" << std::endl;
        std::cout << "📌 STATO ATTUALE DELLA DCEL 📌" << std::endl;
        std::cout << "==============================" << std::endl;
    
        // 📍 Stampa tutti i nodi
        std::cout << "\n🟢 NODI: \n";
        for (auto it = dcel_.nodes_begin(); it != dcel_.nodes_end(); ++it) {
            std::cout << "ID: " << it->id() << " | Coords: (" << it->coords()(0) << ", " << it->coords()(1) << ")"
                      << (it->on_boundary() ? " [BOUNDARY]" : "") << std::endl;
        }
    
        // 🔗 Stampa tutti gli Half-Edges
        std::cout << "\n🔵 HALF-EDGES: \n";
        for (auto it = dcel_.halfedges_begin(); it != dcel_.halfedges_end(); ++it) {
            std::cout << "ID: " << it->id()
                      << " | Nodo Origine: " << (it->node() ? std::to_string(it->node()->id()) : "NULL")
                      << " | Twin: " << (it->twin() ? std::to_string(it->twin()->id()) : "NULL")
                      << " | Next: " << (it->next() ? std::to_string(it->next()->id()) : "NULL")
                      << " | Prev: " << (it->prev() ? std::to_string(it->prev()->id()) : "NULL")
                      << std::endl;
        }
    
        // 🔳 Stampa tutte le Celle
        std::cout << "\n🟠 CELLE: \n";
        for (auto it = dcel_.cells_begin(); it != dcel_.cells_end(); ++it) {
            std::cout << "Cella ID: " << it->id() << " | Half-edge di riferimento: "
                      << (it->halfedge() ? std::to_string(it->halfedge()->id()) : "NULL") << std::endl;
            if (it->halfedge()) {
                halfedge_t* h = it->halfedge();
                std::cout << "  🔗 Half-edges nella cella: ";
                halfedge_t* start = h;
                do {
                    std::cout << h->id() << " ";
                    h = h->next();
                } while (h && h != start);
                std::cout << std::endl;
            }
        }
    
        std::cout << "==============================\n" << std::endl;
    }
    
   private:
    DCEL<local_dim, embed_dim> dcel_;  
    Eigen::Matrix<double, Eigen::Dynamic, embed_dim> boundary_points_;
    Eigen::Matrix<double, Eigen::Dynamic, embed_dim> internal_points_;
};
  
}  // namespace fdapde

#endif // __DELAUNAY_H__