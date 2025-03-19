#ifndef __FDAPDE_DELAUNAY_H__
#define __FDAPDE_DELAUNAY_H__

#include "header_check.h"
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

    Delaunay(const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& boundary,
        const std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>>& holes = {})
    : boundary_points_(boundary),
    hole_points_(holes),
    dcel_(holes.empty() ? DCEL<local_dim, embed_dim>::make_polygon(boundary)
                    : DCEL<local_dim, embed_dim>::make_polygon(boundary, holes)) { }

    // Getter
    const DCEL<local_dim, embed_dim>& dcel() const { return dcel_; }
    DCEL<local_dim, embed_dim>& dcel() { return dcel_; }
    
    const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& boundary_points() const { return boundary_points_; }
    const std::vector<coords_t>& internal_points() const { return internal_points_; }
    // Setter
    void set_boundary_points(const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& points) {
        boundary_points_ = points;
    }

    const cell_t* find_triangle(const coords_t& P) {
        
        for (auto it = dcel_.cells_begin(); it != dcel_.cells_end(); ++it) {
            cell_t* cell = &(*it);  
    
            const coords_t& A = cell->halfedge()->node()->coords();
            const coords_t& B = cell->halfedge()->next()->node()->coords();
            const coords_t& C = cell->halfedge()->prev()->node()->coords();
    
            if (fdapde::internals::point_in_2d_tri(P, A, B, C)) {
                return cell;
            }
            // Checking if one point is on the edge
            if (fdapde::internals::contains(P, A, B) || 
            fdapde::internals::contains(P, B, C) || 
            fdapde::internals::contains(P, C, A)) {
            //i remove from internal_points vector
                auto it = std::remove_if(internal_points_.begin(), internal_points_.end(), 
                [&](const coords_t& point) { return point.isApprox(P); });
                internal_points_.erase(it, internal_points_.end());

                return nullptr;
            }
        }

    return nullptr;
    }

    void dig_cavity(const coords_t& u, halfedge_t* vw) { 
    //if we are on the boundary we add the triangle  
        if(vw->on_boundary()){
            dcel_.add_polygon(vw,u.transpose());
            return;
        }
    //finding the point adjacent to vw
        node_t* x = dcel_.adjacent(vw);
        if (!x) {
            return;
        }
        bool ccw = fdapde::internals::are_2d_counterclockwise_sorted(u, vw->node()->coords(), vw->twin()->node()->coords());
    //test of circumcircle   
        bool inside;
        if (ccw) {
            inside = fdapde::internals::in_circle(u, vw->node()->coords(), vw->twin()->node()->coords(), x->coords());
        } else {
            inside = fdapde::internals::in_circle(u, vw->twin()->node()->coords(), vw->node()->coords(), x->coords());
        }
        if (inside) { 
        //falied the test so remove the triangle and expand the cavity on the remaining edges
            halfedge_t* wv = vw->twin();
            halfedge_t* vx = vw->twin()->next();
            halfedge_t* xw = vw->twin()->prev(); 
            
            dcel_.remove_edge(vw);
            dig_cavity(u, vx);
            dig_cavity(u, xw);
        } else {
        //passed the test,adding the triangle 
            dcel_.add_polygon(vw,u.transpose());
            return;
        } 
    }
    
    
    void insert_vertex(const coords_t& u, const cell_t* triangle) {
        
        halfedge_t* vw = triangle->halfedge();
        halfedge_t* wx = vw->next();
        halfedge_t* xv = vw->prev();
        // expanding cavity
        dig_cavity(u, vw);
        dig_cavity(u, wx);
        dig_cavity(u, xv);
    }


    void build_triangulation(int num_points) {

        double min_x = boundary_points_.col(0).minCoeff();
        double max_x = boundary_points_.col(0).maxCoeff();
        double min_y = boundary_points_.col(1).minCoeff();
        double max_y = boundary_points_.col(1).maxCoeff();
    
        //creating the generator of casual points for the internal_points
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dist_x(min_x, max_x);
        std::uniform_real_distribution<double> dist_y(min_y, max_y);

        //generating internal_points 
        internal_points_.clear();
        internal_points_.reserve(num_points);

        while (static_cast<int>(internal_points_.size()) < num_points) {
            coords_t u;
            u << dist_x(gen), dist_y(gen);
            //verifies if the point is inside the domain (in order to control the concavities)
            if (fdapde::internals::point_in_polygon(boundary_points_, u)) {
                bool in_hole = false;
                for (const auto& hole : hole_points_) {
                    if (fdapde::internals::point_in_polygon(hole, u)) {
                        in_hole = true;
                        break;
                    }
                }
                if (!in_hole) {
                    internal_points_.push_back(u);
                }
            }
        }
    
        if (internal_points_.empty()) {
            return;
        }
    
        //inserting fist node in the domain and creating all the triangles from the boundary edges
        coords_t first_internal = internal_points_.front();
        std::cout<<first_internal<<std::endl;

        auto it = dcel_.halfedges_begin();
        for (int i = 0; i < boundary_points_.rows(); ++i, ++it) {
            halfedge_t* he = &(*it);
            dcel_.add_polygon(he, first_internal.transpose());
        }

        int node_offset = boundary_points_.rows(); 
        for (const auto& hole : hole_points_) {
            if (hole.rows() == 0) continue; 
            
            node_t* first_hole_node = std::addressof(*std::next(dcel_.nodes_begin(), node_offset)); 
            halfedge_t* first_hole_he = first_hole_node->halfedge();
            
            std::vector<halfedge_t*> hole_edges;
            halfedge_t* he = first_hole_he;
            do {
                hole_edges.push_back(he);
                he = he->next();
            } while (he != first_hole_he);  

            for (halfedge_t* he : hole_edges) {
                dcel_.add_polygon(he, first_internal);
            }

            node_offset += hole.rows(); 
        }
    //inserting the num_points-1 inner points in the domain
        for (size_t i = 1; i < internal_points_.size(); ++i) {
            coords_t u = internal_points_[i];
            const cell_t* triangle = find_triangle(u);

            if (!triangle) {
                continue;
            }
            insert_vertex(u, triangle);
        }
    //reordering id of cells and halfedges to cover some jumps between ids after removing
        int cont = 0;
        for (auto it = dcel_.cells_begin(); it != dcel_.cells_end(); ++it) {
            it->set_id(cont);
            cont++;
        }

        int cont_h = 0;
        for (auto it = dcel_.halfedges_begin(); it != dcel_.halfedges_end(); ++it) {
            it->set_id(cont_h);
            cont_h++;
        }
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
    std::vector<coords_t> internal_points_;
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>> hole_points_;
};
  
}  // namespace fdapde

#endif // __DELAUNAY_H__