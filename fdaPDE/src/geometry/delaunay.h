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
    //if the user wants to impose manually the internal points 
    void set_internal_points(const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& points) {
        internal_points_.clear();
        internal_points_.reserve(points.rows());
        for (int i = 0; i < points.rows(); ++i) {
            internal_points_.push_back(points.row(i));
        }
    }    

    halfedge_t* add_triangle(halfedge_t* v, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& node){
        return dcel_.add_polygon(v ,node);
    }

    halfedge_t* add_first_triangle(halfedge_t* v, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& node){                                                                       
        cell_t* c= v->cell();                                                                                    
        std::cout << "-------------------------------------------------------------------------------" <<std::endl;
        auto& cell_begin = *dcel_.cells_begin();
        std::vector<halfedge_t*> ghost_halfedges(3); 
    
        // add nodes and create ghost halfedges
        node_t* n;
        if(dcel_.find_node(node.row(0))==nullptr){  
                n = dcel_.insert_node(node_t(dcel_.n_nodes(), /* boundary = */ false, node.row(0)));
            }
        else {
                n = dcel_.find_node(node.row(0));  //if node is already a vertex of the mesh
        }
        if(dcel_.find_halfedge(n,c)){
            halfedge_t* h = dcel_.find_halfedge(n,c);
            ghost_halfedges[2] = h;
        }
        else if(dcel_.find_halfedge(n,&cell_begin)){
            halfedge_t* h = dcel_.find_halfedge(n,&cell_begin);
            ghost_halfedges[2] = h;
        }
        else{  
            ghost_halfedges[2] = dcel_.emplace_halfedge_(n);
            ghost_halfedges[2]->set_cell(c);
        }
        ghost_halfedges[0] = v;
        if(v->next()->cell()==ghost_halfedges[2]->cell())
            ghost_halfedges[1] = v->next();
        else{
            ghost_halfedges[1] = dcel_.find_halfedge(v->next()->node(),ghost_halfedges[2]->cell());
            if(!ghost_halfedges[1])
                return dcel_.insert_edge(dcel_.find_halfedge(v->prev()->prev()->node(), c),v);
                //ghost_halfedges[1]=dcel_.find_halfedge(v->prev()->node(), c);
        }
            

        for(auto v: ghost_halfedges){
           std::cout << "ghost_halfedges: " << v->id() << std::endl;
        }
        
        // add edges
        for (int i = 0; i < 3 ; ++i) {
            halfedge_t* h1 = ghost_halfedges[i];
            halfedge_t* h2 = ghost_halfedges[(i + 1) % (3)];
            
            // Check for intersection with boundary edges
            coords_t A = h1->node()->coords();
            coords_t B = h2->node()->coords();
            coords_t C;
            coords_t D;
            node_t* node_C;
            node_t* node_D;
            bool flag=false;
            if(!(h1->on_boundary() && h2->on_boundary())) {  //if the edges are not both on the boundary
                for (size_t i = 0; i < boundary_points_.rows(); ++i) {  
                    C = boundary_points_.row(i).transpose();
                    D = boundary_points_.row((i + 1) % boundary_points_.rows()).transpose();  
                    if (A != C && A != D && B != C && B != D) {  
                        if (fdapde::internals::intersect(A, B, C, D)) {  
                            std::cout << "Error: intersection with boundary!" << std::endl;
                            flag = true;
                            node_C = dcel_.find_node(C);
                            node_D = dcel_.find_node(D);
                            std::cout << "node_D: " << node_D->id() << std::endl;
                            break;
                        }
                    }
                }
            }
            if(!flag){
                halfedge_t* h = dcel_.insert_edge(h1, h2); 
                if(h && h!=h1)
                    ghost_halfedges[(i + 1) % (3)] = h->next();  //->next();
            } 
            else{
                //ghost_halfedges[(i + 1) % (3)] = dcel_.insert_edge(h1->prev(), h1->next())->twin();  //non vale sempre
                if(dcel_.find_halfedge(node_D, h1->next()->cell())){
                    std::cout << "qui" << std::endl;
                    ghost_halfedges[(i + 1) % (3)] = dcel_.insert_edge(h1->next(), dcel_.find_halfedge(node_D, h1->next()->cell()))->twin();
                }
                else{
                    std::cout << "C: " <<node_C->id() << std::endl;
                    ghost_halfedges[(i + 1) % (3)] = dcel_.insert_edge(h1->prev(), dcel_.find_halfedge(node_C, h1->next()->cell()))->twin();
                }
            }
        }
        
        //c->set_halfedge(v);

        return c->halfedge();
    }
    

    
    const cell_t* find_triangle(const coords_t& P) {
        
        for (auto it = dcel_.cells_begin(); it != dcel_.cells_end(); ++it) {
            cell_t* cell = &(*it);  
    
            const coords_t& A = cell->halfedge()->node()->coords();
            const coords_t& B = cell->halfedge()->next()->node()->coords();
            const coords_t& C = cell->halfedge()->prev()->node()->coords();

            // Checking if one point is on the edge
            if (fdapde::internals::contains(P, A, B) || 
            fdapde::internals::contains(P, B, C) || 
            fdapde::internals::contains(P, C, A)) {
            //I remove from internal_points vector
                auto it = std::remove_if(internal_points_.begin(), internal_points_.end(), 
                [&](const coords_t& point) { return point.isApprox(P); });
                internal_points_.erase(it, internal_points_.end());

                return nullptr;
            }
    
            if (fdapde::internals::point_in_2d_tri(P, A, B, C)) {
                return cell;
            }
        }

    return nullptr;
    }

    void dig_cavity(const coords_t& u, halfedge_t* vw) { 
    //if we are on the boundary we add the triangle  
        if(vw->on_boundary()){
            add_triangle(vw,u.transpose());
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
            add_triangle(vw,u.transpose());
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
            std::cout << "nuovo punto interno: " << u.transpose() << std::endl;
            //verifies if the point is inside the domain (in order to control the concavities)
            //and discard the points falling on the boundary of the domain
            if (fdapde::internals::point_in_polygon(boundary_points_, u)){
            /* {
                bool in_hole = false;
                for (const auto& hole : hole_points_) {
                    if (fdapde::internals::point_in_polygon(hole, u)) {
                        in_hole = true;
                        break;
                    }
                }
                if (!in_hole) {*/
                std::cout << "nuovo punto inserito: " << u.transpose() << std::endl;
                    internal_points_.push_back(u);}
              //  }
          //  }
        }
    
        if (internal_points_.empty()) {
            return;
        }
    
        //inserting fist node in the domain and creating all the triangles from the boundary edges
        coords_t first_internal = internal_points_.front();
        
        auto it = dcel_.halfedges_begin();
        for (int i = 0; i < boundary_points_.rows(); ++i, ++it) {
            halfedge_t* he = &(*it);
            add_first_triangle(he, first_internal.transpose());
        }
/*
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
                add_first_triangle(he, first_internal);
            }

            node_offset += hole.rows(); 
        }
*/
        flip();
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
        dcel_.set_n_cells_(cont);

        int cont_h = 0;
        for (auto it = dcel_.halfedges_begin(); it != dcel_.halfedges_end(); ++it) {
            it->set_id(cont_h);
            cont_h++;
        }

        Triangulation<local_dim, embed_dim> triangulation = DCEL_to_Triangulation();

        std::string filename = "mesh_output.txt";
        export_triangulation_to_txt(triangulation, filename);

        std::string command = "python3 fdaPDE/src/plot_mesh.py";
        std::system(command.c_str()); 
    }

//overloaded one if user wants to impose internal points manually 
//if one point exceeds the domain find_triangle return nullptr and does not enter in the triangulation
    void build_triangulation() {

        if (internal_points_.empty()) {
            return;
        }
        //inserting fist node in the domain and creating all the triangles from the boundary edges
        coords_t first_internal = internal_points_.front();
        auto it = dcel_.halfedges_begin();
        for (int i = 0; i < boundary_points_.rows(); ++i, ++it) {
            halfedge_t* he = &(*it);
            add_first_triangle(he, first_internal.transpose());
        }
/*
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
                add_first_triangle(he, first_internal);
            }

            node_offset += hole.rows(); 
        }*/

    // flip the initial trinagulation if not Delaunay 
       flip();
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

        dcel_.set_n_cells_(cont);

        int cont_h = 0;
        for (auto it = dcel_.halfedges_begin(); it != dcel_.halfedges_end(); ++it) {
            it->set_id(cont_h);
            cont_h++;
        }

        Triangulation<local_dim, embed_dim> triangulation = DCEL_to_Triangulation();

        std::string filename = "mesh_output.txt";
        export_triangulation_to_txt(triangulation, filename);

        std::string command = "python3 fdaPDE/src/plot_mesh.py";
        std::system(command.c_str()); 
    }


    void flip() {

        // Creating a list of halfedges to check wheter they are locally delaunay or not (in this case flippable)
        std::list<halfedge_t*> halfedges_to_check;
        for (auto it = dcel_.halfedges_begin(); it != dcel_.halfedges_end(); ++it,++it) {
            if(!it->on_boundary())
            halfedges_to_check.push_back(&(*it));
        }
        // flip algorithm
        while (!halfedges_to_check.empty()) {
            halfedge_t* edge = halfedges_to_check.front();
            halfedges_to_check.pop_front();

            cell_t* neighbor = edge->twin()->cell();
            std::cout<<"cella vicina : "<<neighbor->id()<<" ad halfedge: "<<edge->id()<<std::endl;
    
            // obtaining the 4 vertices of the quadrilateral formed by the two adjoining triangles 
            coords_t A = edge->node()->coords();
            coords_t B = edge->twin()->node()->coords();
            coords_t C = edge->prev()->node()->coords();
            coords_t D = edge->twin()->prev()->node()->coords();
    
            if (fdapde::internals::in_circle(A, B, C, D) || fdapde::internals::in_circle(A, D, B, C)) {
                std::cout<<"SONO QUA"<<std::endl;
                std::cout <<"A: "<<edge->node()->id()<<" B: "<<edge->twin()->node()->id()<<" C: "<<edge->prev()->node()->id()<<" D: "<<edge->twin()->prev()->node()->id()<<std::endl;
            
                // we flip since edge is not locally delaunay
                halfedge_t* e = edge;
                dcel_.remove_edge(edge);
                halfedge_t* new_edge = dcel_.insert_edge(e->prev(), e->twin()->prev());
            
                if (new_edge) {
                    // Inserting the new halfedges created by the flip
                    if(!new_edge->prev()->on_boundary())
                    halfedges_to_check.push_back(new_edge->prev());
                    if(!new_edge->next()->on_boundary())
                    halfedges_to_check.push_back(new_edge->next());
                    if(!new_edge->twin()->prev()->on_boundary())
                    halfedges_to_check.push_back(new_edge->twin()->prev());
                    if(!new_edge->twin()->next()->on_boundary())
                    halfedges_to_check.push_back(new_edge->twin()->next());
                }
            
            }
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

    //function to convert the dcel into a triangulation
    Triangulation<local_dim, embed_dim> DCEL_to_Triangulation() {  
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> nodes(dcel_.n_nodes(), embed_dim);
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> cells(dcel_.n_cells(), 3);
    Eigen::Matrix<int, Eigen::Dynamic, 1> boundary_markers(dcel_.n_nodes());

    // Fill nodes matrix
    int node_idx = 0;
    for (auto it = dcel_.nodes_begin(); it != dcel_.nodes_end(); ++it) {
        nodes.row(node_idx) = it->coords().transpose();
        boundary_markers(node_idx) = it->on_boundary() ? 1 : 0;
        node_idx++;
    }

    // Fill cells matrix
    int cell_idx = 0;
    for (auto it = dcel_.cells_begin(); it != dcel_.cells_end(); ++it) {
        halfedge_t* h = it->halfedge();
        for (int i = 0; i < 3; ++i) {
            cells(cell_idx, i) = h->node()->id();
            h = h->next();
        }
        cell_idx++;
    }

    // Create Triangulation object
    Triangulation<local_dim, embed_dim> triangulation(nodes, cells, boundary_markers);
    return triangulation;
    }

    void export_triangulation_to_txt(const Triangulation<LocalDim, EmbedDim>& triangulation, const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            return;
        }

        file << "Nodes:\n";
        for (int i = 0; i < triangulation.n_nodes(); ++i) {
            auto coords = triangulation.node(i);
            int marker = triangulation.is_node_on_boundary(i) ? 1 : 0;
            file << i << " " << coords(0) << " " << coords(1) << " " << marker << "\n";
        }

        file << "\nCells:\n";
        for (int i = 0; i < triangulation.n_cells(); ++i) {
            auto cell = triangulation.cells().row(i);;
            file << i << " " << cell(0) << " " << cell(1) << " " << cell(2) << "\n";
        }

        file.close();
    }

   private:
    DCEL<local_dim, embed_dim> dcel_;  
    Eigen::Matrix<double, Eigen::Dynamic, embed_dim> boundary_points_;
    std::vector<coords_t> internal_points_;
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>> hole_points_;
};
  
}  // namespace fdapde

#endif // __DELAUNAY_H__