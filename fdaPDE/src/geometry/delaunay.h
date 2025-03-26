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

    void remove_triangle(const cell_t* cell){
        dcel_.remove_polygon(cell);
        return;
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
            // obtaining the 4 vertices of the quadrilateral formed by the two adjoining triangles 
            coords_t A = edge->node()->coords();
            coords_t B = edge->twin()->node()->coords();
            coords_t C = edge->prev()->node()->coords();
            coords_t D = edge->twin()->prev()->node()->coords();
    
            if (fdapde::internals::in_circle(A, B, C, D) || fdapde::internals::in_circle(A, D, B, C)) {
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


///////////////////////////////////// NEW STORY ////////////////////////////////////////////////7

    // Function to mark the cavity during insertion
    void mark_cavity(node_t* u, halfedge_t* vw, std::vector<halfedge_t*>& D, std::vector<halfedge_t*>& C) {
        if(vw->on_boundary()){
            C.push_back(vw);  
            std::cout<<"SONO QUA"<<std::endl;
            return;
        }

        node_t* x = dcel_.adjacent(vw);
        if (!x) {
            return;
        }
        std::cout<<"nodo adiacente: "<<x->id()<<std::endl; 
        if(vw->twin()->cell()->visited()) return; 
        bool ccw = fdapde::internals::are_2d_counterclockwise_sorted(u->coords(), vw->node()->coords(), vw->twin()->node()->coords());
        //test of circumcircle   
        bool inside;
        if (ccw) {
            inside = fdapde::internals::in_circle(u->coords(), vw->node()->coords(), vw->twin()->node()->coords(), x->coords());
        } else {
            inside = fdapde::internals::in_circle(u->coords(), vw->twin()->node()->coords(), vw->node()->coords(), x->coords());
        }
        std::cout<<"inside: "<<inside<<std::endl;
        if (inside) {    //test fails so append vw to D and expand the cavity 
            vw->twin()->cell()->mark_visited();
            D.push_back(vw);
            mark_cavity(u, vw->twin()->prev(), D, C);
            mark_cavity(u, vw->twin()->next(), D, C);
        } else {
            C.push_back(vw);  // Add the new triangle to the list (without actually adding it)
            std::cout<<"SONO QUA CON : "<<vw->id()<<std::endl;
            return;
        }
    }

    // Function to insert a vertex handling conflicts
    void insert_vertex_at_conflict(node_t* u) {
        // Retrieve the triangle in conflict with u and we marked as visited 
        cell_t* vwx = u->conflict(); 
        vwx->mark_visited();

        std::vector<halfedge_t*> D;
        std::vector<halfedge_t*> C;

        mark_cavity(u, vwx->halfedge(), D, C);
        mark_cavity(u, vwx->halfedge()->prev(), D, C);
        mark_cavity(u, vwx->halfedge()->next(), D, C);
        
        //passing the conflicts of cthe cavity to a temporary vector 
        std::vector<node_t*> conflict_points_temp; 
        for (halfedge_t* h : D) { 
            cell_t* current_cell = h->cell();
            if (current_cell) {
                auto& conflict_list = current_cell->conflicting_points();
                for (node_t* point : conflict_list) {
                    if (point != u) { 
                        conflict_points_temp.push_back(point);
                        // checking if the point is in the cavity or is in a cell that does not belong to the cavity 
                        if (point->conflict() == current_cell) {
                            point->set_valid_conflict(false); // Invalidating the conflict 
                        }
                    }
                }
                current_cell->clear_conflicts();
            }


            cell_t* twin_cell = h->twin()->cell();
            if (twin_cell) {
                auto& conflict_list = twin_cell->conflicting_points();
                for (node_t* point : conflict_list) {
                    if (point != u) { 
                        conflict_points_temp.push_back(point);
                        // checking if the point is in the cavity or is in a cell that does not belong to the cavity 
                        if (point->conflict() == twin_cell) {
                            point->set_valid_conflict(false); // Invalidating the conflict 
                        }
                    }
                }
                twin_cell->clear_conflicts();
            }
        }
        
        for(node_t* n:conflict_points_temp)
        std::cout<<"nodo nel conflitto: "<<n->id()<<std::endl;
        // removing cells of the cavity 
        for (halfedge_t* h : D) { 
         //   std::cout << h->id() << std::endl;
            dcel_.remove_edge(h);
        }
        // creating the new cells 
        for (halfedge_t* h : C) { 
          //  std::cout << h->id() << std::endl;
            add_triangle(h, u->coords().transpose());
        }
        
        //reassing the conflicts to the new cells 
        for (node_t* y : conflict_points_temp) {
            for (halfedge_t* h : C) { // Ciclyng on the new cells
                cell_t* t = h->cell();
                if (!t) continue;
                const coords_t& t1 = t->halfedge()->prev()->node()->coords();
                const coords_t& t2 = t->halfedge()->node()->coords();
                const coords_t& t3 = t->halfedge()->next()->node()->coords();
    
                bool ccw = fdapde::internals::are_2d_counterclockwise_sorted(t1, t2, t3);
                // Test 1: Verifing if the point is inside the triangle
                //PROBLEMA CONTINUA A TROVARLO ANCHE DOPO AVERLO TROVATO GEOMETRICAMENTE 
                bool found = false;
                if (!y->is_valid_conflict() && !found){
                    bool inside_triangle = ccw ? fdapde::internals::point_in_2d_tri(y->coords(), t1, t2, t3)
                                            : fdapde::internals::point_in_2d_tri(y->coords(), t3, t2, t1);
                    // Assigning principal conflict                            
                    if (inside_triangle) {
                    y->set_conflict(t);
                    found = true;
                    } 
                } 
                // Test 2: Verifing if the point is in the circumcircle 
                bool inside_circumcircle = ccw ? fdapde::internals::in_circle(t1, t2, t3, y->coords())
                                               : fdapde::internals::in_circle(t3, t2, t1, y->coords());
                // adding n to list of conflict with t
                if (inside_circumcircle) t->add_conflict(y);  
            }
        }
    }   

    
    //overloaded one if user wants to impose internal points manually 
    void build_triangulation() {
    //BISOGNA UNIRE I CHECK PUNTO CASCA SUL LATO DELLA FIND_TRIANGLE CHE ORA NON ESISTE PIU
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
    // flip the initial trinagulation if not Delaunay 
        flip();
    
    //costruction of the conflict graph
    std::vector<node_t*> internal_nodes_to_insert; 

    for (coords_t& y : internal_points_) {
        if (y == internal_points_[0]) continue;  //already inserted
    
        node_t* n = dcel_.insert_node(node_t(dcel_.n_nodes(), /* boundary = */ false, y.transpose()));
        internal_nodes_to_insert.push_back(n);
        // finding the conflicts with existing cells 
        for (auto it = dcel_.cells_begin(); it != dcel_.cells_end(); ++it) {
            cell_t* t = &(*it);
            const coords_t& t1 = t->halfedge()->prev()->node()->coords();
            const coords_t& t2 = t->halfedge()->node()->coords();
            const coords_t& t3 = t->halfedge()->next()->node()->coords();

            bool ccw = fdapde::internals::are_2d_counterclockwise_sorted(t1, t2, t3);
            // Test 1: Verifing if the point is inside the triangle
            bool found = false;
            if(!found){
                bool inside_triangle = ccw ? fdapde::internals::point_in_2d_tri(y, t1, t2, t3)
                                        : fdapde::internals::point_in_2d_tri(y, t3, t2, t1);
                // Assigning principal conflict                            
                if (inside_triangle) {
                    n->set_conflict(t);  
                    found = true;
                }
            }
            // Test 2: Verifing if the point is in the circumcircle 
            bool inside_circumcircle = ccw ? fdapde::internals::in_circle(t1, t2, t3, y)
                                           : fdapde::internals::in_circle(t3, t2, t1, y);
            // adding n to list of conflict with t
            if (inside_circumcircle) t->add_conflict(n);  
        }
    }
    //inserting the num_points-1 inner nodes in the domain
    for (node_t* u : internal_nodes_to_insert) {
        insert_vertex_at_conflict(u);
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
/*
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
         //       std::cout << "nuovo punto inserito: " << u.transpose() << std::endl;
         //           internal_points_.push_back(u);}
              //  }
          //  }
/*        }
    
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
    }*/

    //////////////////////////// END OF NEW STORY /////////////////////////////////

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