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
    : dcel_(holes.empty() ? DCEL<local_dim, embed_dim>::make_polygon(boundary)
                    : DCEL<local_dim, embed_dim>::make_polygon(boundary, holes)) { }

    // Getter
    const DCEL<local_dim, embed_dim>& dcel() const { return dcel_; }
    DCEL<local_dim, embed_dim>& dcel() { return dcel_; }
    
    // Setter 
    //if the user wants to pass in input a triangulation with flip in order to make it Delaunay
    void set_dcel(const DCEL<local_dim, embed_dim>& dcel) {
        dcel_ = dcel;
        flip();
    }
    halfedge_t* add_triangle(halfedge_t* v, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& node){
        return dcel_.add_polygon(v ,node);
    }

    void add_first_triangle(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& point,
                            const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& boundary){   
        bool concave=false;
        auto iter = dcel_.halfedges_begin();

        // iterate over all boundary edges to connect to node, if possible
        for (int il = 0; il < boundary.rows(); ++il, ++iter) { 
            halfedge_t* v = &(*iter);                                                                  
            cell_t* c= v->cell();                                                                                    
            auto& cell_begin = *dcel_.cells_begin();
            std::vector<halfedge_t*> ghost_halfedges(3); 
        
            // add nodes and create ghost halfedges
            node_t* n;
            if(dcel_.find_node(point.row(0))==nullptr){  
                    n = dcel_.insert_node(node_t(dcel_.n_nodes(), false, point.row(0)));
                }
            else {  // node is already a vertex of the mesh
                    n = dcel_.find_node(point.row(0));  
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
                if(!ghost_halfedges[1]){
                    dcel_.insert_edge(dcel_.find_halfedge(v->prev()->prev()->node(), c),v);  //NON SICURISSIMA
                    return;
                }
            }
            
            // add edges
            for (int i = 0; i < 3 ; ++i) {
                halfedge_t* h1 = ghost_halfedges[i];
                halfedge_t* h2 = ghost_halfedges[(i + 1) % (3)];
                
                // check for intersection with boundary edges
                coords_t A = h1->node()->coords();
                coords_t B = h2->node()->coords();
                coords_t C;
                coords_t D;
                node_t* node_C;
                node_t* node_D;
                bool intersect=false;
                if(!(h1->on_boundary() && h2->on_boundary())) {  //if the edges are not both on the boundary
                    for (size_t i = 0; i < boundary.rows(); ++i) {  
                        C = boundary.row(i).transpose();
                        D = boundary.row((i + 1) % boundary.rows()).transpose();  
                        if (A != C && A != D && B != C && B != D && fdapde::internals::intersect(A, B, C, D)) {  
                            intersect = true;
                            concave=true;
                            break; 
                        }
                    }
                }
                if(!intersect){
                    halfedge_t* h = dcel_.insert_edge(h1, h2); 
                    if(h && h!=h1)
                        ghost_halfedges[(i + 1) % (3)] = h->next();  
                } 
            }
        }
        
        if(concave){
            // cycle over boundary edges only for safe handling of edges to insert
            for(auto it=dcel_.halfedges_begin(); it!=dcel_.halfedges_end() && (&(*it))->cell(); ++it){
                halfedge_t* h = &(*it);
                int i=0;
                do{
                    h=h->next();
                    i++;
                }while(h!=&(*it));
                if(i>3) //not a triangle yet
                    for(int l=0; l<i-3; ++l){
                      halfedge_t* n=h->next()->next();
                      if(h->node()->on_boundary() && h->next()->next()->node()->on_boundary() && !fdapde::internals::collinear(h->node()->coords(),h->next()->node()->coords(),h->next()->next()->node()->coords()))
                        dcel_.insert_edge(h, n);
                      h=n;   
                    }
            }
        }
        flip();
    }


    void flip() {

        // creating a list of halfedges to check whether they are locally delaunay or not (in this case flippable)
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
                //std::cout << "FLIP" << std::endl;
            
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


    // function to mark the cavity during insertion
    void mark_cavity(node_t* u, halfedge_t* h, std::vector<halfedge_t*>& D, std::vector<halfedge_t*>& C) {
        if(h->on_boundary()){
            C.push_back(h);  
            return;
        }

        node_t* x = dcel_.adjacent(h);
        if (!x) {
            return;
        } 
        bool ccw = fdapde::internals::are_2d_counterclockwise_sorted(u->coords(), h->node()->coords(), h->twin()->node()->coords());
        //test of circumcircle   
        bool inside;
        if (ccw) {
            inside = fdapde::internals::in_circle(u->coords(), h->node()->coords(), h->twin()->node()->coords(), x->coords());
        } else {
            inside = fdapde::internals::in_circle(u->coords(), h->twin()->node()->coords(), h->node()->coords(), x->coords());
        }
        if (inside) {    //test fails so append vw to D and expand the cavity 
            D.push_back(h);
            mark_cavity(u, h->twin()->prev(), D, C);
            mark_cavity(u, h->twin()->next(), D, C);
            return;
        } else {
            C.push_back(h);  // Add the new triangle to the list (without actually adding it)
            return;
        }
    }


    // Function to insert a vertex handling conflicts
    void insert_vertex_at_conflict(node_t* u) {
        // Retrieve the triangle in conflict with u and we marked as visited 
        cell_t* t = u->conflict(); 
        std::vector<halfedge_t*> D;
        std::vector<halfedge_t*> C;

        mark_cavity(u, t->halfedge(), D, C);
        mark_cavity(u, t->halfedge()->prev(), D, C);
        mark_cavity(u, t->halfedge()->next(), D, C);
  
        //invalidating the conflicts node->cell for the point of the cavity
        for (halfedge_t* h : D) {
            cell_t* current_cell = h->cell();
            if (current_cell) {
                for (node_t* point : current_cell->conflicting_points()) {
                    if (point != u && point->conflict() == current_cell) 
                        point->set_conflict(nullptr);  
                }
                current_cell->clear_conflicts();
            }
        
            cell_t* twin_cell = h->twin()->cell();
            if (twin_cell) {
                for (node_t* point : twin_cell->conflicting_points()) {
                    if (point != u && point->conflict() == twin_cell) 
                        point->set_conflict(nullptr); 
                }
                twin_cell->clear_conflicts();
            }
        }
        //if all the new traingles are Delaunay and i am not expanding the cavity 
        if(D.empty()){
            cell_t* current_cell = u->conflict();
            if (current_cell) {
                for (node_t* point : current_cell->conflicting_points()) {
                    if (point != u && point->conflict() == current_cell) 
                        point->set_conflict(nullptr);
                }
                current_cell->clear_conflicts();
            }
        }
        u->remove_conflict();
        
        // removing cells of the cavity 
        for (halfedge_t* h : D) { 
            dcel_.remove_edge(h);
        }
        // creating the new cells 
        for (halfedge_t* h : C) { 
            add_triangle(h, u->coords().transpose());
        }
        
        // Reassigning the conflicts to the new cells  
        for (auto it = dcel_.nodes_begin(); it != dcel_.nodes_end(); ++it) {
            node_t* y = &(*it);
            if (y == u || y->on_boundary() || y->is_inserted())
                continue;

            detect_conflicts(y, C);  
        }
        u->set_inserted(true);
    }  

    void build_triangulation(int N, const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& boundary) {

        // Computing the bounding box
        double min_x = boundary.col(0).minCoeff();
        double max_x = boundary.col(0).maxCoeff();
        double min_y = boundary.col(1).minCoeff();
        double max_y = boundary.col(1).maxCoeff();
    
        // Creating the generator of causal numbers
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dist_x(min_x, max_x);
        std::uniform_real_distribution<double> dist_y(min_y, max_y);
   
        int first_internal_id = -1;
        int generated_points = 0;

        while (generated_points < N) {
            coords_t u;
            u << dist_x(gen), dist_y(gen);
            // Verifiyng if the point is inside the polygon, in order to coorecty dale with concavities
            if (!fdapde::internals::point_in_polygon(boundary, u)) 
                continue; 
            if (generated_points == 0) {
                // Initialize triangulation with the first valid point
                add_first_triangle(u.transpose(), boundary);
                first_internal_id = dcel_.n_nodes() - 1;
            } else {
                // Create the node and detect the conflicts with existing cells 
                node_t* n = dcel_.insert_node(node_t(dcel_.n_nodes(), false, u));
                detect_conflicts(n); 
            }
            ++generated_points;
        }
        // inserting the remaining nodes in the domain 
        for (auto it = dcel_.nodes_begin(); it != dcel_.nodes_end(); ++it) {
            node_t* u = &(*it);
            if (!u->on_boundary() && u->id()!=first_internal_id)   
                insert_vertex_at_conflict(u); 
        }
         /*   
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
        */
        /*mesh generation not included in order to test the efficiency of the triangulation
        Triangulation<local_dim, embed_dim> triangulation = DCEL_to_Triangulation();

        std::string filename = "mesh_output.txt";
        export_triangulation_to_txt(triangulation, filename);

        std::string command = "python3 fdaPDE/src/plot_mesh.py";
        std::system(command.c_str()); 
        */
    }

    //overloaded one if user wants to pass manually the internal points
    //the user must know the passed internal points lie all inside the domain 
    void build_triangulation(const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& internal,
        const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& boundary) {

        add_first_triangle(internal.row(0), boundary);
        int first_internal_id = dcel_.n_nodes() - 1;

        //inserting the remaining internal points in the triangulation
        for (int i = 1; i < internal.rows(); ++i) {
            node_t* n = dcel_.insert_node(node_t(dcel_.n_nodes(), false, internal.row(i).transpose().eval()));
            //detecting the conflicts with existing cells
            detect_conflicts(n);
        }
        //inserting the remaining nodes in the domain
        for (auto it = dcel_.nodes_begin(); it != dcel_.nodes_end(); ++it) {
            node_t* u = &(*it);
            if (!u->on_boundary() && u->id()!=first_internal_id)  
                insert_vertex_at_conflict(u); 
        }
        
        /*   
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
        */
        /*mesh generation not included in order to test the efficiency of the triangulation
        Triangulation<local_dim, embed_dim> triangulation = DCEL_to_Triangulation();

        std::string filename = "mesh_output.txt";
        export_triangulation_to_txt(triangulation, filename);

        std::string command = "python3 fdaPDE/src/plot_mesh.py";
        std::system(command.c_str()); 
        */
    }

    void detect_conflicts(node_t* n, const std::vector<halfedge_t*>& cells_to_check = {}) {
        bool found = false;
    
        if (cells_to_check.empty()) {  // initialization case: scan all the cells
            for (auto it = dcel_.cells_begin(); it != dcel_.cells_end(); ++it) {
                cell_t* t = &(*it);
                const coords_t& t1 = t->halfedge()->prev()->node()->coords();
                const coords_t& t2 = t->halfedge()->node()->coords();
                const coords_t& t3 = t->halfedge()->next()->node()->coords();
    
                bool ccw = fdapde::internals::are_2d_counterclockwise_sorted(t1, t2, t3);
    
                // Test 1: Verifying if the point is inside the triangle
                if (!found) {
                    bool inside_triangle = ccw ? fdapde::internals::point_in_2d_tri(n->coords(), t1, t2, t3)
                                               : fdapde::internals::point_in_2d_tri(n->coords(), t3, t2, t1);
                    if (inside_triangle) {
                        n->set_conflict(t); 
                        found = true;
                    }
                }
    
                // Test 2: Verifying if the point is in the circumcircle 
                bool inside_circumcircle = ccw ? fdapde::internals::in_circle(t1, t2, t3, n->coords())
                                               : fdapde::internals::in_circle(t3, t2, t1, n->coords());
    
                if (inside_circumcircle) t->add_conflict(n);  
            }
        } else {  // normal case: scan the cavity
            for (halfedge_t* h : cells_to_check) {
                cell_t* t = h->cell();
                if (!t) continue;
    
                const coords_t& t1 = t->halfedge()->prev()->node()->coords();
                const coords_t& t2 = t->halfedge()->node()->coords();
                const coords_t& t3 = t->halfedge()->next()->node()->coords();
    
                bool ccw = fdapde::internals::are_2d_counterclockwise_sorted(t1, t2, t3);
    
                if (!found) {
                    bool inside_triangle = ccw ? fdapde::internals::point_in_2d_tri(n->coords(), t1, t2, t3)
                                               : fdapde::internals::point_in_2d_tri(n->coords(), t3, t2, t1);
                    if (inside_triangle) {
                        n->set_conflict(t); 
                        found = true;
                    }
                }
    
                bool inside_circumcircle = ccw ? fdapde::internals::in_circle(t1, t2, t3, n->coords())
                                               : fdapde::internals::in_circle(t3, t2, t1, n->coords());
    
                if (inside_circumcircle) t->add_conflict(n);  
            }
        }
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
};
  
}  // namespace fdapde

#endif // __DELAUNAY_H__