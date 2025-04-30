#ifndef __FDAPDE_DELAUNAY_H__
#define __FDAPDE_DELAUNAY_H__

//Delaunay class build a triangulation and refine (using Ruppert algorithm) a generic connected domain 
//including concavities 
#include "header_check.h"
namespace fdapde {
  
template <int LocalDim, int EmbedDim>
class Delaunay : public TriangulationBase<LocalDim, EmbedDim, Triangulation<2,2>> {
   public:
    static constexpr int local_dim = LocalDim;
    static constexpr int embed_dim = EmbedDim;
    static constexpr int n_nodes_cell = 3;

    using coords_t = Eigen::Matrix<double, embed_dim, 1>;
    using node_t = typename DCEL<local_dim, embed_dim>::node_t;
    using halfedge_t = typename DCEL<local_dim, embed_dim>::halfedge_t;
    using cell_t = typename DCEL<local_dim, embed_dim>::cell_t;
    using triangulation_t = TriangulationBase<local_dim, embed_dim, Triangulation<2,2>>;
    using dcel_t = DCEL<local_dim, embed_dim>;
    using polygon_t = Polygon<local_dim, embed_dim>;

    //decorator pattern from dcel.h for Conflict Graph algorithm (not working)
    /*
    using base_node_t = typename DCEL<local_dim, embed_dim>::node_t;
    using base_cell_t = typename DCEL<local_dim, embed_dim>::cell_t;
    struct node_t;
    struct cell_t;
    struct node_t : public base_node_t{
        using base_node_t::base_node_t;
        cell_t* conflicting_triangle_=nullptr;
        void set_conflict(cell_t* triangle) { conflicting_triangle_ = triangle; }
        cell_t* conflict() const { return conflicting_triangle_; }
        void remove_conflict() { conflicting_triangle_ = nullptr; }
    };
    struct cell_t : public DCEL<local_dim, embed_dim>::cell_t{
        using base_cell_t::base_cell_t;
        std::vector<node_t*> conflicting_points_;
        void add_conflict(node_t* point) { conflicting_points_.push_back(point); }
        std::vector<node_t>& conflicting_points() const{ return conflicting_points_; }
        std::vector<node_t*>& conflicting_points() { return conflicting_points_; }
        void clear_conflicts() { conflicting_points_.clear(); }
    };
    */
    
    
    // constructors 
    // costructor with random generated points
    Delaunay(const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& boundary, int N = 100) :
        triangulation_t(build_triangulation_from_dcel(boundary, N)) {}
    // costructor with given internal points form the user
    // user needs to provide internal points correctly located inside the domain 
    Delaunay(const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& boundary, const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& internal) :
        triangulation_t(build_triangulation_from_dcel(boundary, internal)) {}

    //function running the refinment with Ruppert algorithm 
    //for the moment it can work only from the inside with default parameters in the constructor
    static void Ruppert_refinement(dcel_t& dcel, double rho_bar) {
        //set needed to save in memory edges that encroach a point in the triangulation and another one for the badly shaped traingles 
        std::unordered_set<halfedge_t*> encroached_edges;
        std::unordered_set<cell_t*> bad_triangles;     
    
        // inizialization of the two set
        for (auto it = dcel.halfedges_begin(); it != dcel.halfedges_end(); ++it) {
            halfedge_t* e = &(*it);
            //since in the previuos code we do not touch the boundary we can stop at the twin of the first boundary edge
            if(e->on_boundary() && !e->cell()) break;
            if (e->on_boundary() && check_encroachment(dcel, e)){
                //to make sure the set saves only one copy of the edge
                if (encroached_edges.count(e) == 0)
                    encroached_edges.insert(e);
            }
        }
        for (auto it = dcel.cells_begin(); it != dcel.cells_end(); ++it) {
            cell_t* t = &(*it);
            if (is_bad_triangle(dcel, t, rho_bar)) {
                if (bad_triangles.count(t) == 0)
                    bad_triangles.insert(t);
            }
        }  
        //while keeps running until the two set are empty and every time a bad triangle is exiting the triangulation
        //the test of the encroached edges runs in order to keep track of the newly created triangulation
        while (true) {
            if (split_first_encroached_segment(dcel, encroached_edges, bad_triangles, rho_bar)) {
                continue;
            }
            if (split_first_bad_triangle(dcel, rho_bar, encroached_edges, bad_triangles)) {
                continue;
            }
            break;
        }
        // final reorder to cut no longer existing id of cells and edges
        int cont = 0;
        for (auto it = dcel.cells_begin(); it != dcel.cells_end(); ++it)
            it->set_id(cont++);
        dcel.set_n_cells_(cont);
    
        cont = 0;
        for (auto it = dcel.halfedges_begin(); it != dcel.halfedges_end(); ++it)
            it->set_id(cont++);
        
        //json needed for debug
        dcel.export_to_json("dcel_output.json");
    }
    
    
   private:

    // function performing a Delaunay triangulation with N uniformly random points 
    static triangulation_t build_triangulation_from_dcel(const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& boundary, int N) 
    {
        dcel_t dcel; 
        triangulate(dcel, N, boundary); 

        //clock in order to test the performance of Ruppert refinement
        auto start = high_resolution_clock::now();  
        Ruppert_refinement(dcel, 3.0);
        auto end = high_resolution_clock::now();   
        auto duration = duration_cast<milliseconds>(end - start).count();
        std::cout << "elapsed REFINMENT time for " << dcel.n_nodes() << " points: " << duration << " ms" << std::endl;
       
        //converting the dcel in triangulation 
        return DCEL_to_Triangulation(dcel);
    }

    //overloaded function of above in order to work with user defined internal points 
    static triangulation_t build_triangulation_from_dcel(const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& boundary, const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& internal)
    {
        dcel_t dcel; 
        triangulate(dcel, internal, boundary);  
        Ruppert_refinement(dcel, 3.0);
        return DCEL_to_Triangulation(dcel);
    }

    //function to convert a dcel into a base triangulation
    static triangulation_t DCEL_to_Triangulation(dcel_t& dcel) {
        Eigen::Matrix<double, Eigen::Dynamic, embed_dim> nodes(dcel.n_nodes(), embed_dim);
        Eigen::Matrix<int, Eigen::Dynamic, n_nodes_cell> cells(dcel.n_cells(), n_nodes_cell);
        Eigen::Matrix<int, Eigen::Dynamic, 1> boundary_markers(dcel.n_nodes());

        // Fill the nodes matrix and the boundary marker vector
        int node_idx = 0;
        for (auto it = dcel.nodes_begin(); it != dcel.nodes_end(); ++it) {
            nodes.row(node_idx) = it->coords().transpose();
            boundary_markers(node_idx) = it->on_boundary() ? 1 : 0;
            node_idx++;
        }
        // Fill the cells matrix by extracting the three nodes of each triangle from the halfedge structure
        int cell_idx = 0;
        for (auto it = dcel.cells_begin(); it != dcel.cells_end(); ++it) {
            halfedge_t* h = it->halfedge();
            for (int i = 0; i < n_nodes_cell; ++i) {
                cells(cell_idx, i) = h->node()->id();
                h = h->next();
            }
            cell_idx++;
        }
        
        triangulation_t triangulation(nodes, cells, boundary_markers);
        //lines used to graph the trinagulation
        /*std::string filename = "mesh_output.txt";
        export_triangulation_to_txt(triangulation, filename);
        std::string command = "python3 fdaPDE/src/plot_mesh.py";
        std::system(command.c_str());*/

        return triangulation;
    }


    // Function to convert a base triangulation into a DCEL structure
    static dcel_t Triangulation_to_DCEL(triangulation_t& triangulation) {
        using MatrixType = Eigen::Matrix<double, Eigen::Dynamic, embed_dim>;

        // Get node coordinates and boundary markers from the triangulation
        const auto& coords = triangulation.nodes();
        const auto& markers = triangulation.boundary_nodes();

        // Extract only the boundary nodes
        int n_boundary = triangulation.n_boundary_nodes();  
        MatrixType boundary_nodes(n_boundary, embed_dim);

        int idx = 0;
        for (int i = 0; i < coords.rows(); ++i) {
            if (markers(i, 0) == 1) {
                boundary_nodes.row(idx++) = coords.row(i);
            }
        }

        // Create the initial DCEL polygon using only boundary nodes
        dcel_t dcel = dcel_t::make_polygon(boundary_nodes);

        // Insert internal (non-boundary) nodes into the DCEL
        for (int i = 0; i < coords.rows(); ++i) {
            if (markers(i, 0) == 0) {
                dcel.insert_node(node_t(dcel.n_nodes(), false, coords.row(i)));
            }
        }

        // Reconstruct cells (triangles) by identifying node pointers and adding them to the DCEL
        for (int i = 0; i < triangulation.n_cells(); ++i) {
            int id0 = triangulation.cells()(i, 0);
            int id1 = triangulation.cells()(i, 1);
            int id2 = triangulation.cells()(i, 2);

            // Find the node pointers in the DCEL corresponding to the triangle vertices
            node_t* n0 = std::addressof(*std::find_if(dcel.nodes_begin(), dcel.nodes_end(),
                [&](const node_t& n) { return n.id() == id0; }));

            node_t* n1 = std::addressof(*std::find_if(dcel.nodes_begin(), dcel.nodes_end(),
                [&](const node_t& n) { return n.id() == id1; }));

            node_t* n2 = std::addressof(*std::find_if(dcel.nodes_begin(), dcel.nodes_end(),
                [&](const node_t& n) { return n.id() == id2; }));

            // Ensure the triangle orientation is counterclockwise
            if (!fdapde::internals::are_2d_counterclockwise_sorted(n0->coords(), n1->coords(), n2->coords())) {
                std::swap(n1, n2);
            }

            halfedge_t* h = nullptr;

            // Try to find a suitable half-edge to attach the triangle
            h = dcel.find_halfedge_between(n0, n1);  
            if (h) {
                add_triangle(dcel, h, std::vector<node_t*> {n2});
            } else {
                h = dcel.find_halfedge_between(n1, n2);
                if (h) {
                    add_triangle(dcel, h, std::vector<node_t*> {n0});
                } else {
                    h = dcel.find_halfedge_between(n2, n0);
                    if (h) {
                        add_triangle(dcel, h, std::vector<node_t*> {n1});
                    } else {
                        // If no edge exists, create a new one and add the triangle
                        halfedge_t* h0 = dcel.emplace_halfedge_(n0);
                        halfedge_t* h1 = dcel.emplace_halfedge_(n1);
                        dcel.insert_edge(h0, h1);
                        add_triangle(dcel, h0, std::vector<node_t*> {n2});
                    }
                }
            }
        }

        return dcel;
    }

    //key function for the implementation of the conflict algorithm 
    static void detect_conflicts(dcel_t& dcel ,node_t* n, const std::vector<halfedge_t*>& cells_to_check = {}) {
        bool found = false;
        // initialization case: scan all the cells
        if (cells_to_check.empty()) {  
            for (auto it = dcel.cells_begin(); it != dcel.cells_end(); ++it) {
                cell_t* t = &(*it);
                const coords_t& t1 = t->halfedge()->prev()->node()->coords();
                const coords_t& t2 = t->halfedge()->node()->coords();
                const coords_t& t3 = t->halfedge()->next()->node()->coords();

                bool ccw = fdapde::internals::are_2d_counterclockwise_sorted(t1, t2, t3);

                // Test : Verifying if the point is inside the triangle
                if (!found) {
                    bool inside_triangle = ccw ? fdapde::internals::point_in_2d_tri(n->coords(), t1, t2, t3)
                                            : fdapde::internals::point_in_2d_tri(n->coords(), t3, t2, t1);
                    if (inside_triangle) {
                        //the conflict is from point to triangle and viceversa
                        n->set_conflict(t); 
                        t->add_conflict(n);
                        found = true;
                    }
                }
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
                        t->add_conflict(n);
                        found = true;
                    }
                }
            }
        }
    }

    static halfedge_t* add_triangle(dcel_t& dcel ,halfedge_t* v,const std::vector<node_t*>& node){
        return dcel.add_polygon(v ,node);
    }

    //function performing the first raw triangulation of the domain using the polygon.h class
    static void initialize_triangulation(dcel_t& dcel, const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& boundary){
        polygon_t polygon(boundary);
        auto triangulation = polygon.triangulation();
        const auto& nodes = triangulation.nodes();  
        const auto& cells = triangulation.cells();   
        //we manage our input dcel to be the traslation of the trinagulation object from polygon 
        dcel = Triangulation_to_DCEL(triangulation);
    }

    //function perfoming the flip alghoritm to tranform every non-Delaunay triangulation into a Delaunay one
    //not used for the actual costruction for the O(n^2) complexity
    static void flip(dcel_t& dcel) {

        // creating a list of halfedges to check whether they are locally delaunay or not (in this case flippable)
        std::list<halfedge_t*> halfedges_to_check;
        for (auto it = dcel.halfedges_begin(); it != dcel.halfedges_end(); ++it,++it) {
            //we can already exclude the edges of the triangulation being automatically locally delaunay
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
            //testing if D lies inside the circumcircle of triangle ABC or C lies inside the circumcircle of triangle ABD
            if (fdapde::internals::in_circle(A, B, C, D) || fdapde::internals::in_circle(A, D, B, C)) {
                
                // we flip since edge is not locally delaunay
                halfedge_t* e = edge;
                dcel.remove_edge(edge);
                halfedge_t* new_edge = dcel.insert_edge(e->prev(), e->twin()->prev());
                //std::cout << "FLIP" << std::endl;
            
                if (new_edge) {
                    // Inserting the new halfedges created by the flip into the list to check
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

    // Function to insert a vertex handling conflicts
    static void insert_vertex_at_conflict(dcel_t& dcel, node_t* u) {
        // Retrieve the triangle in conflict with u 
        cell_t* t = u->conflict(); 
        // vector storing the halfedge refering to the cells that need to be deleted(D) and created(C)
        std::vector<halfedge_t*> D;
        std::vector<halfedge_t*> C;

        //removing the edge if the point falls on it 
        const coords_t& t1 = t->halfedge()->node()->coords();
        const coords_t& t2 = t->halfedge()->next()->node()->coords();
        const coords_t& t3 = t->halfedge()->prev()->node()->coords();
        
        bool found_on_edge = false;

        if (fdapde::internals::contains(u->coords(), t1, t2)) {
            D.push_back(t->halfedge());
                //we manually work on the remaining cavity in order to correctly activate the algorithm 
                mark_cavity(dcel, u, t->halfedge()->next(), D, C);
                mark_cavity(dcel, u, t->halfedge()->prev(), D, C);
                mark_cavity(dcel, u, t->halfedge()->twin()->next(), D, C);
                mark_cavity(dcel, u, t->halfedge()->twin()->prev(), D, C);
            found_on_edge = true;
        }
    
        else if (fdapde::internals::contains(u->coords(), t2, t3)) {
            D.push_back(t->halfedge()->next());
                mark_cavity(dcel, u, t->halfedge(), D, C);
                mark_cavity(dcel, u, t->halfedge()->prev(), D, C);
                mark_cavity(dcel, u, t->halfedge()->next()->twin()->next(), D, C);
                mark_cavity(dcel, u, t->halfedge()->next()->twin()->prev(), D, C);
            found_on_edge = true;
        }
    
        else if (fdapde::internals::contains(u->coords(), t3, t1)) {
            D.push_back(t->halfedge()->prev());
                mark_cavity(dcel, u, t->halfedge()->next(), D, C);
                mark_cavity(dcel, u, t->halfedge(), D, C);
                mark_cavity(dcel, u, t->halfedge()->prev()->twin()->next(), D, C);
                mark_cavity(dcel, u, t->halfedge()->prev()->twin()->prev(), D, C);
            found_on_edge = true;
        }
        //if the point did not fall on any edge then we normally try to expand the cavity from the edge of the triangle 
        if (!found_on_edge) {
            mark_cavity(dcel, u, t->halfedge(), D, C);
            mark_cavity(dcel, u, t->halfedge()->next(), D, C);
            mark_cavity(dcel, u, t->halfedge()->prev(), D, C);
        }
    
        //invalidating the conflicts node->cell for the point of the cavity
        std::unordered_set<node_t*> invalidated_nodes;
        for (halfedge_t* h : D) {
            cell_t* current_cell = h->cell();
            if (current_cell) {
                for (node_t* point : current_cell->conflicting_points()) {
                    if (point != u) {
                        point->set_conflict(nullptr);  
                        invalidated_nodes.insert(point);
                    }
                }
                current_cell->clear_conflicts();
            }
            //this must be done also for the cell on the other side of the halfedge 
            cell_t* twin_cell = h->twin()->cell();
            if (twin_cell) {
                for (node_t* point : twin_cell->conflicting_points()) {
                    if (point != u) {
                        point->set_conflict(nullptr);  
                        invalidated_nodes.insert(point);
                    }
                }
                twin_cell->clear_conflicts();
            }
        }
        //if all the new traingles are Delaunay and i am not expanding the cavity 
        if(D.empty()){
            cell_t* current_cell = u->conflict();
            if (current_cell) {
                for (node_t* point : current_cell->conflicting_points()) {
                    if (point != u) {
                        point->set_conflict(nullptr);  
                        invalidated_nodes.insert(point);
                    }
                }
                current_cell->clear_conflicts();
            }
        }
        u->remove_conflict();
        
        // removing cells of the cavity 
        for (halfedge_t* h : D) { 
            dcel.remove_edge(h);
        }
        // creating the new cells 
        for (halfedge_t* h : C) { 
            add_triangle(dcel, h, std::vector<node_t*> {u});
        }
        
        // Reassigning the conflicts to the new cells  
        for (node_t* y : invalidated_nodes) {
            detect_conflicts(dcel, y, C);
        }
    } 

    // function to mark the cavity during insertion
    static void mark_cavity(dcel_t& dcel, node_t* u, halfedge_t* h, std::vector<halfedge_t*>& D, std::vector<halfedge_t*>& C) {
        //if we reach the boundary we automatically create the triangle
        if(h->on_boundary()){
            C.push_back(h);  
            return;
        }

        node_t* x = dcel.adjacent(h);
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
            mark_cavity(dcel, u, h->twin()->prev(), D, C);
            mark_cavity(dcel, u, h->twin()->next(), D, C);
            return;
        } else {
            C.push_back(h);  // Add the new triangle to the list (without actually adding it)
            return;
        }
    }
    
    //trinagulating the domain in a Delaunay fashion with N random points 
    static void triangulate(dcel_t& dcel, int N, const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& boundary) {
        double min_x = boundary.col(0).minCoeff();
        double max_x = boundary.col(0).maxCoeff();
        double min_y = boundary.col(1).minCoeff();
        double max_y = boundary.col(1).maxCoeff();
        //creating the generator for the internal point of the triangulation
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dist_x(min_x, max_x);
        std::uniform_real_distribution<double> dist_y(min_y, max_y);
    
        //we build the first raw triangulation and flip it in order to gain a Delaunay one
        initialize_triangulation(dcel, boundary);
        flip(dcel);
        
        int generated_points = 0;

        while (generated_points < N) {
            coords_t u;
            u << dist_x(gen), dist_y(gen);
            //in case of concavities we fall out of the domain we discard the point
            if (!fdapde::internals::point_in_polygon(boundary, u))
            continue;

            node_t* n = dcel.insert_node(node_t(dcel.n_nodes(), false, u));
            detect_conflicts(dcel, n);
            ++generated_points;
            
        }
    
        for (auto it = dcel.nodes_begin(); it != dcel.nodes_end(); ++it) {
            node_t* u = &(*it);
            //creating the connections for every point not on the boundary
            if (!u->on_boundary())
                insert_vertex_at_conflict(dcel, u);
        }
        //riordering deleting the jumps of id's in halfedges and cells 
        int cont = 0;
        for (auto it = dcel.cells_begin(); it != dcel.cells_end(); ++it)
            it->set_id(cont++);
        dcel.set_n_cells_(cont);
    
        cont = 0;
        for (auto it = dcel.halfedges_begin(); it != dcel.halfedges_end(); ++it)
            it->set_id(cont++);
    }
    
    //overloaded one if user wants to pass manually the internal points
    //the user must know the passed internal points lie all inside the domain 
    static void triangulate(dcel_t& dcel, const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& internal,
        const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& boundary) {
        
        initialize_triangulation(dcel, boundary);
        flip(dcel);

        //inserting the internal points in the triangulation
        for (int i = 0; i < internal.rows(); ++i) {
            node_t* n = dcel.insert_node(node_t(dcel.n_nodes(), false, internal.row(i).transpose().eval()));
            detect_conflicts(dcel, n);
        }
        for (auto it = dcel.nodes_begin(); it != dcel.nodes_end(); ++it) {
            node_t* u = &(*it);
            if (!u->on_boundary()) 
                insert_vertex_at_conflict(dcel, u); 
        }   
        //reordering id of cells and halfedges to cover some jumps between ids after removing
        int cont = 0;
        for (auto it = dcel.cells_begin(); it != dcel.cells_end(); ++it) {
            it->set_id(cont);
            cont++;
        }

        dcel.set_n_cells_(cont);

        int cont_h = 0;
        for (auto it = dcel.halfedges_begin(); it != dcel.halfedges_end(); ++it) {
            it->set_id(cont_h);
            cont_h++;
        }
    }

//now follow a series of function needed to perform the Ruppert refinement algorithm
    //function neeeded to perform the test of encroachment for an edge e of the triangulation
    static bool check_encroachment(dcel_t& dcel, halfedge_t* e) {
        //extracting verticies of the adjacent triangle to the edge e 
        coords_t A = e->node()->coords();
        coords_t B = e->twin()->node()->coords();
        coords_t C = e->prev()->node()->coords();         
        //computing the opposite angle to the edge e 
        double angle = fdapde::internals::angle_between(B, C, A);
        return angle >= 90.0;
    }
    //function needed to perform the test of badly shaped triangle 
    //using the ratio betwwen the radius of the circumcircle and the longest edge
    //the convergence of the algorithm is proved for rho_bar >=sqrt(2)
    static bool is_bad_triangle(dcel_t& dcel, cell_t* t, double rho_bar) {
        coords_t A = t->halfedge()->prev()->node()->coords();
        coords_t B = t->halfedge()->node()->coords();
        coords_t C = t->halfedge()->next()->node()->coords();
        double ratio = fdapde::internals::radius_edge_ratio(A, B, C);
        return ratio > rho_bar;
    }
                               
    // Function to split an encroached subsegment 
    static void split_subsegment(dcel_t& dcel, halfedge_t* e, std::unordered_set<halfedge_t*>& encroached_edges, 
        std::unordered_set<cell_t*>& bad_triangles, double rho_bar) {
        // Extract the two endpoints of the edge e and the third vertex forming the adjacent triangle
        node_t* a = e->node();                
        node_t* b = e->twin()->node();      
        node_t* c = e->prev()->node();      
        coords_t split_pt;                    // the point where the segment will be split

        // Case 1: acute angle at vertex b and e->next is on boundary
        if (e->next()->on_boundary() && fdapde::internals::is_angle_acute(a->coords(), b->coords(), c->coords())) {
            coords_t split_pt_ref = c->coords();
            double r = (split_pt_ref - b->coords()).norm();  
            coords_t ab = a->coords() - b->coords();
            double L = ab.norm();                           
            // Ensure r does not exceed the length of the segment
            if (r > L) r = 0.5 * L;
            double t = r / L;
            // compute split point at distance r from b on segment ab
            split_pt = b->coords() + t * ab;
        } 
        // Case 2: regular midpoint split
        else {
            split_pt = 0.5 * (a->coords() + b->coords());     // midpoint of segment ab
        }
         // Insert the new node into the DCEL
        node_t* m = dcel.insert_node(node_t(dcel.n_nodes(), true, split_pt));
        halfedge_t* prev = e->prev();
        // Subdivide the triangle adjacent to edge e by inserting m
        add_triangle(dcel, e, std::vector<node_t*> {m});
        halfedge_t* h1 = e->next()->twin();  
        halfedge_t* h2 = e->prev()->twin();
        // Remove both affected triangles from the set of bad triangles
        bad_triangles.erase(e->cell());
        bad_triangles.erase(e->twin()->cell());
        // Remove the encroached edge from the DCEL
        dcel.remove_edge(e);
        // Subdivide the opposite triangle (twin) by connecting m
        add_triangle(dcel, prev, std::vector<node_t*> {m});
        // Perform local flips if necessary to maintain Delaunay property
        flip_Ruppert(dcel, encroached_edges, bad_triangles, rho_bar);

        // After inserting m, recheck the two new boundary-adjacent edges for possible encroachment
        if (h1->on_boundary() && h1->cell() && check_encroachment(dcel, h1)) {
            if (encroached_edges.count(e) == 0)
                encroached_edges.insert(h1);
        }
        if (h2->on_boundary() && h2->cell() && check_encroachment(dcel, h2)) {
            if (encroached_edges.count(e) == 0)
                encroached_edges.insert(h2);
        }
    }

    // Attempts to split a bad triangle by inserting its circumcenter.
    // If the circumcenter encroaches a segment of the PLC, it splits that segment instead.
    // Returns true if a refinement was performed.
    static bool split_triangle(dcel_t& dcel, cell_t* t, std::unordered_set<halfedge_t*>& encroached_edges,
        std::unordered_set<cell_t*>& bad_triangles, double rho_bar) {
        // Get triangle vertices A, B, C
        coords_t A = t->halfedge()->prev()->node()->coords();
        coords_t B = t->halfedge()->node()->coords();
        coords_t C = t->halfedge()->next()->node()->coords();

        // Compute the circumcenter of triangle ABC
        coords_t c = fdapde::internals::circumcenter(A, B, C);  

        // Check if a node already exists at the circumcenter (avoid duplicates)
        for (auto it = dcel.nodes_begin(); it != dcel.nodes_end(); ++it) {
            if ((it->coords() - c).norm() < 1e-12) {
                return false;  // Do not insert if a node is already at c
            }
        }

        // Check whether the circumcenter c encroaches any boundary segment
        for (auto it = dcel.halfedges_begin(); it != dcel.halfedges_end(); ++it) {
            halfedge_t* e = &(*it);

            // We only care about edges that are part of the PLC (on the boundary)
            if (e->on_boundary() && e->cell()) {
                coords_t a = e->node()->coords();
                coords_t b = e->twin()->node()->coords();

                // If the circumcenter c encroaches the boundary segment ab
                if (fdapde::internals::is_encroached(c, a, b)) {
                    // Only split if the edge and its neighborhood is not marked as "seditious"
                    if (!is_edge_seditious(e) && !is_edge_seditious(e->twin()) &&
                        !is_edge_seditious(e->next()) && !is_edge_seditious(e->prev()) &&
                        !is_edge_seditious(e->next()->twin()) && !is_edge_seditious(e->prev()->twin())) {
                        
                        // Split the encroached subsegment instead of inserting the circumcenter
                        split_subsegment(dcel, e, encroached_edges, bad_triangles, rho_bar);
                        return true;
                    }

                    // If the segment is seditious, do nothing now (will be retried later)
                    return false;
                }
            }
        }

        // If no encroachment is detected, insert the circumcenter into the mesh
        node_t* circ = dcel.insert_node(node_t(dcel.n_nodes(), false, c));
        // Locate the triangle containing the new point
        cell_t* cf = find_triangle(dcel, c); 
        // Insert the new node into the triangulation (splitting the containing triangle)
        insert_vertex(dcel, circ, cf, encroached_edges, bad_triangles, rho_bar);

        return true;
    }


    // Attempts to split the first non-seditious encroached edge.
    // Returns true if a segment was successfully split.
    static bool split_first_encroached_segment(dcel_t& dcel, std::unordered_set<halfedge_t*>& encroached_edges,
            std::unordered_set<cell_t*>& bad_triangles, double rho_bar) {
        // Iterate over all currently encroached edges
        for (auto it = encroached_edges.begin(); it != encroached_edges.end(); ) {
            halfedge_t* e = *it;
            it = encroached_edges.erase(it);  // Remove the edge from the set to avoid reprocessing
            if (!e) continue;

            // Proceed only if none of the adjacent edges are "seditious" (i.e., dangerous to split now)
            if (!is_edge_seditious(e->next()) && !is_edge_seditious(e->prev()) &&
                !is_edge_seditious(e->next()->twin()) && !is_edge_seditious(e->prev()->twin())) {
                
                // Perform the segment split and update the DCEL
                split_subsegment(dcel, e, encroached_edges, bad_triangles, rho_bar);
                return true;  // Stop after the first successful split
            }
        }

        return false;  // No suitable edge was found for splitting
    }


    // Attempts to split the first bad triangle (with small angle or poor aspect ratio).
    // Returns true if a triangle was successfully split.
    static bool split_first_bad_triangle(dcel_t& dcel, double rho_bar, 
        std::unordered_set<halfedge_t*>& encroached_edges, std::unordered_set<cell_t*>& bad_triangles) {
        // Iterate over all currently bad triangles
        for (auto it = bad_triangles.begin(); it != bad_triangles.end(); ) {
            cell_t* t = *it;
            it = bad_triangles.erase(it);  // Remove the triangle from the set to avoid reprocessing

            // Skip if the triangle or its geometry is invalid
            if (!t || !t->halfedge()) continue;

            // Attempt to split the triangle by inserting its circumcenter
            bool split = split_triangle(dcel, t, encroached_edges, bad_triangles, rho_bar);
            if (split)
                return true;  // Stop after the first successful split
        }

        return false;  // No suitable triangle was split
    }


    // function that checks if edge defined by halfedge h is seditious
    // if it is, it is the triangle's shortest edge since its oppoing angle is < 60 degrees and the triangle is isosceles
    static bool is_edge_seditious(halfedge_t* h) {
        if (h->on_boundary()) return false;  // boundary edges can't be seditious
        if (!(h->prev()->on_boundary() && h->next()->on_boundary())) return false;  

        coords_t a = h->node()->coords();
        coords_t b = h->next()->node()->coords();
        coords_t c = h->prev()->node()->coords();
        // angle in c^ is angle of interest
        if (fdapde::internals::angle_between(b,c,a) >= 60) return false; // angle is not too small
        
        // 2 edges of h's cell need to have same length and to be midpoints of another segment 
        if (fdapde::internals::segment_length(c,a) != fdapde::internals::segment_length(c,b) )  return false;
        coords_t d = h->prev()->twin()->prev()->node()->coords();  
        coords_t e = h->next()->twin()->next()->next()->node()->coords();
        if ( !( fdapde::internals::collinear(c,a, d) && 
                fdapde::internals::collinear(c,b, e) ) )
                return false;
        return true;
    }   

    //same as flip function but keeping track of the encroached_edges and bad_triangles that is creating 
    static void flip_Ruppert(dcel_t& dcel, std::unordered_set<halfedge_t*>& encroached_edges,
        std::unordered_set<cell_t*>& bad_triangles, double rho_bar) {

        std::list<halfedge_t*> halfedges_to_check;
        for (auto it = dcel.halfedges_begin(); it != dcel.halfedges_end(); ++it,++it) {
            if(!it->on_boundary())
            halfedges_to_check.push_back(&(*it));
        }
        while (!halfedges_to_check.empty()) {
            halfedge_t* edge = halfedges_to_check.front();
            halfedges_to_check.pop_front();
            cell_t* neighbor = edge->twin()->cell();

            coords_t A = edge->node()->coords();
            coords_t B = edge->twin()->node()->coords();
            coords_t C = edge->prev()->node()->coords();
            coords_t D = edge->twin()->prev()->node()->coords();
    
            if (fdapde::internals::in_circle(A, B, C, D) || fdapde::internals::in_circle(A, D, B, C)) {
                
                halfedge_t* e = edge;
                bad_triangles.erase(e->cell());
                bad_triangles.erase(e->twin()->cell());
                dcel.remove_edge(edge);
                halfedge_t* new_edge = dcel.insert_edge(e->prev(), e->twin()->prev());
            
                if (new_edge) {
                    if(!new_edge->prev()->on_boundary())
                        halfedges_to_check.push_back(new_edge->prev());
                    if(!new_edge->next()->on_boundary())
                        halfedges_to_check.push_back(new_edge->next());
                    if(!new_edge->twin()->prev()->on_boundary())
                        halfedges_to_check.push_back(new_edge->twin()->prev());
                    if(!new_edge->twin()->next()->on_boundary())
                        halfedges_to_check.push_back(new_edge->twin()->next());

                    if(is_bad_triangle(dcel, new_edge->cell(), rho_bar)){
                        if (bad_triangles.count(new_edge->cell()) == 0)
                            bad_triangles.insert(new_edge->cell());
                    }
                    if(new_edge->prev()->on_boundary() && new_edge->prev()->cell() && check_encroachment(dcel, new_edge->prev())){
                        encroached_edges.insert(new_edge->prev());
                    }
                    if(new_edge->next()->on_boundary() && new_edge->next()->cell() && check_encroachment(dcel, new_edge->next())){
                        encroached_edges.insert(new_edge->next());
                    }
                    if(is_bad_triangle(dcel, new_edge->twin()->cell(), rho_bar)){
                        if (bad_triangles.count(new_edge->twin()->cell()) == 0)
                            bad_triangles.insert(new_edge->twin()->cell());
                    }
                   if(new_edge->twin()->prev()->on_boundary() && new_edge->twin()->prev()->cell() && check_encroachment(dcel, new_edge->twin()->prev())){
                        encroached_edges.insert(new_edge->twin()->prev());
                    }
                    if(new_edge->twin()->next()->on_boundary() && new_edge->twin()->next()->cell() && check_encroachment(dcel, new_edge->twin()->next())){
                        encroached_edges.insert(new_edge->twin()->next());
                    }
                }
            }
        }
    }
    //function that finds the triangle containing P scanning all the triangle of the trinagulation
    static cell_t* find_triangle(dcel_t& dcel, const coords_t& P) {
        
        for (auto it = dcel.cells_begin(); it != dcel.cells_end(); ++it) {  
            cell_t* cell = &(*it);  
            const coords_t& A = cell->halfedge()->node()->coords();
            const coords_t& B = cell->halfedge()->next()->node()->coords();
            const coords_t& C = cell->halfedge()->prev()->node()->coords();
    
            if (fdapde::internals::point_in_2d_tri(P, A, B, C)) {
                return cell;
            }
        }

    return nullptr;
    }
    
    //function working as mark_cavity but does not take into account conflits (since implement the Boyer-Watson algorithm)
    //but keeps track of the encroached edges and bad_triangles is creating 
    static void dig_cavity(dcel_t& dcel, node_t* u, halfedge_t* vw, std::unordered_set<halfedge_t*>& encroached_edges,
        std::unordered_set<cell_t*>& bad_triangles, double rho_bar) { 
        //if we are on the boundary we add the triangle  
        if(vw->on_boundary()){
            add_triangle(dcel, vw,std::vector<node_t*> {u});
            auto it_last = std::prev(dcel.cells_end());  
            cell_t* t = &(*it_last);
            if(is_bad_triangle(dcel, t, rho_bar)){
                if (bad_triangles.count(t) == 0)
                    bad_triangles.insert(t);
            }
            if (check_encroachment(dcel, vw)) {
                if (encroached_edges.count(vw) == 0)
                    encroached_edges.insert(vw);
            }
            return;
        }
        //finding the point adjacent to vw
        node_t* x = dcel.adjacent(vw);
        if (!x) {
            return;
        }
        bool ccw = fdapde::internals::are_2d_counterclockwise_sorted(u->coords(), vw->node()->coords(), vw->twin()->node()->coords());
        //test of circumcircle   
        bool inside;
        if (ccw) {
            inside = fdapde::internals::in_circle(u->coords(), vw->node()->coords(), vw->twin()->node()->coords(), x->coords());
        } else {
            inside = fdapde::internals::in_circle(u->coords(), vw->twin()->node()->coords(), vw->node()->coords(), x->coords());
        }
        if (inside) { 
        //falied the test so remove the triangle and expand the cavity on the remaining edges
            halfedge_t* wv = vw->twin();
            halfedge_t* vx = vw->twin()->next();
            halfedge_t* xw = vw->twin()->prev(); 
            bad_triangles.erase(vw->cell());
            bad_triangles.erase(vw->twin()->cell());
            dcel.remove_edge(vw);
            dig_cavity(dcel, u, vx, encroached_edges, bad_triangles, rho_bar);
            dig_cavity(dcel, u, xw, encroached_edges, bad_triangles, rho_bar);
        } else {
        //passed the test,adding the triangle 
            add_triangle(dcel, vw,std::vector<node_t*> {u});
            auto it_last = std::prev(dcel.cells_end());  
            cell_t* t = &(*it_last);
            if(is_bad_triangle(dcel, t, rho_bar)){
                if (bad_triangles.count(t) == 0)
                    bad_triangles.insert(t);
            }
            return;
        } 
    }
    
    //function working as insert_vertex_at_conflict but does not take into account conflits (since implement the Boyer-Watson algorithm)
    //but keeps track of the encroached edges and bad_triangles is creating 
    static void insert_vertex(dcel_t& dcel, node_t* u, cell_t* triangle, std::unordered_set<halfedge_t*>& encroached_edges,
        std::unordered_set<cell_t*>& bad_triangles, double rho_bar) {
        halfedge_t* vw = triangle->halfedge();
        halfedge_t* wx = vw->next();
        halfedge_t* xv = vw->prev();
    
        const coords_t& v = vw->node()->coords();
        const coords_t& w = wx->node()->coords();
        const coords_t& x = xv->node()->coords();
        bool found_on_edge = false;
    
        if (fdapde::internals::contains(u->coords(), v, w)) {
            halfedge_t* twin_next = vw->twin()->next();
            halfedge_t* twin_prev = vw->twin()->prev();
            bad_triangles.erase(vw->cell());
            bad_triangles.erase(vw->twin()->cell());
            dcel.remove_edge(vw);
            dig_cavity(dcel, u, wx, encroached_edges, bad_triangles, rho_bar);
            dig_cavity(dcel, u, xv, encroached_edges, bad_triangles, rho_bar);
            dig_cavity(dcel, u, twin_next, encroached_edges, bad_triangles, rho_bar);
            dig_cavity(dcel, u, twin_prev, encroached_edges, bad_triangles, rho_bar);
            found_on_edge = true;
        } else if (fdapde::internals::contains(u->coords(), w, x)) {
            halfedge_t* twin_next = wx->twin()->next();
            halfedge_t* twin_prev = wx->twin()->prev();
            bad_triangles.erase(wx->cell());
            bad_triangles.erase(wx->twin()->cell());
            dcel.remove_edge(wx);
            dig_cavity(dcel, u, vw, encroached_edges, bad_triangles, rho_bar);
            dig_cavity(dcel, u, xv, encroached_edges, bad_triangles, rho_bar);
            dig_cavity(dcel, u, twin_next, encroached_edges, bad_triangles, rho_bar);
            dig_cavity(dcel, u, twin_prev, encroached_edges, bad_triangles, rho_bar);
            found_on_edge = true;
        } else if (fdapde::internals::contains(u->coords(), x, v)) {
            halfedge_t* twin_next = xv->twin()->next();
            halfedge_t* twin_prev = xv->twin()->prev();
            bad_triangles.erase(xv->cell());
            bad_triangles.erase(xv->twin()->cell());
            dcel.remove_edge(xv);
            dig_cavity(dcel, u, vw, encroached_edges, bad_triangles, rho_bar);
            dig_cavity(dcel, u, wx, encroached_edges, bad_triangles, rho_bar);
            dig_cavity(dcel, u, twin_next, encroached_edges, bad_triangles, rho_bar);
            dig_cavity(dcel, u, twin_prev, encroached_edges, bad_triangles, rho_bar);
            found_on_edge = true;
        }
    
        if (!found_on_edge) {
            dig_cavity(dcel, u, vw, encroached_edges, bad_triangles, rho_bar);
            dig_cavity(dcel, u, wx, encroached_edges, bad_triangles, rho_bar);
            dig_cavity(dcel, u, xv, encroached_edges, bad_triangles, rho_bar);
        }
        //no need ti flip since Boyer-Watson mantains the Dleaunay
        //flip_Ruppert(dcel, encroached_edges, bad_triangles, rho_bar);
    }

    //helper function that fills a txt in order to then graph the trinagulation
    static void export_triangulation_to_txt(const triangulation_t& triangulation, const std::string& filename) {
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
    
};
  
}  // namespace fdapde

#endif // __DELAUNAY_H__