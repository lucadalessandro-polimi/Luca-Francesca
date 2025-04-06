#ifndef __FDAPDE_DELAUNAY_H__
#define __FDAPDE_DELAUNAY_H__

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
    
    //Constructors 
    // Costructor with random generated points
    Delaunay(const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& boundary,
        int N_internal = 100) :
    triangulation_t(build_triangulation_from_dcel(boundary, N_internal)) {}
    // Costructor with given internal points form the user
    Delaunay(const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& boundary,
            const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& internal) :
        triangulation_t(build_triangulation_from_dcel(boundary, internal)) {}


    void Ruppert_refinement(double rho_bar) {
        // Convert current triangulation to DCEL
        dcel_t dcel = Triangulation_to_DCEL(*this);
    
        while (true) {
            if (split_first_encroached_segment(dcel)) {
                continue;
            }
    
            if (split_first_bad_triangle(dcel, rho_bar)) {
                continue;
            }
    
            break;
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
        dcel.export_to_json("dcel_output.json");
    
        // Update the internal triangulation from the refined DCEL
        this->set_from_triangulation(DCEL_to_Triangulation(dcel));

    }
     

    
   private:
    static triangulation_t build_triangulation_from_dcel(
        const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& boundary, int N) 
    {
        dcel_t dcel = dcel_t::make_polygon(boundary);
        triangulate(dcel, N, boundary);  
        dcel.export_to_json("dcel_output.json");
        return DCEL_to_Triangulation(dcel);
    }


    static triangulation_t build_triangulation_from_dcel(
        const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& boundary,
        const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& internal)
    {
        dcel_t dcel = dcel_t::make_polygon(boundary);
        triangulate(dcel, internal, boundary);  
        dcel.export_to_json("dcel_output.json");
        return DCEL_to_Triangulation(dcel);
    }

    static triangulation_t DCEL_to_Triangulation(dcel_t& dcel) {
        Eigen::Matrix<double, Eigen::Dynamic, embed_dim> nodes(dcel.n_nodes(), embed_dim);
        Eigen::Matrix<int, Eigen::Dynamic, n_nodes_cell> cells(dcel.n_cells(), n_nodes_cell);
        Eigen::Matrix<int, Eigen::Dynamic, 1> boundary_markers(dcel.n_nodes());

        int node_idx = 0;
        for (auto it = dcel.nodes_begin(); it != dcel.nodes_end(); ++it) {
            nodes.row(node_idx) = it->coords().transpose();
            boundary_markers(node_idx) = it->on_boundary() ? 1 : 0;
            node_idx++;
        }

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
        std::string filename = "mesh_output.txt";
        export_triangulation_to_txt(triangulation, filename);

        std::string command = "python3 fdaPDE/src/plot_mesh.py";
        std::system(command.c_str()); 

        return triangulation;
    }


    static dcel_t Triangulation_to_DCEL(triangulation_t& triangulation) {
        using MatrixType = Eigen::Matrix<double, Eigen::Dynamic, embed_dim>;
        const auto& coords = triangulation.nodes();
        const auto& markers = triangulation.boundary_nodes();
    
        int n_boundary = triangulation.n_boundary_nodes();  
        MatrixType boundary_nodes(n_boundary, embed_dim);
    
        int idx = 0;
        for (int i = 0; i < coords.rows(); ++i) {
            if (markers(i,0) == 1) {
                boundary_nodes.row(idx++) = coords.row(i);
            }
        }
        dcel_t dcel = dcel_t::make_polygon(boundary_nodes);
        for (int i = 0; i < coords.rows(); ++i) {
            if (markers(i,0) == 0) {
                dcel.insert_node(node_t(dcel.n_nodes(), false, coords.row(i)));
            }
        }

        for (int i = 0; i < triangulation.n_cells(); ++i) {
            int id0 = triangulation.cells()(i, 0);
            int id1 = triangulation.cells()(i, 1);
            int id2 = triangulation.cells()(i, 2);
    
            node_t* n0 = std::addressof(*std::find_if(dcel.nodes_begin(), dcel.nodes_end(),
                                                      [&](const node_t& n) { return n.id() == id0; }));
            node_t* n1 = std::addressof(*std::find_if(dcel.nodes_begin(), dcel.nodes_end(),
                                                      [&](const node_t& n) { return n.id() == id1; }));
            node_t* n2 = std::addressof(*std::find_if(dcel.nodes_begin(), dcel.nodes_end(),
                                                      [&](const node_t& n) { return n.id() == id2; }));
            if (!fdapde::internals::are_2d_counterclockwise_sorted(n0->coords(), n1->coords(), n2->coords()))
            std::swap(n1, n2);

            halfedge_t* h = nullptr;

            h = dcel.find_halfedge_between(n0, n1);  
            if (h)
                add_triangle(dcel, h, std::vector<node_t*> {n2});
            else {
                h = dcel.find_halfedge_between(n1, n2);
                if (h)
                    add_triangle(dcel, h, std::vector<node_t*> {n0});
                else {
                    h = dcel.find_halfedge_between(n2, n0);
                    if (h)
                        add_triangle(dcel, h, std::vector<node_t*> {n1});
                    else{
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


    static void detect_conflicts(dcel_t& dcel ,node_t* n, const std::vector<halfedge_t*>& cells_to_check = {}) {
        bool found = false;

        if (cells_to_check.empty()) {  // initialization case: scan all the cells
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

    static void add_first_triangle(dcel_t& dcel, node_t* n, const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& boundary){   
        bool concave=false;
        auto iter = dcel.halfedges_begin();
        halfedge_t* first_h= dcel.emplace_halfedge_(n);
        auto& cell_begin = *dcel.cells_begin();
        first_h->set_cell(&cell_begin);

        // iterate over all boundary edges to connect to node, if possible
        for (int il = 0; il < boundary.rows(); ++il, ++iter) { 
            halfedge_t* v = &(*iter);                                                                  
            cell_t* c= v->cell();                                                                                    
            std::vector<halfedge_t*> halfedges_to_call(n_nodes_cell); 

            if(dcel.find_halfedge(n,c)){
                halfedge_t* h = dcel.find_halfedge(n,c);
                halfedges_to_call[2] = h;
            }
            else{  
                halfedges_to_call[2]= first_h;
            }
            halfedges_to_call[0] = v;
            if(v->next()->cell()==halfedges_to_call[2]->cell())
                halfedges_to_call[1] = v->next();
            else{
                halfedges_to_call[1] = dcel.find_halfedge(v->next()->node(),halfedges_to_call[2]->cell());
            }
            // add edges
            for (int i = 0; i < n_nodes_cell ; ++i) {
                halfedge_t* h1 = halfedges_to_call[i];
                halfedge_t* h2 = halfedges_to_call[(i + 1) % (n_nodes_cell)];
                
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
                    halfedge_t* h = dcel.insert_edge(h1, h2); 
                    if(h && h!=h1)
                        halfedges_to_call[(i + 1) % (n_nodes_cell)] = h->next();  
                } 
            }
        }
        
        if(concave){
            // cycle over boundary edges only for safe handling of edges to insert
            for(auto it=dcel.halfedges_begin(); it!=dcel.halfedges_end() && (&(*it))->cell(); ++it){
                halfedge_t* h = &(*it);
                int i=0;
                do{
                    h=h->next();
                    i++;
                }while(h!=&(*it));
                if(i>n_nodes_cell){ //not a triangle yet
                    for(int l=0; l<i-n_nodes_cell; ++l){
                      halfedge_t* n=h->next()->next();
                      if(h->node()->on_boundary() && h->next()->next()->node()->on_boundary() 
                         && !fdapde::internals::collinear(h->node()->coords(),h->next()->node()->coords(),h->next()->next()->node()->coords())
                         && fdapde::internals::are_2d_counterclockwise_sorted(h->node()->coords(),h->next()->node()->coords(),h->next()->next()->node()->coords()) ){
                        dcel.insert_edge(h, n);
                      }
                      h=n;   
                    }
                }
            }
        }
        flip(dcel);
    }

    static void flip(dcel_t& dcel) {

        // creating a list of halfedges to check whether they are locally delaunay or not (in this case flippable)
        std::list<halfedge_t*> halfedges_to_check;
        for (auto it = dcel.halfedges_begin(); it != dcel.halfedges_end(); ++it,++it) {
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
                dcel.remove_edge(edge);
                halfedge_t* new_edge = dcel.insert_edge(e->prev(), e->twin()->prev());
                std::cout << "FLIP" << std::endl;
            
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
    static void mark_cavity(dcel_t& dcel, node_t* u, halfedge_t* h, std::vector<halfedge_t*>& D, std::vector<halfedge_t*>& C) {
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


    // Function to insert a vertex handling conflicts
    static void insert_vertex_at_conflict(dcel_t& dcel, node_t* u) {
        // Retrieve the triangle in conflict with u and we marked as visited 
        cell_t* t = u->conflict(); 
        std::vector<halfedge_t*> D;
        std::vector<halfedge_t*> C;

        mark_cavity(dcel, u, t->halfedge(), D, C);
        mark_cavity(dcel, u, t->halfedge()->prev(), D, C);
        mark_cavity(dcel, u, t->halfedge()->next(), D, C);
  
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

    static void triangulate(dcel_t& dcel, int N, const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& boundary) {

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
                node_t* n = dcel.insert_node(node_t(dcel.n_nodes(), false, u));
                add_first_triangle(dcel, n, boundary);
                first_internal_id = dcel.n_nodes() - 1;
            } else {
                // Create the node and detect the conflicts with existing cells 
                node_t* n = dcel.insert_node(node_t(dcel.n_nodes(), false, u));
                detect_conflicts(dcel, n); 
            }
            ++generated_points;
        }
        // inserting the remaining nodes in the domain 
        for (auto it = dcel.nodes_begin(); it != dcel.nodes_end(); ++it) {
            node_t* u = &(*it);
            if (!u->on_boundary() && u->id()!=first_internal_id)   
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

    //overloaded one if user wants to pass manually the internal points
    //the user must know the passed internal points lie all inside the domain 
    static void triangulate(dcel_t& dcel, const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& internal,
        const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& boundary) {
        
        node_t* n = dcel.insert_node(node_t(dcel.n_nodes(), false, internal.row(0)));
        add_first_triangle(dcel, n, boundary);
        int first_internal_id = dcel.n_nodes() - 1;

        //inserting the remaining internal points in the triangulation
        for (int i = 1; i < internal.rows(); ++i) {
            node_t* n = dcel.insert_node(node_t(dcel.n_nodes(), false, internal.row(i).transpose().eval()));
            //detecting the conflicts with existing cells
            detect_conflicts(dcel, n);
        }
        //inserting the remaining nodes in the domain
        for (auto it = dcel.nodes_begin(); it != dcel.nodes_end(); ++it) {
            node_t* u = &(*it);
            if (!u->on_boundary() && u->id()!=first_internal_id)  
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

    static void split_subsegment(dcel_t& dcel, halfedge_t* e) {
        node_t* a = e->node();
        node_t* b = e->twin()->node();
        coords_t mid = 0.5 * (a->coords() + b->coords());
        // inserting the midpoint
        node_t* m = dcel.insert_node(node_t(dcel.n_nodes(), true, mid));
        // cutting in two the cell of the encroached segment
        halfedge_t* prev = e->prev();
        add_triangle(dcel, e, std::vector<node_t*> {m});
        dcel.remove_edge(e);
        add_triangle(dcel, prev, std::vector<node_t*> {m});
    }

    void split_triangle(dcel_t& dcel, cell_t* t) {
        coords_t A = t->halfedge()->prev()->node()->coords();
        coords_t B = t->halfedge()->node()->coords();
        coords_t C = t->halfedge()->next()->node()->coords();
        coords_t c = fdapde::internals::circumcenter(A, B, C);

        //find if c encroaches some edge of the trinagulation
        for (auto it = dcel.halfedges_begin(); it != dcel.halfedges_end(); ++it) {
            halfedge_t* e = &(*it);
            //in order to evaluate only the segments of the PLC (i.e. the boundary of the domain)
            if (!e->on_boundary()) continue;
            // o si fa cosi o si tiene in memoria nella delaunay.h un contatore di nodi al bordo per gestire meglio il bordo
            if (e->on_boundary() && e->cell()==nullptr) continue;
            
            coords_t a = e->node()->coords();
            coords_t b = e->twin()->node()->coords();

        if (fdapde::internals::is_encroached(c, a, b)) {
                split_subsegment(dcel, e); // if encroaches, then split the subsegment
                return;
            }
        }
        // otherwise c is inserted as node
        //grazie a questa frase: Sì, può cadere sul bordo, e viene accettato se non encroacha un subsegmento 
        //(cioè se sta esattamente sul bordo ma non dentro al disco di diametro di un segmento).
        //dovrebbe essere la prova del fatto che c non puo mai cadere sul bordo del dominio dove ci sono tutti i segmenti della PLC
        node_t* circ = dcel.insert_node(node_t(dcel.n_nodes(), false, c));
        insert_vertex(dcel, circ, find_triangle(dcel,c));
    }
   
    static bool split_first_encroached_segment(dcel_t& dcel) {
        halfedge_t* encroached_edge = nullptr;

        for (auto it = dcel.halfedges_begin(); it != dcel.halfedges_end(); ++it) {
            halfedge_t* e = &(*it);

            if (!e->on_boundary()) continue;
            if (e->on_boundary() && e->cell()==nullptr) continue;

            coords_t a = e->node()->coords();
            coords_t b = e->twin()->node()->coords();

            for (auto nit = dcel.nodes_begin(); nit != dcel.nodes_end(); ++nit) {
                node_t* n = &(*nit);
                if (n == e->node() || n == e->twin()->node()) continue;

                if (fdapde::internals::is_encroached(n->coords(), a, b)) {
                    encroached_edge = e;
                    break;
                }
            }

            if (encroached_edge != nullptr)
                break;
        }

        if (encroached_edge != nullptr) {
            split_subsegment(dcel,encroached_edge);
            return true;
        }

        return false;
    }

    bool split_first_bad_triangle(dcel_t& dcel, double rho_bar) {
        for (auto it = dcel.cells_begin(); it != dcel.cells_end(); ++it) {
            cell_t* t = &(*it);

            coords_t A = t->halfedge()->prev()->node()->coords();
            coords_t B = t->halfedge()->node()->coords();
            coords_t C = t->halfedge()->next()->node()->coords();

            double ratio = fdapde::internals::radius_edge_ratio(A, B, C);

            if (ratio > rho_bar) {
                std::cout << "Splitting triangle with ratio = " << ratio << " > " << rho_bar << "\n";
                split_triangle(dcel, t);
                return true;
            }
        }
        return false;
    }

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


    static cell_t* find_triangle(dcel_t& dcel, const coords_t& P) {
        
        for (auto it = dcel.cells_begin(); it != dcel.cells_end(); ++it) {
            cell_t* cell = &(*it);  
    
            const coords_t& A = cell->halfedge()->node()->coords();
            const coords_t& B = cell->halfedge()->next()->node()->coords();
            const coords_t& C = cell->halfedge()->prev()->node()->coords();

        /*  CAMBIARE IL MODO IN CUI SI VERIFICA SE CASCA O MENO SUL LATO
            if (fdapde::internals::contains(P, A, B) || 
            fdapde::internals::contains(P, B, C) || 
            fdapde::internals::contains(P, C, A)) {
            //I remove from internal_points vector
                auto it = std::remove_if(internal_points_.begin(), internal_points_.end(), 
                [&](const coords_t& point) { return point.isApprox(P); });
                internal_points_.erase(it, internal_points_.end());

                return nullptr;
            }*/
    
            if (fdapde::internals::point_in_2d_tri(P, A, B, C)) {
                return cell;
            }
        }

    return nullptr;
    }

    static void dig_cavity(dcel_t& dcel, node_t* u, halfedge_t* vw) { 
    //if we are on the boundary we add the triangle  
        if(vw->on_boundary()){
            add_triangle(dcel, vw,std::vector<node_t*> {u});
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
            
            dcel.remove_edge(vw);
            dig_cavity(dcel, u, vx);
            dig_cavity(dcel, u, xw);
        } else {
        //passed the test,adding the triangle 
            add_triangle(dcel, vw,std::vector<node_t*> {u});
            return;
        } 
    }
    
    
    static void insert_vertex(dcel_t& dcel, node_t* u, cell_t* triangle) {
        
        halfedge_t* vw = triangle->halfedge();
        halfedge_t* wx = vw->next();
        halfedge_t* xv = vw->prev();
        // expanding cavity
        dig_cavity(dcel, u, vw);
        dig_cavity(dcel, u, wx);
        dig_cavity(dcel, u, xv);
    }

    void set_from_triangulation(const triangulation_t& triang) {
        this->nodes_ = triang.nodes();
        this->cells_ = triang.cells();
        this->boundary_markers_ = triang.boundary_nodes();
        this->neighbors_ = triang.neighbors();
        this->n_nodes_ = triang.n_nodes();
        this->n_cells_ = triang.n_cells();
    }
    
};
  
}  // namespace fdapde

#endif // __DELAUNAY_H__