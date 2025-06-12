#ifndef _FDAPDE_DELAUNAY_H_
#define _FDAPDE_DELAUNAY_H_

//Delaunay class build a triangulation and refine (using Ruppert algorithm) a generic connected domain 
//including concavities 
#include "header_check.h"
namespace fdapde {
  
template <int LocalDim, int EmbedDim>
class Delaunay {
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

    
    /*  DECORATOR
    struct node_t_;
    struct cell_t_;

    struct node_t_ : public DCEL<local_dim, embed_dim>::node_t {

        using Base = typename DCEL<local_dim, embed_dim>::node_t;

        node_t_() noexcept = default;
        template <typename CoordsType>
             requires(internals::is_eigen_dense_xpr_v<CoordsType>)
        node_t_(int id, halfedge_t* halfedge, bool boundary, const CoordsType& coords) : Base(id, halfedge, boundary, coords), conflicting_triangle_(nullptr) {}
        template <typename CoordsType>
             requires(internals::is_eigen_dense_xpr_v<CoordsType>)
        node_t_(int id, bool boundary, const CoordsType& coords) : Base(id, boundary, coords), conflicting_triangle_(nullptr) {}
        template <typename... CoordsType>
            requires(std::is_floating_point_v<CoordsType> && ...) && (sizeof...(CoordsType) == embed_dim)
        node_t_(int id, halfedge_t* halfedge, bool boundary, CoordsType&&... coords) : Base(id, halfedge, boundary, coords...), conflicting_triangle_(nullptr) {}
        template <typename... CoordsType>
            requires(std::is_floating_point_v<CoordsType> && ...) && (sizeof...(CoordsType) == embed_dim)
        node_t_(int id, bool boundary, CoordsType&&... coords) :  Base(id, nullptr, boundary, coords...), conflicting_triangle_(nullptr) { }

        void set_conflict(cell_t_* triangle) { conflicting_triangle_ = triangle; }
        cell_t_* conflict() const { return conflicting_triangle_; }
        void remove_conflict() { conflicting_triangle_ = nullptr; }

       private:
        cell_t_* conflicting_triangle_; 
    };

    struct cell_t_ : public DCEL<local_dim, embed_dim>::cell_t {

        using Base = typename DCEL<local_dim, embed_dim>::cell_t;
        
        cell_t_() : Base() { }
        cell_t_(int id) : Base(id) { }
        cell_t_(int id, halfedge_t* h) : Base(id, h) { }
        cell_t_(const cell_t& base) : Base(base.id(), base.halfedge()) { }
        
        void add_conflict(node_t_* point) { conflicting_points_.push_back(point); }
        const std::vector<node_t_*>& conflicting_points() const{ return conflicting_points_; }
        std::vector<node_t_*>& conflicting_points() { return conflicting_points_; }
        void clear_conflicts() { conflicting_points_.clear(); }

       private:
        std::vector<node_t_*> conflicting_points_;
    };
    */

    
    
    // constructors 
    // costructor with random generated points
    Delaunay(const std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>>& boundaries,  int N=100, const std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>>>& holes = {{}}){
        triangulate(N, boundaries, holes);
    }    
    // costructor with given internal points form the user
    // user needs to provide internal points correctly located inside the domain 
    Delaunay(const std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>>& boundaries, const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& internal, const std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>>>& holes = {{}}) {
        triangulate(internal, boundaries, holes);
    }

    // Getter const
    const dcel_t& dcel() const {
        return dcel_;
    }
    // Getter non-const
    dcel_t& dcel() {
        return dcel_;
    }

    triangulation_t triangulation() {
        return dcel_.template to_triangulation<triangulation_t>();
    }

    //function running the refinment with Ruppert algorithm 
    //for the moment it can work only from the inside with default parameters in the constructor
    void Ruppert_refinement(double rho_bar) {
        //set needed to save in memory edges that encroach a point in the triangulation and another one for the badly shaped traingles 
        std::unordered_set<halfedge_t*> encroached_edges;
        std::unordered_set<cell_t*> bad_triangles;     
        std::unordered_set<halfedge_t*> boundary_edges;

        // inizialization of the two set
        for (auto it = dcel_.halfedges_begin(); it != dcel_.halfedges_end(); ++it) {
            halfedge_t* e = &(*it);
            //since in the previuos code we do not touch the boundary we can stop at the twin of the first boundary edge
            if(e->on_boundary() && !e->cell()) continue;  //otherwise it doesn't do the holes
            //if (e->on_boundary()){
            if(e->is_subsegment()) {
                boundary_edges.insert(e);
                if(check_encroachment(e)){
                    //to make sure the set saves only one copy of the edge
                    if (encroached_edges.count(e) == 0)
                        encroached_edges.insert(e);
                }
            }
        }
        for (auto it = dcel_.cells_begin(); it != dcel_.cells_end(); ++it) {
            cell_t* t = &(*it);
            if (is_bad_triangle(t, rho_bar)) {
                if (bad_triangles.count(t) == 0)
                    bad_triangles.insert(t);
            }
        }  
        //while keeps running until the two set are empty and every time a bad triangle is exiting the triangulation
        //the test of the encroached edges runs in order to keep track of the newly created triangulation
        int cont1 = 0;
        while (true) {
            if (split_first_encroached_segment(boundary_edges, encroached_edges, bad_triangles, rho_bar)) {
                cont1++;
                //if(cont1==3) break;
                continue;
            }
            if (split_first_bad_triangle(rho_bar, boundary_edges, encroached_edges, bad_triangles)) {
                continue;
            }
            break;
        }
        // final reorder to cut no longer existing id of cells and edges
        int cont = 0;
        for (auto it = dcel_.cells_begin(); it != dcel_.cells_end(); ++it)
            it->set_id(cont++);
        dcel_.set_n_cells_(cont);
    
        cont = 0;
        for (auto it = dcel_.halfedges_begin(); it != dcel_.halfedges_end(); ++it)
            it->set_id(cont++);
        
        //json needed for debug
        dcel_.export_to_json("dcel_output.json");
    }


   private:
    dcel_t dcel_;


    void triangulate(int N, const std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>>& boundaries_entry, const std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>>>& holes_entry ) {
        
        // compute bounding box of the input polygon
        double min_x = boundaries_entry[0].col(0).minCoeff();
        double max_x = boundaries_entry[0].col(0).maxCoeff();
        double min_y = boundaries_entry[0].col(1).minCoeff();
        double max_y = boundaries_entry[0].col(1).maxCoeff();
        int generated_points = 0;

        // Estimate grid resolution based on desired number of interior points
        int points_per_row = static_cast<int>(std::ceil(std::sqrt(N)));
        double dx = (max_x - min_x) / (points_per_row + 1);
        double dy = (max_y - min_y) / (points_per_row + 1);

        // Use average dimension to scale the perturbation uniformly in all directions
        double width = max_x - min_x;
        double height = max_y - min_y;
        double avg_dim = (width + height) / 2.0;

        // Random perturbation to avoid aligned grid artifacts
        double perturbation_scale = avg_dim / (points_per_row + 1);
        std::mt19937 gen(42);  // Deterministic seed for reproducibility
        std::uniform_real_distribution<double> pert_x(-perturbation_scale / 2, perturbation_scale / 2);
        std::uniform_real_distribution<double> pert_y(-perturbation_scale / 2, perturbation_scale / 2);

        // check if boundaries are counterclockwise sorted
        std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>> boundaries;
        for(const auto& bd: boundaries_entry){
            if (!internals::are_2d_counterclockwise_sorted(bd)) {
                Eigen::Matrix<double, Eigen::Dynamic, embed_dim> reversed_bd(bd.rows(), embed_dim);
                for (int i = 0; i < bd.rows(); ++i)
                    reversed_bd.row(i) = bd.row(bd.rows() - 1 - i);
                boundaries.push_back(reversed_bd);
            }
            else{
                boundaries.push_back(bd);
            }
        }
        std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>>> holes;
        for(const auto& hole_vect : holes_entry){
            std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>> new_vector;
            for (const auto& hole : hole_vect) {
                if (internals::are_2d_counterclockwise_sorted(hole)) {
                    Eigen::Matrix<double, Eigen::Dynamic, embed_dim> reversed_hole(hole.rows(), embed_dim);
                    for (int i = 0; i < hole.rows(); ++i)
                        reversed_hole.row(i) = hole.row(hole.rows() - 1 - i);
                    new_vector.push_back(reversed_hole);
                } else {
                    new_vector.push_back(hole);
                }
            }
            holes.push_back(new_vector);
        }
        
        std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>> boundary_vertices=boundaries; 
        for(int i=1; i< boundaries.size(); ++i){   
            boundary_vertices[i] = split_boundary_points(boundaries[i]);
        }
        std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>>> holes_vertices(holes.size());
        for(int j=0; j< holes.size(); ++j){
            std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>> holes_vertices_i(holes[j].size());
            for(int i=0; i< holes[j].size(); ++i){
                holes_vertices_i[i] = split_boundary_points(holes[j][i]);
            }
            holes_vertices[j] = holes_vertices_i;
        }
        if(boundaries.size()==1)
            boundary_vertices[0] = split_boundary_points(boundaries[0]);
        else
            boundary_vertices[0] = split_boundary_points_with_attachments(boundary_vertices);
        
        // Initialize the triangulation with the boundary polygon
        initialize_triangulation(boundary_vertices, holes_vertices);
        
        if(boundaries.size()==1) complete_boundary(boundaries[0], boundary_vertices[0]);
        // add boundary vertices to the DCEL
        for(int i=1; i< boundaries.size(); ++i){
            complete_boundary(boundaries[i], boundary_vertices[i]);
        }
        for(int j=1; j< holes.size(); ++j){
            for(int i=0; i< holes[j].size(); ++i){
                complete_boundary(holes[j][i], holes_vertices[j][i]);
            }
        }
        
        flip();  // Ensure boundary triangulation satisfies Delaunay property
        
        int n_nodes_boundaries = dcel_.n_nodes();
        // Generate and insert interior points into the DCEL
        for (int i = 1; i <= points_per_row; ++i) {
            for (int j = 1; j <= points_per_row; ++j) {
                coords_t u;
                u << min_x + i * dx + pert_x(gen),
                    min_y + j * dy + pert_y(gen);
                // Keep only points that lie inside the polygonal domain
                if (!fdapde::internals::point_safely_in_polygon(boundaries[0], holes, u, perturbation_scale/16)){
                    //std::cout<<"SCARTO PUNTO CON COORDINATE: "<<u<<std::endl;
                    continue;
                }

                node_t* n = dcel_.insert_node(node_t(dcel_.n_nodes(), false, u));
                detect_conflicts(n);
                ++generated_points;
            }
        }
        
        int cont_pt=0;
        // Insert all internal points into the triangulation using the conflict graph
        for (auto it = dcel_.nodes_begin(); it != dcel_.nodes_end(); ++it, cont_pt++) {
            if(cont_pt < n_nodes_boundaries) continue; // Skip boundary nodes
            node_t* u = &(*it);
            insert_vertex_at_conflict(u);
        }
        
        // Reassign consecutive IDs to all cells and half-edges for consistency
        int cont = 0;
        for (auto it = dcel_.cells_begin(); it != dcel_.cells_end(); ++it)
            it->set_id(cont++);
        dcel_.set_n_cells_(cont);

        cont = 0;
        for (auto it = dcel_.halfedges_begin(); it != dcel_.halfedges_end(); ++it)
            it->set_id(cont++);
        
        // Export the resulting DCEL structure to a JSON file for visualization
        dcel_.export_to_json("dcel_output.json");
    }
    
    //overloaded one if user wants to pass manually the internal points
    //the user must know the passed internal points lie all inside the domain 
    void triangulate(const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& internal,
        const std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>>& boundaries_entry, const std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>>>& holes_entry) {
        
        // check if boundaries are counterclockwise sorted
        std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>> boundaries;
        for(const auto& bd: boundaries_entry){
            if (!internals::are_2d_counterclockwise_sorted(bd)) {
                Eigen::Matrix<double, Eigen::Dynamic, embed_dim> reversed_bd(bd.rows(), embed_dim);
                for (int i = 0; i < bd.rows(); ++i)
                    reversed_bd.row(i) = bd.row(bd.rows() - 1 - i);
                boundaries.push_back(reversed_bd);
            }
            else{
                boundaries.push_back(bd);
            }
        }
        std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>>> holes;
        for(const auto& hole_vect : holes_entry){
            std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>> new_vector;
            for (const auto& hole : hole_vect) {
                if (internals::are_2d_counterclockwise_sorted(hole)) {
                    Eigen::Matrix<double, Eigen::Dynamic, embed_dim> reversed_hole(hole.rows(), embed_dim);
                    for (int i = 0; i < hole.rows(); ++i)
                        reversed_hole.row(i) = hole.row(hole.rows() - 1 - i);
                    new_vector.push_back(reversed_hole);
                } else {
                    new_vector.push_back(hole);
                }
            }
            holes.push_back(new_vector);
        }

        std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>> boundary_vertices=boundaries; 
        for(int i=1; i< boundaries.size(); ++i){   
            boundary_vertices[i] = split_boundary_points(boundaries[i]);
        }
        std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>>> holes_vertices(holes.size());
        for(int j=0; j< holes.size(); ++j){
            std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>> holes_vertices_i(holes[j].size());
            for(int i=0; i< holes[j].size(); ++i){
                holes_vertices_i[i] = split_boundary_points(holes[j][i]);
            }
            holes_vertices[j] = holes_vertices_i;
        }
        if(boundaries.size()==1)
            boundary_vertices[0] = split_boundary_points(boundaries[0]);
        else
            boundary_vertices[0] = split_boundary_points_with_attachments(boundary_vertices);
         
        // Initialize the triangulation with the boundary polygon
        initialize_triangulation(boundary_vertices, holes_vertices);
        
        if(boundaries.size()==1) complete_boundary(boundaries[0], boundary_vertices[0]);
        // add boundary vertices to the DCEL
        for(int i=1; i< boundaries.size(); ++i){
            complete_boundary(boundaries[i], boundary_vertices[i]);
        }
        for(int j=1; j< holes.size(); ++j){
            for(int i=0; i< holes[j].size(); ++i){
                complete_boundary(holes[j][i], holes_vertices[j][i]);
            }
        }
        
        flip();

        int n_nodes_boundaries = dcel_.n_nodes();

        //inserting the internal points in the triangulation
        for (int i = 0; i < internal.rows(); ++i) {
            node_t* n = dcel_.insert_node(node_t(dcel_.n_nodes(), false, internal.row(i).transpose().eval()));
            detect_conflicts(n);
        }

        int cont_pt=0;
        for (auto it = dcel_.nodes_begin(); it != dcel_.nodes_end(); ++it, cont_pt++) {
            if(cont_pt < n_nodes_boundaries ) continue;
            node_t* u = &(*it);
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
        dcel_.export_to_json("dcel_output.json");
    }

    //function performing the first raw triangulation of the domain using the polygon.h class
    void initialize_triangulation(const std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>>& boundaries, std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>>>& holes) {
        if(holes.size()!=0)
            dcel_= dcel_t::make_polygon(boundaries[0], holes[0]);
        else{
            dcel_= dcel_t::make_polygon(boundaries[0]);
            holes.push_back({});
        }
       

        if(boundaries.size() == 1 ){
            polygon_t polygon(boundaries[0], holes[0]);
            auto triangulation = polygon.triangulation();
            const auto& nodes = triangulation.nodes(); 
            dcel_.from_triangulation(triangulation, holes[0]);
        }
        else for(int i=1; i< boundaries.size(); ++i){
            cell_t* c= &(*std::prev(dcel_.cells_end()));   //è GIUSTO ?????
            for(int j=0; j< boundaries[i].rows(); ++j){
                    coords_t co= boundaries[i].row(j).transpose();
                    if(dcel_.find_node(co)) continue;
                    node_t* n1 = dcel_.insert_node(node_t(dcel_.n_nodes(), false, co));
                    halfedge_t* h1 = dcel_.emplace_halfedge_(n1, true);
                    n1->set_halfedge(h1);
                    h1->set_cell(c);
            } 
            cell_t* c_holes= nullptr;
            for (int j = 0; j < boundaries[i].rows(); ++j) {
                    coords_t co1= boundaries[i].row(j).transpose();
                    coords_t co2= boundaries[i].row( (j+1) % boundaries[i].rows() ).transpose();
                    node_t* n1= dcel_.find_node(co1);
                    node_t* n2= dcel_.find_node(co2);
                    std::cout << "n1: " << n1->id() << " n2: " << n2->id() << std::endl;
                    if(!dcel_.find_halfedge_between(co1,co2)){
                        /*halfedge_t* h1 = n1->halfedge();
                        halfedge_t* h2 = n2->halfedge();
                        halfedge_t* h_new=dcel_.insert_edge(h1, h2);
                        h_new->set_subsegment(true);
                        h_new->twin()->set_subsegment(true);
                        c_holes = h_new->cell();  */
                        insert_collinear_chain(n1, n2);
                        //c_holes= n1->halfedge()->cell();
                    }
                    if(c_holes == nullptr) {
                        c_holes = dcel_.find_halfedge_between(co1, co2)->cell();
                    }
            } 
            if(holes.size() == i)
                holes.push_back({});
            polygon_t polygon(boundaries[i], holes[i]);
            auto triangulation = polygon.triangulation();
            const auto& nodes = triangulation.nodes(); 
            // update holes[i] edges to the cell of boundary i
            for (int j = 0; j < holes[i].size(); ++j) {
                coords_t co = holes[i][j].row(0).transpose();
                if (!dcel_.find_node(co)) continue;
                node_t* n1 = dcel_.find_node(co);
                halfedge_t* h_hole= n1->halfedge();   //make_polygon creates holes' nodes s.t. its own halfedge is defined 
                h_hole->set_cell(c_holes);
                halfedge_t* next= h_hole->next();
                do{
                    next->set_cell(c_holes);
                    next = next->next();
                }while(next!= h_hole);
            }
            dcel_.from_triangulation(triangulation, holes[i]);
            fix_edges_over_collinear_nodes();
        
        }
        // check on the DCEL to connect any nodes that might appear in only some regions and not in others (vertices in a region and collinear in another)
        for (auto it = dcel_.cells_begin(); it != dcel_.cells_end(); ++it) {
            cell_t* c = &(*it);
            halfedge_t* h = c->halfedge();
            halfedge_t* h_nn = h->next()->next();
            while(true){
                if(!fdapde::internals::collinear(h->node()->coords(), h->next()->node()->coords(), h_nn->node()->coords())){
                    dcel_.insert_edge(h, h_nn);
                    break;
                }
                h=h->next();
                h_nn= h_nn->next();
            }
        }
        
    }

    void fix_edges_over_collinear_nodes() {
        std::unordered_set<halfedge_t*> edges_to_replace;

        // find all edges crossing existing collinear edges
        for (auto it = dcel_.cells_begin(); it != dcel_.cells_end(); ++it) {
            cell_t* c = &(*it);
            halfedge_t* h_start = c->halfedge();
            if (!h_start) continue;

            halfedge_t* h = h_start;
            do {
                node_t* n1 = h->node();
                node_t* n2 = h->twin()->node();

                if (!n1 || !n2) continue;

                coords_t A = n1->coords();
                coords_t B = n2->coords();
                coords_t AB = B - A;
                double ab2 = AB.squaredNorm();
                if (ab2 < 1e-12) continue;

                for (auto nit = dcel_.nodes_begin(); nit != dcel_.nodes_end(); ++nit) {
                    node_t* P = &(*nit);
                    if (P == n1 || P == n2) continue;

                    coords_t p = P->coords();
                    if (!fdapde::internals::collinear(A, p, B)) continue;

                    coords_t AP = p - A;
                    double t = AB.dot(AP) / ab2;
                    if (t > 1e-6 && t < 1.0 - 1e-6) {
                        if (!edges_to_replace.count(h) && !edges_to_replace.count(h->twin())) {
                            edges_to_replace.insert(h);
                        }
                        break;
                    }
                }

                h = h->next();
            } while (h && h != h_start);
        }

        // substitute problematic edges with collinear chains
        for (halfedge_t* h : edges_to_replace) {
            node_t* A = h->node();
            node_t* B = h->twin()->node();
            if (!A || !B) continue;
            std::cout << "h removed: " << h->id() << std::endl;
            dcel_.remove_edge(h);
        }
    }


    void insert_collinear_chain(node_t* A, node_t* B) {
        coords_t pA = A->coords();
        coords_t pB = B->coords();
        coords_t AB = pB - pA;
        double ab_norm2 = AB.squaredNorm();

        std::vector<std::pair<double, node_t*>> intermediate;

        for (auto it = dcel_.nodes_begin(); it != dcel_.nodes_end(); ++it) {
            node_t* P = &(*it);
            if (P == A || P == B) continue;
            coords_t p = P->coords();
            if (!fdapde::internals::collinear(pA, p, pB)) continue;
            coords_t AP = p - pA;
            double t = AB.dot(AP) / ab_norm2;
            if (t > 1e-6 && t < 1.0 - 1e-6) {
                intermediate.emplace_back(t, P);
            }
        }

        std::sort(intermediate.begin(), intermediate.end(),
                [](const std::pair<double, node_t*>& a, const std::pair<double, node_t*>& b) {
                    return a.first < b.first;
                });

        // Construct chain A → P1 → ... → Pk → B
        node_t* prev = A;
        for (auto& [_, curr] : intermediate) {
            halfedge_t* h1 = prev->halfedge();
            halfedge_t* h2 = curr->halfedge();
            if (h1 && h2) {
                halfedge_t* h_new = dcel_.insert_edge(h1, h2);
                h_new->set_subsegment(true);
                h_new->twin()->set_subsegment(true);
            }
            prev = curr;
        }
        // Connect last point to B
        halfedge_t* h1 = prev->halfedge();
        halfedge_t* h2 = B->halfedge();
        if (h1 && h2) {
            halfedge_t* h_new = dcel_.insert_edge(h1, h2);
            h_new->set_subsegment(true);
            h_new->twin()->set_subsegment(true);
        }
    }


    // function to divide the boundary points into vertex points and collinear points 
    Eigen::Matrix<double, Eigen::Dynamic, 2>  split_boundary_points(const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& boundary) {
        std::list<Eigen::Matrix<double, 1, 2>> vertex_points_list;

        int n = boundary.rows();

        // iterate through the boundary points
        for (int i = 0; i < n; ++i) {
            const auto& point = boundary.row(i);
            const auto& point_prev = boundary.row((i-1+n)%n);
            const auto& point_next = boundary.row((i+1)%n);
            if (!fdapde::internals::collinear(point_prev, point, point_next)) 
                vertex_points_list.push_back(point);
        }

        Eigen::Matrix<double, Eigen::Dynamic, 2> vertex_points(vertex_points_list.size(), 2);
        int i=0;
        for (const auto & point : vertex_points_list) {
            vertex_points.row(i++) = point;
        }
        
        return vertex_points;
    }

    Eigen::Matrix<double, Eigen::Dynamic, 2> split_boundary_points_with_attachments(const std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>>& boundary_vertices) 
    {
        std::list<Eigen::Matrix<double, 1, 2>> vertex_points_list;

        // Definizione locale del comparatore per righe
        auto row_matrix_less = [](const Eigen::Matrix<double, 1, 2>& a,const Eigen::Matrix<double, 1, 2>& b) -> bool 
        {
            if (std::abs(a(0) - b(0)) > 1e-10) return a(0) < b(0);
            return a(1) < b(1) - 1e-10;
        };

        std::set<Eigen::Matrix<double, 1, 2>, decltype(row_matrix_less)> internal_points_set(row_matrix_less);

        // Costruisce il set dei punti interni
        for (int j=1; j< boundary_vertices.size(); ++j) 
            for (int i = 0; i < boundary_vertices[j].rows(); ++i) 
                internal_points_set.insert(boundary_vertices[j].row(i));

        int n = boundary_vertices[0].rows();
        // Scorre i punti del bordo esterno
        for (int i = 0; i < n; ++i) {
            const auto& point = boundary_vertices[0].row(i);
            const auto& point_prev = boundary_vertices[0].row((i - 1 + n) % n);
            const auto& point_next = boundary_vertices[0].row((i + 1) % n);
            bool is_vertex = !fdapde::internals::collinear(point_prev, point, point_next);
            bool is_attachment = internal_points_set.find(point) != internal_points_set.end();
            if (is_vertex || is_attachment)
                vertex_points_list.push_back(point);
        }

        Eigen::Matrix<double, Eigen::Dynamic, 2> vertex_points(vertex_points_list.size(), 2);
        int i = 0;
        for (const auto& point : vertex_points_list) {
            vertex_points.row(i++) = point;
        }
        return vertex_points;
    }


    void complete_boundary(const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& boundary, const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& boundary_vertices ) {
        if(boundary.rows()==boundary_vertices.rows()) return; //if the boundary is already complete we do not need to do anything

        std::cout << "in COMPLETE!!" << std::endl;
        auto row_in_matrix = [](const Eigen::Matrix<double, 1, embed_dim>& row, const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& mat) -> bool {
            for (int i = 0; i < mat.rows(); ++i) {
                if (mat.row(i).isApprox(row))
                    return true;
            }
            return false;
        };
        
        /*int j=0;
        int k=0;
        for(int i=0; i < boundary.rows(); i=j){
            if(!row_in_matrix(boundary.row(i), boundary_vertices)){
                j=i;
                coords_t n1= boundary_vertices.row((k-1+ boundary_vertices.rows())%boundary_vertices.rows());
                coords_t n2= boundary_vertices.row(k%boundary_vertices.rows());
                std::cout << "in for " << boundary.row(j).transpose() << std::endl;
                if(dcel_.find_node(boundary.row(j))){
                    j++;
                    k++;
                    continue;
                }

                halfedge_t* e= dcel_.find_halfedge_between(n1, n2);
                node_t* m;
                std::cout << "n1: " << n1 << " n2: " << n2 << std::endl;
                if(!e)
                    e= dcel_.find_halfedge_between(n1, boundary.row(0));
                if(e->on_boundary()){
                    m = dcel_.insert_node(node_t(dcel_.n_nodes(), true, boundary.row(j)));
                }
                else
                    m = dcel_.insert_node(node_t(dcel_.n_nodes(), false, boundary.row(j)));
                halfedge_t* prev = e->prev();
                halfedge_t* next = e->next();
                halfedge_t* twin_prev = e->twin()->prev();
                halfedge_t* h1 = dcel_.emplace_halfedge_(m, true);
                h1->set_cell(next->cell());
                dcel_.insert_edge(next, h1);
                h1->twin()->set_subsegment(true);
                halfedge_t* h2 = dcel_.insert_edge(prev->next(), h1);
                h2->set_subsegment(true);
                h2->twin()->set_subsegment(true);
                dcel_.insert_edge(prev,h1);
                // Remove the encroached edge from the DCEL
                dcel_.remove_edge(e);
                
                if(!h1->on_boundary()) {
                    dcel_.insert_edge(twin_prev, h2->twin());
                }

                ++j;
                
                while(j<boundary.rows() && !row_in_matrix(boundary.row(j), boundary_vertices)){
                    e = h1;
                    std::cout << "in while " << boundary.row(j).transpose() << std::endl;
                    if(dcel_.find_node(boundary.row(j))){
                        j++;
                        continue;
                    }
                    else{
                    if(e->on_boundary())
                        m = dcel_.insert_node(node_t(dcel_.n_nodes(), true, boundary.row(j)));
                    else
                        m = dcel_.insert_node(node_t(dcel_.n_nodes(), false, boundary.row(j)));
                    
                    prev = e->prev();
                    next = e->next();
                    twin_prev = e->twin()->prev();
                    h1 = dcel_.emplace_halfedge_(m, true);
                    h1->set_cell(next->cell());
                    dcel_.insert_edge(next, h1);
                    h1->twin()->set_subsegment(true);
                    h2 = dcel_.insert_edge(prev->next(), h1);
                    h2->set_subsegment(true);
                    h2->twin()->set_subsegment(true);
                    dcel_.insert_edge(prev,h1);
                    dcel_.remove_edge(e);
                    if(!h1->on_boundary()) {
                        dcel_.insert_edge(twin_prev, h2->twin());
                    }

                    ++j;
                    }
                }

            }
            else {++j; ++k;}
        }*/
        int j = 0;
        int k = 0;
        coords_t last_coords = boundary_vertices.row((boundary_vertices.rows() - 1)).transpose();

        for (int i = 0; i < boundary.rows(); i = j) {
            if (!row_in_matrix(boundary.row(i), boundary_vertices)) {
                j = i;

                // n1 = precedente punto del bordo completo
                coords_t n1 = last_coords;

                // n2 = punto successivo (già nel bordo)
                coords_t n2 = boundary_vertices.row(k % boundary_vertices.rows()).transpose();
                std::cout << "n1: " << n1.transpose() << " n2: " << n2.transpose() << std::endl;
                // p = punto intermedio da inserire
                coords_t p = boundary.row(j).transpose();
                std::cout << "p: " << p.transpose() << std::endl;

                if (dcel_.find_node(p)) {
                    last_coords = p;
                    ++j;
                    //if(!dcel_.find_node(p)->on_boundary())
                    //    ++k;
                    continue;
                }

                // Trova e rimuovi l’edge che stava tra n1 e n2
                halfedge_t* e = dcel_.find_halfedge_between(n1, n2);
                if (!e) e = dcel_.find_halfedge_between(n1, boundary.row(0)); // fallback (anche se dovrebbe essere orientato)
                if (!e) {
                    std::cerr << "Errore: edge mancante tra " << n1.transpose() << " e " << n2.transpose() << std::endl;
                    return;
                }

                node_t* m = dcel_.insert_node(node_t(dcel_.n_nodes(), e->on_boundary(), p));

                // Salva i puntatori
                halfedge_t* prev = e->prev();
                halfedge_t* next = e->next();
                halfedge_t* twin_prev = e->twin()->prev();

                // Inserisci i due nuovi lati
                halfedge_t* h1 = dcel_.emplace_halfedge_(m, true);
                h1->set_cell(next->cell());
                dcel_.insert_edge(next, h1);
                h1->twin()->set_subsegment(true);

                halfedge_t* h2 = dcel_.insert_edge(prev->next(), h1);
                h2->set_subsegment(true);
                h2->twin()->set_subsegment(true);

                dcel_.insert_edge(prev, h1);
                dcel_.remove_edge(e);

                if (!h1->on_boundary()) {
                    dcel_.insert_edge(twin_prev, h2->twin());
                }

                last_coords = p;
                ++j;

                // gestisci altri punti consecutivi da inserire tra n1 e n2
                while (j < boundary.rows() && !row_in_matrix(boundary.row(j), boundary_vertices)) {
                    coords_t p = boundary.row(j).transpose();
                    if (dcel_.find_node(p)) {
                        last_coords = p;
                        ++j;
                        continue;
                    }
                    
                    std::cout << "p: " << p.transpose() << std::endl;

                    node_t* m = dcel_.insert_node(node_t(dcel_.n_nodes(), h1->on_boundary(), p));

                    e=h1;
                    prev = e->prev();
                    next = e->next();
                    twin_prev = e->twin()->prev();

                    h1 = dcel_.emplace_halfedge_(m, true);
                    h1->set_cell(next->cell());
                    dcel_.insert_edge(next, h1);
                    h1->twin()->set_subsegment(true);

                    h2 = dcel_.insert_edge(prev->next(), h1);
                    h2->set_subsegment(true);
                    h2->twin()->set_subsegment(true);
                    dcel_.insert_edge(prev, h1);
                    dcel_.remove_edge(e);

                    if (!h1->on_boundary()) {
                        dcel_.insert_edge(twin_prev, h2->twin());
                    }

                    last_coords = p;
                    ++j;
                }

            } else {
                last_coords = boundary.row(i).transpose();
                ++j;
                ++k;
            }
        }

    }


    // Function to insert a vertex handling conflicts
    void insert_vertex_at_conflict(node_t* u) {
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
                mark_cavity(u, t->halfedge()->next(), D, C);
                mark_cavity(u, t->halfedge()->prev(), D, C);
                mark_cavity(u, t->halfedge()->twin()->next(), D, C);
                mark_cavity(u, t->halfedge()->twin()->prev(), D, C);
            found_on_edge = true;
        }
    
        else if (fdapde::internals::contains(u->coords(), t2, t3)) {
            D.push_back(t->halfedge()->next());
                mark_cavity(u, t->halfedge(), D, C);
                mark_cavity(u, t->halfedge()->prev(), D, C);
                mark_cavity(u, t->halfedge()->next()->twin()->next(), D, C);
                mark_cavity(u, t->halfedge()->next()->twin()->prev(), D, C);
            found_on_edge = true;
        }
    
        else if (fdapde::internals::contains(u->coords(), t3, t1)) {
            D.push_back(t->halfedge()->prev());
                mark_cavity(u, t->halfedge()->next(), D, C);
                mark_cavity(u, t->halfedge(), D, C);
                mark_cavity(u, t->halfedge()->prev()->twin()->next(), D, C);
                mark_cavity(u, t->halfedge()->prev()->twin()->prev(), D, C);
            found_on_edge = true;
        }
        //if the point did not fall on any edge then we normally try to expand the cavity from the edge of the triangle 
        if (!found_on_edge) {
            mark_cavity(u, t->halfedge(), D, C);
            mark_cavity(u, t->halfedge()->next(), D, C);
            mark_cavity(u, t->halfedge()->prev(), D, C);
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
            dcel_.remove_edge(h);
        }
        // creating the new cells 
        for (halfedge_t* h : C) { 
            add_triangle(h, std::vector<node_t*> {u});
        }
        
        // Reassigning the conflicts to the new cells  
        for (node_t* y : invalidated_nodes) {
            detect_conflicts(y, C);
        }
    } 

    //key function for the implementation of the conflict algorithm 
    void detect_conflicts(node_t* n, const std::vector<halfedge_t*>& cells_to_check = {}) {
        bool found = false;
        // initialization case: scan all the cells
        if (cells_to_check.empty()) {  
            for (auto it = dcel_.cells_begin(); it != dcel_.cells_end(); ++it) {
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

    // function to mark the cavity during insertion
    void mark_cavity(node_t* u, halfedge_t* h, std::vector<halfedge_t*>& D, std::vector<halfedge_t*>& C) {
        //if we reach the boundary we automatically create the triangle
        //if(h->on_boundary()){
        if(h->is_subsegment()){
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


    halfedge_t* add_triangle(halfedge_t* v,const std::vector<node_t*>& node){
        return dcel_.add_polygon(v ,node);
    }


    //function perfoming the flip alghoritm to tranform every non-Delaunay triangulation into a Delaunay one
    //not used for the actual costruction for the O(n^2) complexity
    void flip() {

        // creating a list of halfedges to check whether they are locally delaunay or not (in this case flippable)
        std::unordered_set<halfedge_t*> halfedges_to_check;
        for (auto it = dcel_.halfedges_begin(); it != dcel_.halfedges_end(); ++it,++it) {
            //we can already exclude the edges of the triangulation being automatically locally delaunay
            if(!it->is_subsegment()) {
                halfedges_to_check.insert(&(*it));
                std::cout << "Halfedge to check: " << it->id() << std::endl;
            }
        }
        // flip algorithm
        while (!halfedges_to_check.empty()) {
            halfedge_t* edge = *halfedges_to_check.begin();
            halfedges_to_check.erase(edge);
            if(halfedges_to_check.count(edge->twin()) > 0) 
                halfedges_to_check.erase(edge->twin()); 
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
                dcel_.remove_edge(edge);
                halfedge_t* new_edge = dcel_.insert_edge(e->prev(), e->twin()->prev());
                //std::cout << "FLIP" << std::endl;
            
                if (new_edge) {
                    // Inserting the new halfedges created by the flip into the list to check
                    if(!new_edge->prev()->is_subsegment())
                        halfedges_to_check.insert(new_edge->prev());
                    if(!new_edge->next()->is_subsegment())
                        halfedges_to_check.insert(new_edge->next());
                    if(!new_edge->twin()->prev()->is_subsegment())
                        halfedges_to_check.insert(new_edge->twin()->prev());
                    if(!new_edge->twin()->next()->is_subsegment())
                        halfedges_to_check.insert(new_edge->twin()->next());
                }
            }
        }
    }

    
    

//the following are a series of function needed to perform the Ruppert refinement algorithm
    
    //function neeeded to perform the test of encroachment for an edge e of the triangulation
    bool check_encroachment(halfedge_t* e) {
        //extracting verticies of the adjacent triangle to the edge e 
        coords_t A = e->node()->coords();
        coords_t B = e->twin()->node()->coords();
        coords_t C = e->prev()->node()->coords();         
        //computing the opposite angle to the edge e 
        double angle = fdapde::internals::angle_between(B, C, A);
        return angle >= 90.0;
    }
    
    //modified function needed to test the encorachment with circumencenter without actually inserting it
    bool check_encroachment(halfedge_t* e, const coords_t& c) {
        //extracting verticies of the adjacent triangle to the edge e 
        coords_t A = e->node()->coords();
        coords_t B = e->twin()->node()->coords();         
        //computing the opposite angle to the edge e 
        double angle = fdapde::internals::angle_between(B, c, A);
        return angle >= 90.0;
    }

    // function that checks if edge defined by halfedge h is seditious
    // if it is, it is the triangle's shortest edge since its oppoing angle is < 60 degrees and the triangle is isosceles
    bool is_edge_seditious(halfedge_t* h) {
        if (h->is_subsegment()) return false;  // boundary edges can't be seditious
        if (!(h->prev()->is_subsegment() && h->next()->is_subsegment())) return false;  

        coords_t a = h->node()->coords();
        coords_t b = h->next()->node()->coords();
        coords_t c = h->prev()->node()->coords();
        // angle in c^ is angle of interest
        if (fdapde::internals::angle_between(b,c,a) >= 60) return false; // angle is not too small
        
        // 2 edges of h's cell need to have same length and to be midpoints of another segment 
        if (fdapde::internals::segment_length(c,a) != fdapde::internals::segment_length(c,b) )  return false;
        coords_t d = h->prev()->twin()->prev()->node()->coords();  
        coords_t e = h->next()->twin()->next()->next()->node()->coords();
        if ( !( fdapde::internals::collinear(c, a, d) && 
                fdapde::internals::collinear(c, b, e) ) )
                return false;
        
        std::cout << "SEDITIOUS EDGE: " << h->id() << std::endl;
        return true;
    }   

    // Attempts to split the first non-seditious encroached edge.
    // Returns true if a segment was successfully split.
    bool split_first_encroached_segment(std::unordered_set<halfedge_t*>& boundary_edges, std::unordered_set<halfedge_t*>& encroached_edges,
            std::unordered_set<cell_t*>& bad_triangles, double rho_bar) {
        // Iterate over all currently encroached edges
        for (auto it = encroached_edges.begin(); it != encroached_edges.end(); ) {
            halfedge_t* e = *it;
            it = encroached_edges.erase(it);  // remove the edge from the set to avoid reprocessing
            encroached_edges.erase(e->twin());  // also remove the twin edge (subsegment case)
            if (!e) continue;
            std::cout << "Checking encroached edge: " << e->id() << std::endl;
            // Proceed only if none of the adjacent edges are "seditious" (i.e., dangerous to split)
            if (!is_edge_seditious(e->next()) && !is_edge_seditious(e->prev()) &&
                !is_edge_seditious(e->next()->twin()) && !is_edge_seditious(e->prev()->twin())) {
                
                // Perform the segment split and update the DCEL
                split_subsegment(e, boundary_edges, encroached_edges, bad_triangles, rho_bar);
                return true;  // Stop after the first successful split
            }
        }

        return false;  // no suitable edge was found for splitting
    }

    // function to split an encroached subsegment 
    void split_subsegment(halfedge_t* e, std::unordered_set<halfedge_t*>& boundary_edges, std::unordered_set<halfedge_t*>& encroached_edges, 
        std::unordered_set<cell_t*>& bad_triangles, double rho_bar) {
        
        boundary_edges.erase(e);  // remove the edge from the boundary edges set
        boundary_edges.erase(e->twin());  // also remove the twin edge (subsegment case)

        // Extract the two endpoints of the edge e and the third vertex forming the adjacent triangle
        node_t* a = e->node();                
        node_t* b = e->twin()->node();      
        node_t* c = e->prev()->node();      
        coords_t split_pt;                    // the point where the segment will be split

        // Case 1: acute angle at vertex b and e->next is on boundary
        if (e->next()->is_subsegment() && fdapde::internals::is_angle_acute(a->coords(), b->coords(), c->coords())) {
            coords_t split_pt_ref = c->coords();
            double r = (split_pt_ref - b->coords()).norm();  
            coords_t ab = a->coords() - b->coords();
            double L = ab.norm();                           
            // Ensure r does not exceed the length of the segment
            if (r > L) r = 0.5 * L;
            double t = r / L;
            // compute split point at distance r from b on segment ab
            split_pt = b->coords() + t * ab;
            std::cout << " Splitting segment at acute angle: " << split_pt.transpose() << std::endl;
        } 
        // Case 2: regular midpoint split
        else {
            split_pt = 0.5 * (a->coords() + b->coords());     // midpoint of segment ab
        }
        // Insert the new node into the DCEL
        node_t* m;
        if(e->on_boundary()) 
            m = dcel_.insert_node(node_t(dcel_.n_nodes(), true, split_pt));
        else
            m = dcel_.insert_node(node_t(dcel_.n_nodes(), false, split_pt));

        halfedge_t* prev = e->prev();
        halfedge_t* twin_prev = e->twin()->prev();
        halfedge_t* next = e->next();
        
        // Remove both affected triangles from the set of bad triangles
        bad_triangles.erase(e->cell());
        bad_triangles.erase(e->twin()->cell());

        halfedge_t* h1 = dcel_.emplace_halfedge_(m, true);
        h1->set_cell(next->cell());
        dcel_.insert_edge(next, h1);
        h1->twin()->set_subsegment(true);
        halfedge_t* h2 = dcel_.insert_edge(prev->next(), h1);
        h2->set_subsegment(true);
        h2->twin()->set_subsegment(true);
        dcel_.insert_edge(prev,h1);
        // Remove the encroached edge from the DCEL
        dcel_.remove_edge(e);
        
        if(is_bad_triangle(h1->cell(), rho_bar)) {
            bad_triangles.insert(h1->cell());
        }
        if(is_bad_triangle(h2->cell(), rho_bar)) {
            bad_triangles.insert(h2->cell());
        }
        
        if(!h1->on_boundary()) {
            dcel_.insert_edge(twin_prev, h2->twin());
            if(is_bad_triangle(h1->twin()->cell(), rho_bar)) {
                bad_triangles.insert(h1->twin()->cell());
            }
            if(is_bad_triangle(h2->twin()->cell(), rho_bar)) {
                bad_triangles.insert(h2->twin()->cell());
            } 
        }
        
        // Perform local flips if necessary to maintain Delaunay property
        flip_Ruppert(encroached_edges, bad_triangles, rho_bar);

        // After inserting m, recheck the two new boundary-adjacent edges for possible encroachment
        boundary_edges.insert(h1);
        boundary_edges.insert(h2);
        if (!h1->on_boundary())  {
            boundary_edges.insert(h1->twin());
        }
        if (!h2->on_boundary()) {
            boundary_edges.insert(h2->twin());
        }
        
        for(auto it = boundary_edges.begin(); it != boundary_edges.end(); ++it) {
            halfedge_t* e = *it;
            if (check_encroachment(e)) {
                if(encroached_edges.count(e) == 0)
                    encroached_edges.insert(e);
            }
        }

    }


    //function that finds the triangle containing P scanning all the triangles of the trinagulation (not used in the end)
    cell_t* find_triangle(const coords_t& P) {
        
        for (auto it = dcel_.cells_begin(); it != dcel_.cells_end(); ++it) {  
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

    cell_t* find_triangle_local(const coords_t& P, cell_t* start) {
        std::unordered_set<cell_t*> visited; //storing the cell visited since it will be pushed in the queue and not analyzed immeadiatly 
        std::queue<cell_t*> queue; //storing the cell that needs to be checked
        queue.push(start);
        visited.insert(start);

        while (!queue.empty()) {
            cell_t* current = queue.front();
            queue.pop();

            const coords_t& A = current->halfedge()->node()->coords();
            const coords_t& B = current->halfedge()->next()->node()->coords();
            const coords_t& C = current->halfedge()->prev()->node()->coords();

            if (fdapde::internals::point_in_2d_tri(P, A, B, C)) {
                return current;
            }

            for (int i = 0; i < 3; ++i) {
                halfedge_t* e = current->halfedge();
                for (int j = 0; j < i; ++j) e = e->next();
                    cell_t* neighbor = e->twin()->cell();
                if (neighbor && visited.count(neighbor) == 0) { //it exists and never visited
                    queue.push(neighbor);
                    visited.insert(neighbor);
                }
            }
        }
        return nullptr;  // fallback if not found
    }

    //function needed to perform the test of badly shaped triangle 
    //using the ratio betwwen the radius of the circumcircle and the longest edge
    //the convergence of the algorithm is proved for rho_bar >=sqrt(2)
    bool is_bad_triangle(cell_t* t, double rho_bar) {
        coords_t A = t->halfedge()->prev()->node()->coords();
        coords_t B = t->halfedge()->node()->coords();
        coords_t C = t->halfedge()->next()->node()->coords();
        double ratio = fdapde::internals::radius_edge_ratio(A, B, C);
        return ratio > rho_bar;
    }

    // Attempts to split the first bad triangle (with small angle or poor aspect ratio).
    // Returns true if a triangle was successfully split.
    bool split_first_bad_triangle(double rho_bar, std::unordered_set<halfedge_t*>& boundary_edges,
        std::unordered_set<halfedge_t*>& encroached_edges, std::unordered_set<cell_t*>& bad_triangles) {
        // Iterate over all currently bad triangles
        for (auto it = bad_triangles.begin(); it != bad_triangles.end(); ) {
            cell_t* t = *it;
            it = bad_triangles.erase(it);  // Remove the triangle from the set to avoid reprocessing

            // Skip if the triangle or its geometry is invalid
            if (!t || !t->halfedge()) continue;

            // Attempt to split the triangle by inserting its circumcenter
            bool split = split_triangle(t, boundary_edges, encroached_edges, bad_triangles, rho_bar);
            if (split)
                return true;  // Stop after the first successful split
        }

        return false;  // No suitable triangle was split
    }

    // Attempts to split a bad triangle by inserting its circumcenter.
    // If the circumcenter encroaches a segment of the PLC, it splits that segment instead.
    // Returns true if a refinement was performed.
    bool split_triangle(cell_t* t, std::unordered_set<halfedge_t*>& boundary_edges, std::unordered_set<halfedge_t*>& encroached_edges,
        std::unordered_set<cell_t*>& bad_triangles, double rho_bar) {
        // Get triangle vertices A, B, C
        coords_t A = t->halfedge()->prev()->node()->coords();
        coords_t B = t->halfedge()->node()->coords();
        coords_t C = t->halfedge()->next()->node()->coords();

        // Compute the circumcenter of triangle ABC
        coords_t c = fdapde::internals::circumcenter(A, B, C); 
        
        // Understanding where the circumcenter is actually falling 
        std::cout << "Circumcenter: " << c.transpose() << std::endl;
        std::cout << "Triangle edges: " << t->halfedge()->id() << " " << t->halfedge()->next()->id() << " " << t->halfedge()->prev()->id() << std::endl;
        cell_t* cf = find_triangle_local(c, t); 
        halfedge_t* e1 = cf->halfedge()->prev();
        halfedge_t* e2 = cf->halfedge();
        halfedge_t* e3 = cf->halfedge()->next();
        
        // Check whether the circumcenter c encroaches any boundary segment of the triangle cf
        for (halfedge_t* e : {e1, e2, e3}) {

            //if (e->on_boundary() && e->cell()) {
            if(e->is_subsegment() && e->cell()) {
                if (check_encroachment(e, c)) {
                    // Only split if the edge and its neighborhood is not marked as "seditious"
                    if (!is_edge_seditious(e->next()) && !is_edge_seditious(e->prev()) &&
                        !is_edge_seditious(e->next()->twin()) && !is_edge_seditious(e->prev()->twin())) {
                        
                        // Split the encroached subsegment instead of inserting the circumcenter
                        split_subsegment(e, boundary_edges, encroached_edges, bad_triangles, rho_bar);
                        return true;
                    }

                    // If the segment is seditious, do nothing now (will be retried later)
                    return false;
                }
            }
        }

        // If no encroachment is detected, insert the circumcenter into the mesh
        node_t* circ = dcel_.insert_node(node_t(dcel_.n_nodes(), false, c));
        // Insert the new node into the triangulation (splitting the containing triangle)
        insert_vertex(circ, cf, encroached_edges, bad_triangles, rho_bar);

        return true;
    }

     
    //same as flip function but keeping track of the encroached_edges and bad_triangles that is creating 
    void flip_Ruppert(std::unordered_set<halfedge_t*>& encroached_edges,
        std::unordered_set<cell_t*>& bad_triangles, double rho_bar) {

        std::unordered_set<halfedge_t*> halfedges_to_check;
        for (auto it = dcel_.halfedges_begin(); it != dcel_.halfedges_end(); ++it,++it) {
            //if(!it->on_boundary())
            if(!it->is_subsegment())
                halfedges_to_check.insert(&(*it));
        }

        while (!halfedges_to_check.empty()) {
            halfedge_t* edge = *(halfedges_to_check.begin());
            halfedges_to_check.erase(edge);
            if(halfedges_to_check.count(edge->twin()) > 0)
                halfedges_to_check.erase(edge->twin());

            cell_t* neighbor = edge->twin()->cell();

            coords_t A = edge->node()->coords();
            coords_t B = edge->twin()->node()->coords();
            coords_t C = edge->prev()->node()->coords();
            coords_t D = edge->twin()->prev()->node()->coords();
    
            if (fdapde::internals::in_circle(A, B, C, D) || fdapde::internals::in_circle(A, D, B, C)) {
                
                halfedge_t* e = edge;
                std::cout << "edge: " << edge->id() << std::endl;
                std::cout << "twin: " << edge->twin()->id() << std::endl;
                
                bad_triangles.erase(edge->cell());
                bad_triangles.erase(edge->twin()->cell());
                halfedge_t* prev= edge->prev();
                halfedge_t* twin_prev= edge->twin()->prev();
                
                dcel_.remove_edge(edge);
                halfedge_t* new_edge = dcel_.insert_edge(prev, twin_prev);
                std::cout <<"new edge: " << new_edge->id() << std::endl;
                std::cout << "new edge twin: " << new_edge->twin()->id() << std::endl;
                std::cout << "new edge cell: " << new_edge->cell()->id() << std::endl;
                std::cout << "new edge prev: " << new_edge->prev()->id() << std::endl;
                std::cout << "new edge prev cell: " << new_edge->prev()->cell()->id() << std::endl;
                std::cout << "FLIP RUPPERT" << std::endl;
            
                if (new_edge) {
                    if(!new_edge->prev()->is_subsegment()){
                        halfedges_to_check.insert(new_edge->prev());
                        std::cout << "prev: " << new_edge->prev()->id() << std::endl;
                    }
                    if(!new_edge->next()->is_subsegment()){
                        halfedges_to_check.insert(new_edge->next());
                        std::cout << "next: " << new_edge->next()->id() << std::endl;
                        std::cout << "next cell: " << new_edge->next()->cell()->id() << std::endl;
                    }
                    if(!new_edge->twin()->prev()->is_subsegment()){
                        halfedges_to_check.insert(new_edge->twin()->prev());
                        std::cout << "twin_prev: " << new_edge->twin()->prev()->id() << std::endl;
                    }
                    if(!new_edge->twin()->next()->is_subsegment()){
                        halfedges_to_check.insert(new_edge->twin()->next());
                        std::cout << "twin_next: " << new_edge->twin()->next()->id() << std::endl;
                    }

                    if(is_bad_triangle(new_edge->cell(), rho_bar)){
                        if (bad_triangles.count(new_edge->cell()) == 0)
                            bad_triangles.insert(new_edge->cell());
                    }
                    if(new_edge->prev()->is_subsegment() && new_edge->prev()->cell() && check_encroachment(new_edge->prev())){
                        encroached_edges.insert(new_edge->prev());
                    }
                    if(new_edge->next()->is_subsegment() && new_edge->next()->cell() && check_encroachment(new_edge->next())){
                        encroached_edges.insert(new_edge->next());
                    }
                    if(is_bad_triangle(new_edge->twin()->cell(), rho_bar)){
                        if (bad_triangles.count(new_edge->twin()->cell()) == 0)
                            bad_triangles.insert(new_edge->twin()->cell());
                    }
                   if(new_edge->twin()->prev()->is_subsegment() && new_edge->twin()->prev()->cell() && check_encroachment(new_edge->twin()->prev())){
                        encroached_edges.insert(new_edge->twin()->prev());
                    }
                    if(new_edge->twin()->next()->is_subsegment() && new_edge->twin()->next()->cell() && check_encroachment(new_edge->twin()->next())){
                        encroached_edges.insert(new_edge->twin()->next());
                    }
                }
            }
        }
    }
    
    //function working as mark_cavity but does not take into account conflits (since implement the Boyer-Watson algorithm)
    //but keeps track of the encroached edges and bad_triangles is creating 
    void dig_cavity(node_t* u, halfedge_t* vw, std::unordered_set<halfedge_t*>& encroached_edges,
        std::unordered_set<cell_t*>& bad_triangles, double rho_bar) { 
        //if we are on the boundary we add the triangle  
        if(vw->is_subsegment()){
            add_triangle(vw,std::vector<node_t*> {u});
            auto it_last = std::prev(dcel_.cells_end());  
            cell_t* t = &(*it_last);
            if(is_bad_triangle(t, rho_bar)){
                if (bad_triangles.count(t) == 0)
                    bad_triangles.insert(t);
            }
            if (check_encroachment(vw)) {
                if (encroached_edges.count(vw) == 0)
                    encroached_edges.insert(vw);
            }
            return;
        }
        //finding the point adjacent to vw
        node_t* x = dcel_.adjacent(vw);
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
            dcel_.remove_edge(vw);
            dig_cavity(u, vx, encroached_edges, bad_triangles, rho_bar);
            dig_cavity(u, xw, encroached_edges, bad_triangles, rho_bar);
        } else {
        //passed the test,adding the triangle 
            add_triangle(vw,std::vector<node_t*> {u});
            auto it_last = std::prev(dcel_.cells_end());  
            cell_t* t = &(*it_last);
            if(is_bad_triangle(t, rho_bar)){
                if (bad_triangles.count(t) == 0)
                    bad_triangles.insert(t);
            }
            return;
        } 
    }
    
    //function working as insert_vertex_at_conflict but does not take into account conflits (since implement the Boyer-Watson algorithm)
    //but keeps track of the encroached edges and bad_triangles is creating 
    void insert_vertex(node_t* u, cell_t* triangle, std::unordered_set<halfedge_t*>& encroached_edges,
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
            dcel_.remove_edge(vw);
            dig_cavity(u, wx, encroached_edges, bad_triangles, rho_bar);
            dig_cavity(u, xv, encroached_edges, bad_triangles, rho_bar);
            dig_cavity(u, twin_next, encroached_edges, bad_triangles, rho_bar);
            dig_cavity(u, twin_prev, encroached_edges, bad_triangles, rho_bar);
            found_on_edge = true;
        } else if (fdapde::internals::contains(u->coords(), w, x)) {
            halfedge_t* twin_next = wx->twin()->next();
            halfedge_t* twin_prev = wx->twin()->prev();
            bad_triangles.erase(wx->cell());
            bad_triangles.erase(wx->twin()->cell());
            dcel_.remove_edge(wx);
            dig_cavity(u, vw, encroached_edges, bad_triangles, rho_bar);
            dig_cavity(u, xv, encroached_edges, bad_triangles, rho_bar);
            dig_cavity(u, twin_next, encroached_edges, bad_triangles, rho_bar);
            dig_cavity(u, twin_prev, encroached_edges, bad_triangles, rho_bar);
            found_on_edge = true;
        } else if (fdapde::internals::contains(u->coords(), x, v)) {
            halfedge_t* twin_next = xv->twin()->next();
            halfedge_t* twin_prev = xv->twin()->prev();
            bad_triangles.erase(xv->cell());
            bad_triangles.erase(xv->twin()->cell());
            dcel_.remove_edge(xv);
            dig_cavity(u, vw, encroached_edges, bad_triangles, rho_bar);
            dig_cavity(u, wx, encroached_edges, bad_triangles, rho_bar);
            dig_cavity(u, twin_next, encroached_edges, bad_triangles, rho_bar);
            dig_cavity(u, twin_prev, encroached_edges, bad_triangles, rho_bar);
            found_on_edge = true;
        }
    
        if (!found_on_edge) {
            dig_cavity(u, vw, encroached_edges, bad_triangles, rho_bar);
            dig_cavity(u, wx, encroached_edges, bad_triangles, rho_bar);
            dig_cavity(u, xv, encroached_edges, bad_triangles, rho_bar);
        }
        //no need to flip since Boyer-Watson mantains the Dleaunay
        //flip_Ruppert(dcel, encroached_edges, bad_triangles, rho_bar);
    }
    
};
  
}  // namespace fdapde

#endif // _DELAUNAY_H_