// This file is part of fdaPDE, a C++ library for physics-informed
// spatial and functional data analysis.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __FDAPDE_ADAPTIVITY_H__
#define __FDAPDE_ADAPTIVITY_H__

#include "header_check.h"

namespace fdapde {


//template <int LocalDim, int EmbedDim, AdaptiveStrategy Strategy>
template <int LocalDim, int EmbedDim>
class Adaptivity {
   public:
    static constexpr int local_dim = LocalDim;
    static constexpr int embed_dim = EmbedDim;

    using coords_t = Eigen::Matrix<double, 1, embed_dim>;
    using node_t = typename DCEL<local_dim, embed_dim>::node_t;
    using halfedge_t = typename DCEL<local_dim, embed_dim>::halfedge_t;
    using cell_t = typename DCEL<local_dim, embed_dim>::cell_t;
    using dcel_t = DCEL<local_dim, embed_dim>;
    using delaunay_t = Delaunay<local_dim, embed_dim>;
    using polygon_t = Polygon<local_dim, embed_dim>;

    //using storage_t = typename StrategyData<embed_dim, Strategy>::storage_t;
    using storage_t = typename Eigen::Matrix<double, Eigen::Dynamic, embed_dim>;

    template<int local_dim, int embed_dim>
    struct Metric {
        using halfedge_t = typename fdapde::DCEL<local_dim, embed_dim>::halfedge_t;
        using coords_t = Eigen::Matrix<double, 1, embed_dim>;
        using dcel_t = DCEL<local_dim, embed_dim>;

        Metric(std::function<Eigen::Matrix<double, embed_dim, embed_dim>(const coords_t&)> metric_fun, dcel_t& dcel):  dcel_(dcel) {
            if (!metric_fun) {throw std::invalid_argument("Metric function not provided.");}
            metric_fun_ = metric_fun;
        }

        Metric(std::unordered_map<node_t*, Eigen::Matrix<double, embed_dim, embed_dim>>& node_metrics, dcel_t& dcel):  dcel_(dcel), node_metrics_(node_metrics) {
            if (node_metrics.empty()) {throw std::invalid_argument("Map not provided.");}
        }

        double metric_edge_length(halfedge_t* e) const{
            node_t* A = e->node();
            node_t* B = e->twin()->node();
            const coords_t AB = B->coords() - A->coords();

            Eigen::Matrix<double,2,2> M = 0.5*metric_at_node(A) + 0.5*metric_at_node(B);
            double len2 = (AB * M * AB.transpose())(0, 0);
            double lenM = std::sqrt(std::max(0.0, len2));
            return lenM; 
        }

        double compute_metric_error() const {
            int total_edges = 0;
            int non_compliant_edges = 0;

            for (auto it = dcel_.halfedges_begin(); it != dcel_.halfedges_end(); ++it) {
                halfedge_t* e = &(*it);
                if (e->id() >= e->twin()->id()) continue;
                total_edges++;
                double lenM = metric_edge_length(e);
                if (lenM < Lmin_ || lenM > Lmax_) {
                    non_compliant_edges++;
                }
            }

            return total_edges? static_cast<double>(non_compliant_edges) / total_edges : 0.0;
        }

        void delete_nodes(std::unordered_set<node_t*> deleted_nodes){
            for(auto n: deleted_nodes){
                node_metrics_.erase(n);
            }
        }

        void add_node(node_t* n_new, node_t* n1, node_t* n2, node_t* n3){
            if(!node_metrics_.empty()){
                auto bar_coords = barycentric_coordinates_(n1->coords(), n2->coords(), n3->coords(), n_new->coords());
                node_metrics_.insert({n_new, interpolate_log_euclidean_(n_new, n1, n2, n3, bar_coords )});
            }
        }

        void update_node(node_t* n_new, coords_t old_coords, node_t* n1, node_t* n2){
            if(!node_metrics_.empty()){
                auto bar_coords = barycentric_coordinates_(old_coords, n1->coords(), n2->coords(), n_new->coords());
                node_metrics_[n_new] = interpolate_log_euclidean_(n_new, n_new, n1, n2, bar_coords );
            }
        }

        Eigen::Matrix<double, embed_dim, embed_dim> metric_at_node(node_t* n) const {
            if(!node_metrics_.empty()) return node_metrics_.at(n);
            else return metric_fun_(n->coords());
        }

        double Lmax() const { return Lmax_; }
        double Lmin() const { return Lmin_; }
        std::function<Eigen::Matrix<double, embed_dim, embed_dim>(const coords_t&)> metric_fun()  {return metric_fun_;}

        
    private:
        std::function<Eigen::Matrix<double, embed_dim, embed_dim>(const coords_t&)> metric_fun_;
        dcel_t& dcel_;
        double Lmax_ = 1.3;
        double Lmin_ = 0.8;
        std::unordered_map<node_t*, Eigen::Matrix<double, embed_dim, embed_dim>> node_metrics_ = {};

        Eigen::Matrix<double, embed_dim, embed_dim> matrix_log_(Eigen::Matrix<double, embed_dim, embed_dim> M){
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, embed_dim, embed_dim>> es(M);
            Eigen::Matrix<double, embed_dim, embed_dim> R = es.eigenvectors();
            Eigen::Matrix<double, embed_dim, 1> l = es.eigenvalues();
            Eigen::Matrix<double, embed_dim, embed_dim> L;
            L << std::log(l(0)), 0.0,
                 0.0,         std::log(l(1));
            return R * L * R.transpose();
        }

        Eigen::Matrix<double, embed_dim, embed_dim> matrix_exp_(Eigen::Matrix<double, embed_dim, embed_dim> M){
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, embed_dim, embed_dim>> es(M);
            Eigen::Matrix<double, embed_dim, embed_dim> R = es.eigenvectors();
            Eigen::Matrix<double, embed_dim, 1> e = es.eigenvalues();
            Eigen::Matrix<double, embed_dim, embed_dim> E;
            E << std::exp(e(0)), 0.0,
                 0.0,         std::exp(e(1));
            return R * E * R.transpose();
        }

        std::vector<double> barycentric_coordinates_(const coords_t& v1,const coords_t& v2,const coords_t& v3,const coords_t& x_new) {
            
            double D_tot = fdapde::internals::signed_measure_2d_tri(v1, v2, v3);
            std::vector<double> alpha(3, 0.0);
            double D1 = fdapde::internals::signed_measure_2d_tri(x_new, v2, v3);
            double D2 = fdapde::internals::signed_measure_2d_tri(x_new, v3, v1);
            double D3 = fdapde::internals::signed_measure_2d_tri(x_new, v1, v2);

            alpha[0] = D1 / D_tot; // alpha1
            alpha[1] = D2 / D_tot; // alpha2
            alpha[2] = D3 / D_tot; // alpha3

            return alpha;
        }

        // split point n_new falls on old edge v1v2
        Eigen::Matrix<double, embed_dim, embed_dim> interpolate_log_euclidean_(node_t* n_new, node_t* v1, node_t* v2, node_t* v3, std::vector<double>& alpha_coords) {
            
            double alpha1 = alpha_coords[0]; 
            double alpha2 = alpha_coords[1]; 
            double alpha3 = alpha_coords[2]; 
            
            // M(x_i) = M_i
            Eigen::Matrix<double, embed_dim, embed_dim> M1 = node_metrics_.at(v1); 
            Eigen::Matrix<double, embed_dim, embed_dim> M2 = node_metrics_.at(v2);
            Eigen::Matrix<double, embed_dim, embed_dim> M3 = node_metrics_.at(v3);
            
            Eigen::Matrix<double, embed_dim, embed_dim> L1 = matrix_log_(M1);
            Eigen::Matrix<double, embed_dim, embed_dim> L2 = matrix_log_(M2);
            Eigen::Matrix<double, embed_dim, embed_dim> L3 = matrix_log_(M3);

            // S = sum(alpha_i * ln(M_i))
            Eigen::Matrix<double, embed_dim, embed_dim> S = alpha1 * L1 + alpha2 * L2 + alpha3 * L3;
            // M(x) = exp(S)
            Eigen::Matrix<double, embed_dim, embed_dim> M_interpolated = matrix_exp_(S);
            
            return M_interpolated;
        }


    };



    Adaptivity(delaunay_t& mesh, const storage_t& data_points, std::function<Eigen::Matrix<double, embed_dim, embed_dim>(const coords_t&)> metric_fun): 
        mesh_(mesh), data_points_(data_points), M_(metric_fun, mesh_.dcel()) {   
    }

    Adaptivity(delaunay_t& mesh, std::unordered_map<node_t*, Eigen::Matrix<double, embed_dim, embed_dim>>& node_metrics): 
        mesh_(mesh), M_(node_metrics, mesh_.dcel()) {}


    // method to adapt the mesh according to a metric M, by splitting and collapsing edges
    void adaptivity_cycle(const double min_angle, const double max_area){
        
        std::cout << "Initial metric error: " << M_.compute_metric_error() << std::endl;

        // nodes-based collapse, more efficient
        int it1 = 0;
        for(auto it= mesh_.dcel_.nodes_begin(); it != mesh_.dcel_.nodes_end();){
            node_t* n = &(*it);
            it= std::next(it);
            std::unordered_set <cell_t*> modified={};
            std::unordered_set <node_t*> deleted_nodes={};
            if(n->on_boundary() || n->halfedge()->is_segment()) {
                auto star = halfedges_from_node_(n);
                halfedge_t* e;
                for(auto h: star)
                    if(h->on_boundary() || h->is_segment()){
                        e = h;
                        break;
                    }
                if(M_.metric_edge_length(e)<M_.Lmin() && !e->prev()->is_segment() && fdapde::internals::collinear(e->prev()->twin()->prev()->node()->coords(), n->coords(), e->next()->node()->coords() )){
                    deleted_nodes.insert(n);
                    modified = collapse_segment_(e, n);

                }
            }else{
                halfedge_t* he = best_collapse_target_(n);
                if(he){
                    deleted_nodes.insert(n);
                    modified = collapse_edge_(he, n);
                }
            }
            if(!modified.empty()){
                ++it1;
                M_.delete_nodes(deleted_nodes);
                //refinenement_local_(min_angle, max_area, modified);
            }

        }
        std::cout << "it1: " << it1 << std::endl;
        
        auto it  = mesh_.dcel_.halfedges_begin();
        auto end = mesh_.dcel_.halfedges_end();
        bool split_done = false;
        int it2=0;
        do {
            ++it2;
            split_done = false;
            for (auto it = mesh_.dcel_.halfedges_begin(); it != mesh_.dcel_.halfedges_end(); ++it) {
                halfedge_t* e = &(*it);
                if (e->id() >= e->twin()->id()) continue;
                double id = M_.metric_edge_length(e);
                std::unordered_set <cell_t*> modified={};
                if (M_.metric_edge_length(e) > M_.Lmax()) {
                    if(!e->is_segment())
                        modified = split_edge_(e);
                    else
                        modified = split_segment_(e);
                    if(!modified.empty()){
                        split_done = true;
                        //refinenement_local_(min_angle, max_area, modified);
                        break; 
                    }
                }
            }
        } while (split_done); 
        std::cout << "it2: " << it2 << std::endl;

        smoothing_(1.8, 3);  // smoothing routine
        
        std::cout << "Error: "<< M_.compute_metric_error() << std::endl;   // after refinement metric is not updated
        mesh_.check_quality_(min_angle, max_area);
        mesh_.refinement(min_angle, max_area);

        return;
    }

    struct CollapseSimulate {
        std::optional<delaunay_t> temp_dcel;
        std::unordered_map<int, std::pair<bool,int>> cavity_info;          // <halfedge on boundary, vertex degree>
        std::unordered_map<node_t*, Eigen::Matrix<double, embed_dim, embed_dim>> temp_node_metrics={}; // metrics at nodes in the temporary dcel
        bool valid() const { return !cavity_info.empty(); }
    };

    dcel_t&  dcel()  { return mesh_.dcel_; }
    delaunay_t& delaunay() { return mesh_; }


   private:
    delaunay_t& mesh_;
    const storage_t data_points_; // data points
    Metric<local_dim, embed_dim> M_;


    void refinenement_local_(double min_angle, double max_area, std::unordered_set<cell_t*> modified){
        std::multimap<double, cell_t*> bad_triangles = {};
        std::unordered_set<halfedge_t*> segments = {}, encroached = {};
        
        for(auto it = modified.begin(); it != modified.end();++it){
            auto c = *it;

            auto h_start = c->halfedge();
            auto h = h_start;
            do{
                if(h->is_segment()){
                    segments.insert(h);
                    if(mesh_.check_encroachment_(h)){
                        encroached.insert(h);
                    }
                }
                h=h->next();
            }while(h!=h_start);

            auto pr = mesh_.is_bad_triangle_(c, min_angle, max_area);
            if(pr>=0 && std::find_if(bad_triangles.begin(), bad_triangles.end(), [c](const auto& entry) { return entry.second == c; }) == bad_triangles.end()){ 
                bad_triangles.insert({pr,c});
            }
        }

        while(true){  // local refinement cycle
            if (mesh_.split_first_encroached_segment_(segments, encroached, bad_triangles, min_angle, max_area))
                continue;
            // otherwise, attempt to split the worst triangle in order of priority
            if (mesh_.split_first_bad_triangle_(segments, encroached, bad_triangles, min_angle, max_area)) 
                continue;
            break;
        }
    }

    std::list<halfedge_t*> halfedges_from_node_(node_t* n) {
        std::list<halfedge_t*> star;
        if(n->halfedge()->cell()==nullptr) n->set_halfedge(n->halfedge()->twin()->next()); // in case the boundary node starts from an external halfedge
        halfedge_t* h_start = n->halfedge();
        halfedge_t* h = h_start;
        do {
            star.push_back(h);
            h = h->prev()->twin();
            if(!h->cell() && h->on_boundary()) {h=h->prev()->twin(); } // in case of boundary node don't count the external halfedge
        } while (h != h_start);
        return star;
    }

    halfedge_t* best_collapse_target_(node_t* n) {
        std::list<halfedge_t*> star = halfedges_from_node_(n); 
        if(star.empty()) return nullptr;

        // order halfedges based on metric length
        std::multimap<double, halfedge_t*> sorted = {};
        for (auto h : star) {
            double length = M_.metric_edge_length(h);
            if(length<M_.Lmin()) sorted.insert({length, h});
        }

        for(auto pair_l_h : sorted) {
            auto result = simulate_collapse_(pair_l_h.second, n);
            if(!result.valid()) continue; 
            if(!collapse_cost_(result)) continue;
            else return pair_l_h.second;
        }

        return nullptr;
    }

    bool collapse_cost_(CollapseSimulate& result) {
        std::unordered_set<halfedge_t*> elems_modified;
        for(auto it = result.temp_dcel->dcel_.halfedges_begin(); it != result.temp_dcel->dcel_.halfedges_end(); ++it, ++it)
            elems_modified.insert(const_cast<halfedge_t*>(&*it));

        double meanL = 0.0;
        double maxL = 0.0;
        Metric<local_dim, embed_dim> M_temp(result.temp_node_metrics, result.temp_dcel->dcel_);
        for(auto h : elems_modified){
            double L = M_temp.metric_edge_length(h);
            meanL += L;
            maxL = std::max(maxL, L);
        }

        // reject if violates upper length bound
        if(maxL > M_.Lmax())
            return false;

        return true;
    }


    // collapse edge e into node; if node==nullptr, insert the midpoint of the edge
    std::unordered_set<cell_t*> collapse_edge_(halfedge_t* e, node_t* node) {
        std::unordered_set<cell_t*> elems_modified;

        if ((node && node->on_boundary()) || (!node && (e->node()->on_boundary() || e->twin()->node()->on_boundary()))) 
            return elems_modified;
        if (e->is_segment() && !node) // don't collapse segment into its midpoint
            return elems_modified;
        if (e->is_segment() && 
               (node == e->node() && !e->prev()->is_segment() && fdapde::internals::collinear(e->prev()->twin()->prev()->node()->coords(), node->coords(), e->next()->node()->coords()) ) ||
               (node == e->twin()->node() && !e->next()->is_segment()&& fdapde::internals::collinear(e->node()->coords(), node->coords(), e->next()->twin()->prev()->node()->coords() )  ) )   
                return collapse_segment_(e,node);

                
        //std::cout << "Collapsing edge " << e->id() << std::endl;

        // identify the cavity edges
        std::list<halfedge_t*> polygonal_cavity;
        std::list<halfedge_t*> edges_to_remove;

        if (node) { // one of e's extrema
            halfedge_t* h_start = node->halfedge();
            halfedge_t* h = h_start;
            do {
                edges_to_remove.push_back(h);
                polygonal_cavity.push_back(h->next());
                h = h->prev()->twin();
            } while (h != h_start);

        } else { // midpoint case          NON LO FACCIO PIU IN REALTA, DA TOGLIERE..........
            edges_to_remove.push_back(e);
            halfedge_t* h = e->prev()->twin();
            halfedge_t* limit1 = e->twin()->next();
            halfedge_t* limit2 = e->next();
            do {
                edges_to_remove.push_back(h);
                polygonal_cavity.push_back(h->next());
                h = h->prev()->twin();
            } while(h != limit1);
            edges_to_remove.push_back(limit1);  // otherwise it would be skipped
            h = h->next()->twin();
            do {
                edges_to_remove.push_back(h);
                polygonal_cavity.push_back(h->next());
                h = h->prev()->twin();
            } while(h != limit2);
            edges_to_remove.push_back(limit2);  // otherwise it would be skipped
        }

        node_t* u = e->node();
        node_t* v = e->twin()->node();
        coords_t midpoint = 0.5 * (u->coords() + v->coords());
        // empty the cavity
        for (auto he : edges_to_remove){
            mesh_.dcel_.remove_edge(he);
        }
        if (node) {
            mesh_.dcel_.remove_node(node);
        } else {
            mesh_.dcel_.remove_node(u);
            mesh_.dcel_.remove_node(v);
        }

        if(polygonal_cavity.size()>3){
            // triangulate the cavity using the polygon.h algorithm
            Eigen::Matrix<double, Eigen::Dynamic, 2> cavity_coords(polygonal_cavity.size(), 2);
            int i = 0;
            for (auto h : polygonal_cavity){
                cavity_coords.row(i++) = h->node()->coords();
            }
            delaunay_t temp_dcel({cavity_coords});  // need to create a dcel that deals with collinear cavity nodes cases
            auto triangulation = temp_dcel.dcel_.template to_triangulation<Triangulation<LocalDim, EmbedDim>>();
            mesh_.dcel_.from_triangulation(triangulation, {});
        }

        if (!node) { // add the midpoint to the triangulation
            auto node_triangle = mesh_.find_triangle_local_(midpoint, polygonal_cavity.front()->cell());   // VEDI SE CAMBIARE CON SOLO CELLE DELLA CAVITà
            node = mesh_.dcel_.insert_node(node_t(std::prev(mesh_.dcel_.nodes_end())->id() + 1, false, midpoint));
            mesh_.insert_vertex_(node, node_triangle);
        }

        // collect modified elements; since the mesh is Delaunay, only the internal cavity cells are modified (similar logic to Bowyer-Watson)
        for (auto cavity_edge: polygonal_cavity) {  // cavity edges are not modified    
            if (cavity_edge->cell()) {
                elems_modified.insert(cavity_edge->cell());
            }
            if (cavity_edge->prev()->twin()->cell()) {
                elems_modified.insert(cavity_edge->prev()->twin()->cell());
            }
        }
        return elems_modified;
    }


    // collapse edge e into node (which is one of e's extremes), when e is a boundary segment
    std::unordered_set<cell_t*> collapse_segment_(halfedge_t* e, node_t* node) {
        std::unordered_set<cell_t*> elems_modified = {};
        if (!e->is_segment()) return elems_modified;
        if (!e->cell()) e = e->twin();
        bool on_boundary = e->on_boundary();

        // identify the cavity edges
        std::vector<std::list<halfedge_t*>> polygonal_cavity(2);
        std::list<halfedge_t*> edges_to_remove;

        halfedge_t* h = node->halfedge();
        bool b = 0;
        int segment_count = 0;
        do {
            if(!h->is_segment())
                edges_to_remove.push_back(h);
            else segment_count++;
            polygonal_cavity[b].push_back(h->next());
            h = h->prev()->twin();
            if(h->is_segment() && !on_boundary){
                polygonal_cavity[b].push_back(h->next());
                b = !b;
            }
        } while (h != node->halfedge() );
        // cannot collapse if more than the two collinear segments depart from node
        if(segment_count > 2) return elems_modified;
        //std::cout << "Collapsing segment " << e->id() << std::endl;

        for(auto re: edges_to_remove){
            mesh_.dcel_.remove_edge(re);
        }
        
        halfedge_t* prev = e->prev();
        halfedge_t* next = e->next();
        node_t* u = nullptr;
        
        cell_t* new_cell;
        if(node==e->node()){   
            u = prev->node();
            h = mesh_.dcel_.insert_edge(prev, next);
            mesh_.dcel_.remove_edge(e);
            mesh_.dcel_.remove_edge(prev);
            new_cell = next->cell();
            
        }
        else{ // node = e->twin()->node()
            u = e->node();
            h = mesh_.dcel_.insert_edge(e, next->next());
            mesh_.dcel_.remove_edge(e);
            mesh_.dcel_.remove_edge(next);
            new_cell = prev->cell();
        }
        mesh_.dcel_.remove_node(node);
        u->set_halfedge(h);
        h->set_segment(true);
        h->twin()->set_segment(true);
        h->set_cell(new_cell);
        
        if(on_boundary)
            h->twin()->set_cell(nullptr);
        else
            h->twin()->set_cell(h->twin()->next()->cell());

        for(auto& polygonal_cavity_i : polygonal_cavity){
            if(polygonal_cavity_i.size()>3){
                // triangulate the cavity using the polygon.h algorithm     // DEVO CAPIRE COME TENERE CONTO DEL BORDO !!!! 
                Eigen::Matrix<double, Eigen::Dynamic, 2> cavity_coords(polygonal_cavity_i.size(), 2);
                int j = 0;
                for (auto h : polygonal_cavity_i){
                    cavity_coords.row(j++) = h->node()->coords();
                }
                delaunay_t temp_dcel({cavity_coords});  
                auto triangulation = temp_dcel.dcel_.template to_triangulation<Triangulation<LocalDim, EmbedDim>>();
                mesh_.dcel_.from_triangulation(triangulation, {});
            }
            // collect modified elements; since the mesh is Delaunay, only the internal cavity cells are modified (similar logic to Bowyer-Watson)
            for (auto cavity_edge: polygonal_cavity_i) {  // cavity edges are not modified
                if (cavity_edge->cell()) 
                    elems_modified.insert(cavity_edge->cell());
            }
        }

        return elems_modified;
    }


    // simulate the collapse of an edge and return the modified elements (in a temporary DCEL) and cavity info
    CollapseSimulate simulate_collapse_(halfedge_t* e, node_t* node) {
        CollapseSimulate result;

        // check if collapse is topologically possible                    
        // inserting the midpoint is not possible if one of e's vertices is at the boundary (boundary would be squeezed)
        if(e->is_segment() || (node && node->on_boundary()) || (!node && (e->node()->on_boundary() || e->twin()->node()->on_boundary())) )
            return result;

        // identify the polygonal cavity and the elements' cavity info
        std::list<halfedge_t*> polygonal_cavity;
        std::unordered_map<int, std::pair<bool,int>> cavity_info;  // <halfedge on boundary, vertex degree>

        if (node) { // extreme node option
            halfedge_t* h_start = node->halfedge();  // don't know if e is the edge going out of node or its twin
            halfedge_t* h = h_start;
            do{
                polygonal_cavity.push_back(h->next());
                cavity_info.insert({h->next()->id(), {h->next()->on_boundary(), count_triangles_from_vertex_(h->next()->node())}});
                h = h->prev()->twin();
            }while(h!=h_start);
        } else { // midpoint option
            halfedge_t* h = e->prev()->twin();
            halfedge_t* limit1= e->twin()->next();
            halfedge_t* limit2= e->next();
            do{
                polygonal_cavity.push_back(h->next());
                cavity_info.insert({h->next()->id(), {h->next()->on_boundary(), count_triangles_from_vertex_(h->next()->node())}});
                h = h->prev()->twin();
            }while(h!= limit1);

            h = h->next()->twin();
            do{
                polygonal_cavity.push_back(h->next());
                cavity_info.insert({h->next()->id(), {h->next()->on_boundary(), count_triangles_from_vertex_(h->next()->node())}});
                h = h->prev()->twin();
            }while(h!= limit2);
        } 
        
        // create a small DCEL cavity to simulate the collapse
        Eigen::Matrix<double, Eigen::Dynamic, 2> cavity_coords(polygonal_cavity.size(),2); 
        int i=0;
        for(auto h: polygonal_cavity){
            cavity_coords.row(i++) = h->node()->coords();
        }
        // min size, otherwise collapse is not possible
        if ((!node && cavity_coords.rows() < 4) || (node && cavity_coords.rows() < 3))
            return result;

        // if the cavity has duplicate points, it means some incident edges to u,v share a common node besides the edges on the two sides of e (the ones forming the 2 adjacent cells to u) 
        // therefore, the collapse is not possible since it would delete another mesh node besides u and v
        if (!node) {
            auto has_duplicates = [&](const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& m) {
                std::set<std::pair<double,double>> seen;
                for (int i=0; i<m.rows(); ++i) {
                    Eigen::Vector2d pi = m.row(i);
                    auto key = std::make_pair(pi(0), pi(1));
                    if (!seen.insert(key).second) return true;
                }
                return false;
            };
            if (has_duplicates(cavity_coords)) return result;  
        }
        
        // build a temporary DCEL cavity to simulate the collapse
        std::unordered_set<cell_t*> elems_modified;
        if(node){
            result.temp_dcel.emplace({cavity_coords});  // Delaunay cavity
            //make the cavity information coherent with the actual mesh information (create a mapping)
            node_t* start_node = result.temp_dcel->dcel_.find_node(cavity_coords.row(0));
            auto star = halfedges_from_node_(start_node);
            halfedge_t* h ;  // external boundary (its cell is null)
            for(auto he: star){
                if(he->on_boundary()){
                    h = he->twin();
                    break;
                }
            }
            for(auto cavity_edge: polygonal_cavity){
                h->twin()->set_id(cavity_edge->id());
                result.temp_node_metrics.insert({h->twin()->node(), M_.metric_at_node(cavity_edge->node())});
                h = h->prev(); //to go in counterclockwise order
            }        
        }
        else{
            coords_t midpoint = 0.5*(e->node()->coords() + e->twin()->node()->coords());
            Eigen::Matrix<double, 1, 2> cavity_coords_mid(1, 2);
            cavity_coords_mid.row(0) = midpoint;
            result.temp_dcel.emplace({cavity_coords}, cavity_coords_mid);  // Delaunay cavity

            // make the cavity information coherent with the actual mesh information (create a mapping)
            node_t* start_node = result.temp_dcel->dcel_.find_node(cavity_coords.row(0));
            auto star = halfedges_from_node_(start_node);
            halfedge_t* h ;  // external boundary (its cell is null)
            for(auto he: star){
                if(he->on_boundary()){
                    h = he->twin();
                    break;
                }
            }
            for(auto cavity_edge: polygonal_cavity){
                h->twin()->set_id(cavity_edge->id());
                result.temp_node_metrics.insert({h->twin()->node(), M_.metric_at_node(cavity_edge->node())});
                h = h->prev(); // to go in counterclockwise order
            } 
        }        

        result.cavity_info = cavity_info;
        return result;
    }

    // count the number of triangles incident to a vertex (same as the number of edges radiating from said vertex)
    int count_triangles_from_vertex_(node_t* v) const {
        int count=0;
        halfedge_t* h_start = v->halfedge();
        halfedge_t* h = h_start;
        do {
            count++;
            h = h->prev()->twin();
            if(!h->cell()){ // on boundary
               if(h==h_start || h->twin()==h_start) break; 
               h = h->prev()->twin();  // skip the boundary (no cell in between halfedges)
            }
        } while (h != h_start);
        return count;
    }


    // split edge e into two edges by inserting a midpoint
    std::unordered_set<cell_t*> split_edge_(halfedge_t* e){

        //std::cout << "Splitting edge " << e->id() << std::endl;
        node_t* A = e->node();
        node_t* B = e->twin()->node();
        node_t* C = e->prev()->node();

        // coordinates of the point to insert to split e
        //coords_t split_pt = 0.5 * (A->coords() + B->coords());
        coords_t AB = B->coords() - A->coords();
        double lA = std::sqrt(AB * M_.metric_at_node(A) * AB.transpose());
        double lB = std::sqrt(AB * M_.metric_at_node(B) * AB.transpose());
        double w= 1/(1+std::sqrt(lB/lA));
        coords_t split_pt = A->coords() + w*AB;
        if(e->is_segment())
            return split_segment_(e);
        // insert the split point into the dcel_
        node_t* m = mesh_.dcel_.insert_node(node_t(std::prev(mesh_.dcel_.nodes_end())->id()+1, e->is_segment(), split_pt));
        
        halfedge_t* prev = e->prev(); // CA
        halfedge_t* next = e->next();  // BC  
        halfedge_t* twin_next = e->twin()->next();  // AD          
        halfedge_t* twin_prev = e->twin()->prev();  // DB
        
                       
        // mimic mesh_.insert_vertex_ logic (already know the point falls on edge e)
        mesh_.dcel_.remove_edge(e);
        mesh_.dig_cavity_(m, prev);
        mesh_.dig_cavity_(m, next);
        mesh_.dig_cavity_(m, twin_next);
        mesh_.dig_cavity_(m, twin_prev);

        std::unordered_set<cell_t*> modified_cells;
        halfedge_t* h = m->halfedge();
        do {
                if(h->cell())  // boundary case
                   modified_cells.insert(h->cell());
                h = h->prev()->twin();
        } while(h!=m->halfedge());

        M_.add_node(m, A, B, C);
        return modified_cells;
    }

    // split segment e into two edges by inserting a midpoint
    std::unordered_set<cell_t*> split_segment_(halfedge_t* e){

        std::unordered_set<cell_t*> modified_cells;
        if (!e->is_segment()) return modified_cells;
        if(!e->cell()) e = e->twin();
        bool on_boundary = e->on_boundary();
        //std::cout << "Splitting segment " << e->id() << std::endl;

        node_t* A = e->node();
        node_t* B = e->twin()->node();
        node_t* C = e->prev()->node();

        // coordinates of the point to insert to split e
        //coords_t split_pt = 0.5 * (A->coords() + B->coords());
        coords_t AB = B->coords() - A->coords();
        double lA = std::sqrt(AB* M_.metric_at_node(A) * AB.transpose());
        double lB = std::sqrt(AB* M_.metric_at_node(B) * AB.transpose());
        double w= 1/(1+std::sqrt(lB/lA));
        coords_t split_pt = A->coords() + w*AB;
        // insert the split point into the dcel_
        node_t* M = mesh_.dcel_.insert_node(node_t(std::prev(mesh_.dcel_.nodes_end())->id() + 1, e->is_segment(), split_pt));
        
        halfedge_t* prev = e->prev(); // CA
        halfedge_t* next = e->next();  // BC  

        // create new connections between am, mb, bc before removing e
        halfedge_t* h1 = mesh_.dcel_.emplace_halfedge(M, e->is_segment());
        h1->set_cell(next->cell());
        mesh_.dcel_.insert_edge(h1, next); // MB

        M->set_halfedge(h1);
        h1->twin()->set_segment(e->is_segment());
        halfedge_t* h2 = mesh_.dcel_.insert_edge(prev->next(), h1);  // AM
        h2->set_segment(e->is_segment());
        h2->twin()->set_segment(e->is_segment());

        mesh_.dcel_.remove_edge(e);

        mesh_.dig_cavity_(M,prev);
        mesh_.dig_cavity_(M,next);

        if(!on_boundary){
            mesh_.dig_cavity_(M, h2->twin()->next());  // AD
            mesh_.dig_cavity_(M, h1->twin()->prev());  // BD
        }

        halfedge_t* h = M->halfedge();
        do {
                if(h->cell()){  // boundary case 
                   modified_cells.insert(h->cell());
                }
                h = h->prev()->twin();
        } while(h!=M->halfedge());

        M_.add_node(M, A, B, C);
        return modified_cells;
    }


    // reallocates vertexes
    void smoothing_(double omega, int nb_iter) {

        for (int iter = 0; iter < nb_iter; ++iter) {
            std::vector<coords_t> new_positions(mesh_.dcel_.n_nodes());
            
            int i=0;
            for (auto it = mesh_.dcel_.nodes_begin(); it!= mesh_.dcel_.nodes_end(); ++it, ++i) {
                node_t* n = &(*it);
                if (n->on_boundary()) {  // ignore boundary nodes so that shape doesn't change
                    new_positions[i] = n->coords();
                    continue;
                }

                // calculate the metric baricenter
                Eigen::Matrix<double, embed_dim, embed_dim> W_total = Eigen::Matrix<double, embed_dim, embed_dim>::Zero();
                Eigen::Matrix<double, embed_dim, 1> P_weighted = Eigen::Matrix<double, embed_dim, 1>::Zero();
                const Eigen::Matrix<double, embed_dim, embed_dim>& M_n = M_.metric_at_node(n);
                
                // iterate over n adjacent nodes
                std::list<node_t*> adjacent_nodes;
                halfedge_t* h = n->halfedge();
                bool b = 0;
                do {
                    adjacent_nodes.push_back(h->next()->node());
                    h = h->prev()->twin();
                } while (h != n->halfedge() );
                for (auto& nn : adjacent_nodes) {
                    node_t* N = nn;
                    
                    const Eigen::Matrix<double, embed_dim, embed_dim> W_i = M_n; 
                    W_total += W_i;
                    P_weighted += W_i * nn->coords().transpose(); // 2x1
                }

                coords_t n_optimal;  
                if (W_total.determinant() > 1e-12) {  // metric baricenter
                    n_optimal = (W_total.inverse() * P_weighted).transpose();
                } else {
                    // if W-total is not invertible, use Euclidean baricenter
                    n_optimal = P_weighted.transpose() / W_total.sum(); 
                }

                // move n towards n_optimal, through the relaxation parameter omega (weighted laplacian)
                coords_t n_new = (1.0 - omega) * n->coords() + omega * n_optimal;

                if (is_valid_movement_(n, n_new)) { 
                    //new_positions[i] = n_new;
                    coords_t old_coords = n->coords();
                    n->set_coords(n_new);
                    halfedge_t* h_start = n->halfedge();
                    halfedge_t* h = h_start;
                    do{
                        if(fdapde::internals::point_in_2d_tri(n_new, old_coords, h->next()->node()->coords(), h->prev()->node()->coords()) ){
                            M_.update_node(n, old_coords, h->next()->node(), h->prev()->node() );
                            break;
                        }
                        h = h->prev()->twin();
                    }while(h!=h_start);  // n_new always falls inside the cavity if the movement is valid
                    
                    mesh_.flip();
                }
            }
        }
    }

    bool is_valid_movement_(node_t* n, const coords_t& new_pos) const {
        // check that no edges will be inverted by moving n to new_pos
        halfedge_t* h_start = n->halfedge();
        halfedge_t* h = h_start;
        do {
            coords_t A = h->next()->node()->coords();
            coords_t B = h->prev()->node()->coords();

            double area_after = fdapde::internals::signed_measure_2d_tri(A, B, new_pos);

            if (area_after <= 0) 
                return false;  // edge would be inverted or degenerate

            h = h->prev()->twin();
        } while (h != h_start);

        return true; // all edges valid
    }




};

}  // namespace fdapde

#endif // __ADAPTIVITY_H__








// CICLO COLLAPSE SUGLI EDGE CON OPZIONALE VALUTAZIONE DEL COSTO
                /*do {
                collapsed_in_pass = false; 
                for (auto it = mesh_.dcel_.halfedges_begin(); it != mesh_.dcel_.halfedges_end(); ++it) {
                    halfedge_t* e = &(*it);
                    if (e->id() >= e->twin()->id()) continue;  
                    std::unordered_set <cell_t*> modified={};
                    std::unordered_set <node_t*> deleted_nodes={};
                    if (M_.metric_edge_length(e) < M_.Lmin()) {
                        if(!e->is_segment()){
                            //auto node_cost = best_collapse_cost_(e, max_area);  
                            //if(node_cost.second!= std::numeric_limits<double>::max()){
                            //    deleted_nodes.insert(node_cost.first);
                            //    modified = collapse_edge_(e, node_cost.first);
                            //}
                            deleted_nodes.insert(e->node());
                            modified = collapse_edge_(e, e->node());
                        }
                        else {
                            if(!e->prev()->is_segment() && fdapde::internals::collinear(e->prev()->twin()->prev()->node()->coords(), e->node()->coords(), e->next()->node()->coords() )){
                                deleted_nodes.insert(e->node());
                                modified = collapse_segment_(e, e->node());
                            }
                            else if (!e->next()->is_segment()  && fdapde::internals::collinear(e->node()->coords(), e->twin()->node()->coords(), e->next()->twin()->next()->twin()->node()->coords() )){
                                    deleted_nodes.insert(e->twin()->node());
                                    modified = collapse_segment_(e, e->twin()->node());
                            }
                        }
                        if(!modified.empty()){
                            structural_change = true;
                            collapsed_in_pass = true;
                            ++it1;
                            M_.delete_nodes(deleted_nodes);
                            //refinenement_local_(min_angle, max_area, modified);
                            break; 
                        }
                    }
                }
            } while (collapsed_in_pass);*/




// CICLO VECCHIO CHE FA COLLAPSE/SPLIT TUTTO INSIEME E FUNZIONE CON LA LOGICA DELLA METRICA-STRATEGY
            /*while (it != end && iter < max_iter) {
            halfedge_t* e = &(*it);
            it = std::next(it);  
            if (e->id() >= e->twin()->id())    continue;  // each edge considered only once

            double lM = M_.metric_edge_length(e); 
            std::unordered_set<cell_t*> modified = {} ;  
            std::unordered_set<node_t*> deleted_nodes = {};

            if (std::isfinite(lM)) { 
                //else if (lM < M_.metPar().Lmin){ 
                //mesh_.dcel().export_to_json("Meshes/Delaunay/delaunay_output.json");
                if (lM < M_.Lmin()){  
                    if(!e->is_segment()){
                        auto node_cost = best_collapse_cost_(e, max_area);  // find best node to collapse between endpoint and midpoint
                        if(node_cost.second!= std::numeric_limits<double>::max()){
                            deleted_nodes.insert(node_cost.first);
                            modified = collapse_edge_(e, node_cost.first);
                        }
                    }
                    else {
                        if(!e->prev()->is_segment() && fdapde::internals::collinear(e->prev()->twin()->prev()->node()->coords(), e->node()->coords(), e->next()->node()->coords() )){
                            deleted_nodes.insert(e->node());
                            modified = collapse_segment_(e, e->node());
                        }
                        else if (!e->next()->is_segment()  && fdapde::internals::collinear(e->node()->coords(), e->twin()->node()->coords(), e->next()->twin()->next()->twin()->node()->coords() )){
                                deleted_nodes.insert(e->twin()->node());
                                modified = collapse_segment_(e, e->twin()->node());
                        }
                    }
                }
                //if (lM > M_.metPar().Lmax) {
                if (lM > M_.Lmax()) {
                    if(!e->is_segment())
                        modified = split_edge_(e);
                    else
                        modified = split_segment_(e);
                }
                bool did = !modified.empty();
                ++iter; 
                it  = mesh_.dcel_.halfedges_begin();   // restart without risking invalidations
                end = mesh_.dcel_.halfedges_end();

                /*if (did) {
                    std::cout << "Edge length in metric: " << lM << std::endl;
                    auto last_node = std::prev(mesh_.dcel_.nodes_end());
                    std::multimap<double, cell_t*> bad_triangles = {};
                    std::unordered_set<halfedge_t*> segments = {}, encroached = {};
                    std::unordered_set<node_t*> nodes_modified;  
                    for (auto* c : modified) {
                        if (!c || !c->halfedge()) continue;
                        nodes_modified.insert(c->halfedge()->node());
                        nodes_modified.insert(c->halfedge()->next()->node());
                        nodes_modified.insert(c->halfedge()->prev()->node());
                    }

                    refinenement_local_(min_angle, max_area, modified);
                    for(auto it= std::next(last_node); it != mesh_.dcel_.nodes_end(); ++it){
                        auto n = &(*it);
                        nodes_modified.insert(n); // refinement doesn't delete nodes, only adds new ones
                    }
                    
                    std::cout << "Adaptivity cycle iteration " << iter << std::endl;
                    /*M_.update_metric(nodes_modified, deleted_nodes); 
                    auto stats = M_.length_stats(tol);
                    double mean_diff = std::abs(stats.mean - prev_mean);
                    double frac_diff = std::abs(stats.frac_in_band - prev_frac);
                    bool lengths_stable =  (mean_diff / stats.mean < 0.1 * tol);
                    bool fraction_stable = (frac_diff < tol*1e-2);

                    if (lengths_stable && fraction_stable){
                        stable_count++;
                    }
                    else
                        stable_count = 0;

                    // convergence if stable for enough iterations
                    if (stable_count >= stable_iters_required) {
                        std::cout << "Adaptivity converged after " << iter << " iterations." << std::endl;
                        break;
                    }
                    // alternative stopping criterion if metric edge lengths are almost 1
                    bool no_more_splits = (stats.frac_in_band > 0.95);
                    if (no_more_splits && iter > 100) {
                        std::cout << "No more edges violating metric. Stopping after " << iter << " iterations." << std::endl;
                        break;
                    }
                    
                    prev_mean = stats.mean;
                    prev_frac = stats.frac_in_band;
                    
                    ++iter; 
                    it  = mesh_.dcel_.halfedges_begin();   // restart without risking invalidations
                    end = mesh_.dcel_.halfedges_end();
                }
            }     
        }*/

    
    // VECCHIO COLLAPSE COST
    /*std::pair<node_t*, double> best_collapse_cost_(halfedge_t* e, double max_area = 0.0) {
        node_t* u = e->node();
        node_t* v = e->twin()->node();
        std::vector<node_t*> candidates = {u, v, nullptr};  // nullptr stands for the midpoint (not insterted in the triangulation yet)

        double best_cost = std::numeric_limits<double>::max();
        node_t* best_node = nullptr;

        for (auto n: candidates) {
            auto result = simulate_collapse_(e, n);   
            if (!result.valid()) continue;  // collapse not possible
            //fdapde::Metric<LocalDim, EmbedDim, Strategy> M_temp(result.temp_dcel->dcel_, data_points_, max_area);
            std::unordered_set<cell_t*> elems_modified;
            for(auto it = result.temp_dcel->dcel_.cells_begin(); it != result.temp_dcel->dcel_.cells_end(); ++it)
                elems_modified.insert(const_cast<cell_t*>(&*it));

            double mean_length = 0.0;
            std::unordered_set<halfedge_t*> new_edges = {};
            for(auto cell: elems_modified){
                auto h_start = cell->halfedge();
                auto h = h_start;
                do{
                    if(new_edges.find(h) == new_edges.end() || new_edges.find(h->twin()) == new_edges.end()){
                        new_edges.insert(h);
                        new_edges.insert(h->twin());
                        mean_length += M_.metric_edge_length(h);
                    }
                    h=h->next();
                }while(h!=h_start);
            }
            mean_length /= (2*new_edges.size());   //PERCHE HO MOLTIPLICATO PER DUE???

            double difference = std::abs(1 - mean_length);
            if(difference<best_cost){
                best_cost = difference;
                best_node = n;
            }
        }
        return {best_node,best_cost};
    }*/




// VECCHIO COSTRUTTORE E CICLO ADATTIVO CHE FACEVA SOLO COLLAPSE FINO A N_POINTS CON COST-OBJECT VECCHIO
    /*Adaptivity(delaunay_t& mesh, const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& data_points, const int n_points, std::vector<double> costs_weights = {}): 
        mesh_(mesh), data_points_(data_points), cost_obj_(mesh_.dcel_, data_points_) {
        adaptivity_cycle(n_points);
    }
    void adaptivity_cycle(const int n_points){
        while(mesh_.dcel_.n_nodes() > n_points) {
            // find worst edge to collapse
            double best_cost = std::numeric_limits<double>::max();
            halfedge_t* best_edge = nullptr;
            node_t* best_node = nullptr;
            for (auto it = mesh_.dcel_.halfedges_begin(); it != mesh_.dcel_.halfedges_end(); ++it) {
                halfedge_t* e = &(*it);
                if (e->id() < e->twin()->id()) {  // each edge considered only once
                    auto [n, cost] = best_collapse_cost_(e);
                    if (cost < best_cost) {
                        best_cost = cost;
                        best_edge = e;
                        best_node = n;
                    }
                }
            }
            if (best_edge) {
                auto modified_cells = collapse_edge_(best_edge, best_node);
                cost_obj_.update(modified_cells);  
            }
            else break;
        }
    }*/



// VECCHIO BEST COLLAPSE COST CHEUSAVA IL COST_OBJ VECCHIO
    // find the best node to delete between the two edge extremes or the midpoint (not inserted yet)
    // it is the one that minimizes the cost
    /*std::pair<node_t*, double> best_collapse_cost_(halfedge_t* e) {
        node_t* u = e->node();
        node_t* v = e->twin()->node();
        std::vector<node_t*> candidates = {u, v, nullptr};  // nullptr stands for the midpoint (not insterted in the triangulation yet)

        double best_cost = std::numeric_limits<double>::max();
        node_t* best_node = nullptr;

        for (auto n: candidates) {
            auto result = simulate_collapse_(e, n);   
            if (!result.valid()) continue;  // collapse not possible
            std::unordered_set<cell_t*> elems_modified;
            for(auto it = result.temp_dcel->dcel_.cells_begin(); it != result.temp_dcel->dcel_.cells_end(); ++it)
                elems_modified.insert(const_cast<cell_t*>(&*it));
            double cost = cost_obj_.get_cost({elems_modified}, result.cavity_info);
            if (cost < best_cost) {
                best_cost = cost;
                best_node = n;
            }
        }
        return {best_node,best_cost};
    }*/


    // ADAPTIVITY COST per plottare percorsi geodesici
        // Optional: save geodesic distance field for debugging / plotting
        //std::vector<std::pair<coords_t, double>> dist_field;
        //dist_field.reserve(1024);

            // (optional) store finalized distances for visualization  (nel ciclo)
            //dist_field.push_back({v->coords(), dist});

        /*// --- Reconstruct TRUE shortest path v0 -> v_k ---
        std::vector<node_t*> path_nodes;
        path_nodes.reserve(256);
        node_t* cur = v_k;
        while (cur)
        {
            path_nodes.push_back(cur);
            if (cur == v0) break;
            auto itp = parent.find(cur);
            if (itp == parent.end()) break; // safety
            cur = itp->second;
        }
        std::reverse(path_nodes.begin(), path_nodes.end());

        // --- Export TRUE path as CSV (no spaghetti) ---
        std::string fname = "Meshes/Delaunay/sequenze/path_" +std::to_string(v0->id()) + "_to_" +std::to_string(v_k->id()) + ".csv";
        std::ofstream out(fname);
        if (out.is_open()){
            out << "x,y\n";
            for (auto* n : path_nodes){
                auto c = n->coords();
                out << c[0] << "," << c[1] << "\n";
            }
            out.close();
        }*/


    
    // ANISOTROPA --> prova
    /*Eigen::Matrix<double,2,2> MA = node_metrics_.at(A);
    Eigen::Matrix<double,2,2> MB = node_metrics_.at(B);
    coords_t u = AB/std::sqrt(AB.squaredNorm());
    // h at endpoints, projected along edge direction
    double lambdaA = (u * MA * u.transpose())(0,0);
    double lambdaB = (u * MB * u.transpose())(0,0);
    double hA = 1.0 / sqrt(lambdaA);
    double hB = 1.0 / sqrt(lambdaB);
    // geometric interpolation
    double a = hB / hA;
    // avoid log singularity
    if (fabs(a - 1.0) < 1e-12)
        return std::sqrt(AB.squaredNorm()) / hA;
    double la = std::sqrt(AB.squaredNorm()) / hA;
    double L = la * (a - 1.0) / (a*log(a));
    return L; */   
            