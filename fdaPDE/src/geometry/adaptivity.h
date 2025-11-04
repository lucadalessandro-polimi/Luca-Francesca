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


template <int LocalDim, int EmbedDim, AdaptiveStrategy Strategy>
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

    using storage_t = typename StrategyData<embed_dim, Strategy>::storage_t;


    /*Adaptivity(delaunay_t& mesh, const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& data_points, const int n_points, std::vector<double> costs_weights = {}): 
        mesh_(mesh), data_points_(data_points), cost_obj_(mesh_.dcel_, data_points_) {
        adaptivity_cycle(n_points);
    }*/

    Adaptivity(delaunay_t& mesh, const storage_t& data_points, const double min_angle, const double max_area, const double tol = 0.1): 
        mesh_(mesh), data_points_(data_points), M_(mesh_.dcel_, data_points_, max_area) {
        //adaptivity_cycle_two(min_angle, max_area, tol);
    }


    /*void adaptivity_cycle(const int n_points){
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

    // method to adapt the mesh according to a metric M, by splitting and collapsing edges
    int adaptivity_cycle_two(const double min_angle, const double max_area, const double tol){
        //fdapde::Metric M(mesh_.dcel_, data_points_, max_area);
        int iter = 0;
        auto it  = mesh_.dcel_.halfedges_begin();
        auto end = mesh_.dcel_.halfedges_end();
        
        double prev_mean = 1.0;
        double prev_frac = 0.0;
        int stable_count = 0;
        const int stable_iters_required = 10;
        int max_iter = 5000;

        while (it != end && iter < max_iter) {
            halfedge_t* e = &(*it);
            it = std::next(it);  
            if (e->id() >= e->twin()->id())    continue;  // each edge considered only once

            double lM = M_.metric_edge_length(e); 
            std::unordered_set<cell_t*> modified = {} ;  
            std::unordered_set<node_t*> deleted_nodes = {};

            if (std::isfinite(lM)) {
                if (lM > M_.metPar().Lmax) {
                    if(!e->is_segment())
                        modified = split_edge_(e);
                    else
                        modified = split_segment_(e);
                } 
                else if (lM < M_.metPar().Lmin){  
                    if(!e->is_segment()){
                        //std::cout << "Prima di collapse di " << e->id() << std::endl;
                        auto node_cost = best_collapse_cost_(e, max_area);  // find best node to collapse between endpoint and midpoint
                        //std::cout << "Best collapse cost for edge " << e->id() << std::endl;
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
                bool did = !modified.empty();

                if (did) {

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

                    for(auto it= std::next(last_node); it != mesh_.dcel_.nodes_end(); ++it){
                        auto n = &(*it);
                        nodes_modified.insert(n); // refinement doesn't delete nodes, only adds new ones
                    }
                    
                    std::cout << "Adaptivity cycle iteration " << iter << std::endl;
                    //mesh_.dcel_.export_to_json("Meshes/Delaunay/delaunay_output.json");
                    M_.update_metric(nodes_modified, deleted_nodes); 
                    auto stats = M_.length_stats(tol);
                    double mean_diff = std::abs(stats.mean - prev_mean);
                    double frac_diff = std::abs(stats.frac_in_band - prev_frac);
                    //std::cout << "fraction: " << stats.frac_in_band << std::endl;
                    //std::cout << "frac diff: " << frac_diff << std::endl;
                    //std::cout << "mean diff: " << mean_diff << std::endl;
                    bool lengths_stable =  (mean_diff / stats.mean < 0.1 * tol);
                    bool fraction_stable = (frac_diff < tol*1e-2);

                    if (lengths_stable && fraction_stable){
                        stable_count++;
                        //std::cout << "stable: " << stable_count << std::endl;
                    }
                    else
                        stable_count = 0;

                    //if (std::abs(stats.mean - 1) > tol) {
                        //M_.rescale_c(stats.mean);
                        //std::cout << "Rescale c at iter " << iter << std::endl;
                    //}

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
        }

        mesh_.check_quality_(min_angle, max_area);
        return iter;

    }

    struct CollapseSimulate {
        std::optional<delaunay_t> temp_dcel;
        std::unordered_map<int, std::pair<bool,int>> cavity_info;          // <halfedge on boundary, vertex degree>
        bool valid() const { return !cavity_info.empty(); }
    };

    dcel_t&  dcel()  { return mesh_.dcel_; }

    fdapde::Metric<local_dim, embed_dim, Strategy>& metric() { return M_; }



   private:
    delaunay_t& mesh_;
    const storage_t data_points_; // data points
    fdapde::Metric<local_dim, embed_dim, Strategy> M_;
    //DataEquiCost<local_dim, embed_dim> cost_obj_;  // object to compute costs   PER ORA SOLO DATAEQUICOST

    bool modifiable_edge_(halfedge_t* e) const {
        double lM = M_.metric_edge_length(e);
        if (std::isfinite(lM) && (lM > M_.metPar().Lmax || lM < M_.metPar().Lmin) )
            return true;
        return false;
    }

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

    std::pair<node_t*, double> best_collapse_cost_(halfedge_t* e, double max_area = 0.0) {
        node_t* u = e->node();
        node_t* v = e->twin()->node();
        std::vector<node_t*> candidates = {u, v};  // nullptr stands for the midpoint (not insterted in the triangulation yet)

        double best_cost = std::numeric_limits<double>::max();
        node_t* best_node = nullptr;

        for (auto n: candidates) {
            auto result = simulate_collapse_(e, n);   
            if (!result.valid()) continue;  // collapse not possible
            fdapde::Metric<LocalDim, EmbedDim, Strategy> M_temp(result.temp_dcel->dcel_, data_points_, max_area);
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
                        mean_length += M_temp.metric_edge_length(h);
                    }
                    h=h->next();
                }while(h!=h_start);
            }
            mean_length /= (2*new_edges.size());

            double difference = std::abs(1 - mean_length);
            if(difference<best_cost){
                best_cost = difference;
                best_node = n;
            }
        }
        return {best_node,best_cost};
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
            //std::cout << cavity_coords << std::endl;
            result.temp_dcel.emplace({cavity_coords});  // Delaunay cavity
            // make the cavity information coherent with the actual mesh information (create a mapping)
            node_t* start_node = result.temp_dcel->dcel_.find_node(cavity_coords.row(0));
            halfedge_t* h = start_node->halfedge()->twin();  // external boundary (its cell is null)
            for(auto cavity_edge: polygonal_cavity){
                h->twin()->set_id(cavity_edge->id());
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
            halfedge_t* h = start_node->halfedge()->twin();  // external boundary (its cell is null)
            for(auto cavity_edge: polygonal_cavity){
                h->twin()->set_id(cavity_edge->id());
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
        coords_t split_pt = 0.5 * (A->coords() + B->coords());
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
        coords_t split_pt = 0.5 * (A->coords() + B->coords());
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
        return modified_cells;
    }



    /*void smooth_node_(node_t* n, const coords_t new_coords){
        n->coords() = new_coords;

        // local flip
        // ...
    }*/

    

    // only interior vertices that do not lie on segments or boundaries can be deleted
    // SERVE segment_ per node_t ??
    /*void delete_vertex_(node_t* n){
        if(!n->halfedge())
            return;
        if(n->halfedge()->is_segment())   // NON SONO SICURA SIA SUFFICIENTE
            return;                       // ECCEZIONE: SE NODO è A METà DI UN SEGMENTO (PUNTI PRIMA E DOPO SONO COLLINEARI), SI PUO FONDERE SEGMENTO -->2 CAVITà
        
        // identify the polygonal cavity and the edges to remove
        std::list<halfedge_t*> polygonal_cavity;
        std::list<halfedge_t*> edges_to_remove;
        
        halfedge_t* h_start = n->halfedge();
        halfedge_t* h = h_start;
        do{
            edges_to_remove.push_back(h);
            polygonal_cavity.push_back(h->next());
            h = h->prev()->twin();
        }while(h!=h_start);
        
        // remove edges incident to n
        for(auto e: edges_to_remove){
            mesh_.dcel_.remove_edge(e);
        }
        mesh_.dcel_.remove_node(n);     
        // triangulate cavity using polygon.h triangulation algorithm
        if(polygonal_cavity.size()>3){
            Eigen::Matrix<double, Eigen::Dynamic, 2> cavity_coords(polygonal_cavity.size(),2);
            int i=0;
            for(auto h: polygonal_cavity){
                cavity_coords.row(i) = h->node()->coords();
                ++i;
            }
            polygon_t polygon(cavity_coords);
            auto triangulation = polygon.triangulation();
            mesh_.dcel_.from_triangulation(triangulation, std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>>{});
        }
        // local flip
        // n-3 internal new edges added inside the cavity
        // n-2 new triangles created
        std::unordered_set<halfedge_t*> new_edges;
        int cont=0;
        auto rb = std::make_reverse_iterator(mesh_.dcel_.halfedges_end());
        auto re = std::make_reverse_iterator(mesh_.dcel_.halfedges_begin());
        for (auto it = rb; it != re && cont< polygonal_cavity.size()-3 ; ++it,++it) { //halfedge and its twin
            new_edges.insert(&(*it));
            cont++;
        }
        for(auto e: polygonal_cavity)
            if(!e->is_segment()){
                new_edges.insert(e);
            }

        mesh_.flip(new_edges);
        
    }*/

    



};

}  // namespace fdapde

#endif // __ADAPTIVITY_H__

