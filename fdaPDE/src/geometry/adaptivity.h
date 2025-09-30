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

    using Indicator = std::function<double(const cell_t&)>;
    using RefineLambda = std::function<void(delaunay_t& mesh)>;
    // o  using RefineLambda = std::function<std::optional<coords_t>>(const cell_t&)>;  (rstituisce punto in piu)



    Adaptivity(delaunay_t& mesh, RefineLambda adaptive_strategy) : mesh_(mesh), adaptive_strategy_(adaptive_strategy) {
        // adaptive cycle ..........
        // ...............
    }

    Adaptivity(delaunay_t& mesh, const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& data_points, const int n_points, std::vector<double> costs_weights = {}): 
        mesh_(mesh), data_points_(data_points), cost_obj_(mesh_.dcel_, data_points_) {
        adaptivity_cycle(n_points);
        std::cout << "Final number of nodes: " << mesh_.dcel_.n_nodes() << std::endl;
    }

    void adaptivity_cycle(const int n_points){
        while(mesh_.dcel_.n_nodes() > n_points) {
            std::cout << "Adaptive cycle: " << mesh_.dcel_.n_nodes() << " nodes" << std::endl;
            mesh_.dcel_.export_to_json("Meshes/Delaunay/delaunay_output.json");
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
                std::cout << "Best edge to collapse: " << best_edge->id() << " with cost: " << best_cost << std::endl;
                if(!best_node) std::cout << "midpoint collapse" << std::endl;
                auto modified_cells = collapse_edge_(best_edge, best_node);
                cost_obj_.update(modified_cells);  
                std::cout << "Fatto il collapse" << std::endl;
            }
            else break;

        }
    } 

    struct CollapseSimulate {
        std::optional<delaunay_t> temp_dcel;
        std::unordered_map<int, std::pair<bool,int>> cavity_info;          // <halfedge on boundary, vertex degree>
        bool valid() const { return !cavity_info.empty(); }
    };

    dcel_t&  dcel()  { return mesh_.dcel_; }



   private:
    delaunay_t& mesh_;
    RefineLambda adaptive_strategy_;
    const Eigen::Matrix<double, Eigen::Dynamic, embed_dim> data_points_; // data points
    DataEquiCost<local_dim, embed_dim> cost_obj_;  // object to compute costs   PER ORA SOLO DATAEQUICOST


    // find the best node to delete between the two edge extremes or the midpoint (not inserted yet)
    // it is the one that minimizes the cost
    std::pair<node_t*, double> best_collapse_cost_(halfedge_t* e) {
        node_t* u = e->node();
        node_t* v = e->twin()->node();
        std::vector<node_t*> candidates = {u, v, nullptr};  // nullptr stands for the midpoint (not insterted in the triangulation yet)

        double best_cost = std::numeric_limits<double>::max();
        node_t* best_node = nullptr;

        for (auto n: candidates) {
            std::cout << "inizio simulate collapse di " << e->id() << std::endl;
            auto result = simulate_collapse_(e, n);   
            std::cout << "finisco simulate collapse di " << e->id() << std::endl;
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
    }


    // collapse edge e into node; if node==nullptr, insert the midpoint of the edge
    std::unordered_set<cell_t*> collapse_edge_(halfedge_t* e, node_t* node) {
        std::unordered_set<cell_t*> elems_modified;

        if (e->is_segment() || (node && node->on_boundary()) || (!node && (e->node()->on_boundary() || e->twin()->node()->on_boundary()))) 
            return elems_modified;

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

        } else { // midpoint case
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
        for (auto he : edges_to_remove)
            mesh_.dcel_.remove_edge(he);

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
            mesh_.dcel_.export_to_json("Meshes/Delaunay/delaunay_output.json");
        }

        // local flip
        // n cavity edges, n-3 internal new edges added inside the cavity, n-2 new triangles created
        /*std::unordered_set<halfedge_t*> new_edges;
        int cont = 0;
        auto rb = std::make_reverse_iterator(mesh_.dcel_.halfedges_end());
        auto re = std::make_reverse_iterator(mesh_.dcel_.halfedges_begin());
        for (auto it = rb; it != re && cont < polygonal_cavity.size() - 3; ++it, ++it) {
            new_edges.insert(&(*it));
            cont++;
        }
        mesh_.flip(new_edges); */  // ensure the Delaunay property   // SERVE IN TEORIA NO....

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
        // cavity_coords are in counterclockwise order, so each node id and edge id corresponds to the corresponding index in cavity_info
        // NO NEI CASI IN CUI I NODI SONO COLLINEARI E VENGONO AGGIUNTI DOPO ...
        std::unordered_set<cell_t*> elems_modified;
        if(node){
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
                h = h->prev(); //to go in counterclockwise order
            } 
        }        

        result.cavity_info = cavity_info;
        return result;
    }

    // count the number of triangles incident to a vertex (same as the number of edges radiating from said vertex)
    int count_triangles_from_vertex_(node_t* v) const {
        int count=0;
        if(!v->halfedge()->cell()) v->set_halfedge(v->halfedge()->twin());  // ensure starting from an internal halfedge
        halfedge_t* h_start = v->halfedge();
        halfedge_t* h = h_start;
        do {
            count++;
            h = h->prev()->twin();
            if(!h->cell()){ // on boundary
               if(h->twin()==h_start) break; 
               h = h->prev()->twin();  // skip the boundary (no cell in between halfedges)
            }
        } while (h != h_start);
        return count;
    }
    


    // LOCAL FLIP in ogni operazione o nel MASTER ALGORITHM
    // +++ EDGE COLLAPSE

    /*void smooth_node_(node_t* n, const coords_t new_coords){
        n->coords() = new_coords;

        // local flip
        // ...
    }

    void split_edge_(halfedge_t* e){

        node_t* a = e->node();
        node_t* b = e->twin()->node();
        node_t* c = e->prev()->node();

        // coordinates of the point to insert to split e
        coords_t split_pt; 

        // CRITERION TO SPLIT EDGE
        // ......
        if (e->next()->is_segment() &&
            (fdapde::internals::segment_length(a->coords(), b->coords()) - fdapde::internals::segment_length(b->coords(), c->coords()) ) / fdapde::internals::segment_length(a->coords(), b->coords()) > 0.1 &&
            fdapde::internals::angle_between(a->coords(), b->coords(), c->coords()) <=45) {
            // project c onto the line ab s.t. b_split = bc
            coords_t split_pt_ref = c->coords();
            double r = (split_pt_ref - b->coords()).norm();
            coords_t ab = a->coords() - b->coords();
            double L = ab.norm();
            if (r > L) r = 0.5 * L;
            double t = r / L;
            split_pt = b->coords() + t * ab;
        } else {  // split e in the middle
            split_pt = 0.5 * (a->coords() + b->coords());
        }

        // insert the split point into the dcel_
        node_t* m = dcel_.insert_node(node_t(dcel_.n_nodes(), e->on_boundary(), split_pt));

        halfedge_t* prev = e->prev(); //ca
        halfedge_t* twin_prev = e->twin()->prev();   
        halfedge_t* next = e->next();  //bc

        // create new connections between am, mb, bc before removing e
        halfedge_t* h1 = dcel_.emplace_halfedge(m, true);
        h1->set_cell(next->cell());
        dcel_.insert_edge(next, h1); //mb
        h1->twin()->set_segment(true);
        halfedge_t* h2 = dcel_.insert_edge(prev->next(), h1);  //am
        h2->set_segment(true);
        h2->twin()->set_segment(true);
        dcel_.insert_edge(prev, h1);  //mc
        dcel_.remove_edge(e);

        // local flip
        // ...
    }*/

    // only interior vertices that do not lie on segments or boundaries can be deleted
    // SERVE segment_ per node_t ??
    void delete_vertex_(node_t* n){
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
        
    }

    



};

}  // namespace fdapde

#endif // __ADAPTIVITY_H__

