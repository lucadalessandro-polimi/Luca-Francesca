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

#ifndef __FDAPDE_DCEL_H__
#define __FDAPDE_DCEL_H__

#include "header_check.h"
#include <chrono>
using namespace std::chrono;

namespace fdapde {

// implementation of the Double Connected Edge List data structure (also known as DCEL or half-edge)
template <int LocalDim, int EmbedDim> class DCEL {
   public:
    static constexpr int local_dim = LocalDim;
    static constexpr int embed_dim = EmbedDim;
    // forward decl
    struct node_t;
    struct halfedge_t;
    struct cell_t;
    using coords_t = Eigen::Matrix<double, embed_dim, 1>;
    // internal data structures
    struct node_t {
        private:
         int id_;                  // global node index
         halfedge_t* halfedge_;    // any edge having this node as its origin
         bool boundary_;           // asserted true if node is on boundary
         coords_t coords_;
         // code needed for conflict graph algorithm
         cell_t* conflicting_triangle_=nullptr;
 
         public:
 
         node_t() : coords_(), halfedge_(nullptr), boundary_(false) { }
 
         template <typename CoordsType>
             requires(internals::is_eigen_dense_xpr_v<CoordsType>)
         node_t(int id, halfedge_t* halfedge, bool boundary, const CoordsType& coords) :
             id_(id), halfedge_(halfedge), boundary_(boundary), coords_() {
             fdapde_assert(
               (coords.rows() == 1 && coords.cols() == embed_dim) || (coords.rows() == embed_dim && coords.cols() == 1));
             if (coords.rows() == 1) {
                 coords_ = coords.transpose();
             } else {
                 coords_ = coords;
             }
         }
 
         template <typename CoordsType>
             requires(internals::is_eigen_dense_xpr_v<CoordsType>)
         node_t(int id, bool boundary, const CoordsType& coords) : node_t(id, nullptr, boundary, coords) { }
        
         template <typename... CoordsType>
             requires(std::is_floating_point_v<CoordsType> && ...) && (sizeof...(CoordsType) == embed_dim)
         node_t(int id, halfedge_t* halfedge, bool boundary, CoordsType&&... coords) :
             id_(id), halfedge_(halfedge), boundary_(boundary), coords_(coords...) { }
        
             template <typename... CoordsType>
             requires(std::is_floating_point_v<CoordsType> && ...) && (sizeof...(CoordsType) == embed_dim)
         node_t(int id, bool boundary, CoordsType&&... coords) :
             node_t(id, nullptr, boundary, coords...) { }
 
 
         // observers and modifiers
         const Eigen::Matrix<double, embed_dim, 1>& coords() const { return coords_; }
         halfedge_t* halfedge() const { return halfedge_; }
         void set_halfedge(halfedge_t* halfedge) { halfedge_ = halfedge; }
         int id() const { return id_; }
         bool on_boundary() const { return boundary_; }
         void set_boundary(bool boundary) { boundary_ = boundary; }
         node_t* next() const { return halfedge_->next()->node(); }
         node_t* prev() const { return halfedge_->prev()->node(); }
 
         // code for conflict graph
         void set_conflict(cell_t* triangle) { conflicting_triangle_ = triangle; }
         cell_t* conflict() const { return conflicting_triangle_; }
         void remove_conflict() { conflicting_triangle_ = nullptr; }
 

    };
    struct halfedge_t {
       private:
        int id_;   // global halfedge index
        halfedge_t *prev_, *next_, *twin_;
        node_t* node_;
        cell_t* cell_;   // cell to which this halfedge belongs to
        std::list<halfedge_t>::iterator it_;   // iterator to the halfedge in the list
        bool subsegment_ ;  // true if the halfedge needs to be mantained in the mesh
       public:
        halfedge_t() : node_(nullptr), prev_(nullptr), next_(nullptr), twin_(nullptr), subsegment_(false) { }
        halfedge_t(int id, halfedge_t* prev, halfedge_t* next, halfedge_t* twin, node_t* node, bool sub=false) :
            id_(id), prev_(prev), next_(next), twin_(twin), node_(node), subsegment_(sub) { }
        // no twin constructors
        halfedge_t(int id, halfedge_t* prev, halfedge_t* next, node_t* node, bool sub=false) :
            halfedge_t(id, prev, next, nullptr, node, sub) { }
        // minimal constructor
        halfedge_t(int id, node_t* node, bool sub= false) : halfedge_t(id, nullptr, nullptr, nullptr, node, sub) { }

        // observers
        halfedge_t* prev() const { return prev_; }
        halfedge_t* next() const { return next_; }
        halfedge_t* twin() const { return twin_; }
        node_t* node() const { return node_; }
        cell_t* cell() const { return cell_; }
        int id() const { return id_; }
        bool on_boundary() const { return (node_->on_boundary() && twin_->node()->on_boundary() && (cell()==nullptr || twin()->cell()==nullptr)); }
        std::list<halfedge_t>::iterator it() const { return it_; }
        bool is_subsegment() const { return subsegment_; }
        // modifiers
        void set_prev(halfedge_t* prev) { prev_ = prev; }
        void set_next(halfedge_t* next) { next_ = next; }
        void set_twin(halfedge_t* twin) { twin_ = twin; }
        void set_node(node_t* node) { node_ = node; }
        void set_cell(cell_t* cell) { cell_ = cell; }
        void set_id(int id) {id_=id;}
        void set_it(std::list<halfedge_t>::iterator it) { it_ = it; }
        void set_subsegment(bool sub) { subsegment_ = sub; }

        // iterator (follows the chain of directed edges until no next valid edge or this edge is found)
        struct circulator {
            using value_type = halfedge_t;
            using pointer = std::add_pointer_t<value_type>;
            using reference = std::add_lvalue_reference_t<value_type>;
            using size_type = std::size_t;
            using difference_type = std::ptrdiff_t;
            using iterator_category = std::forward_iterator_tag;

            circulator(halfedge_t* halfedge) :
                halfedge_(halfedge), end_(halfedge == nullptr ? nullptr : halfedge->prev()) { }
            circulator& operator++() {
                if (last_) { [[unlikely]]
                    end_ = nullptr;
                } else {
                    halfedge_ = halfedge_->next();
                    if (halfedge_ == end_) { last_ = true; }   // implement cyclic structure
                }
                return *this;
            }
            // access
            pointer operator->() { return halfedge_; }
            const pointer operator->() const { return halfedge_; }
            reference operator*() { return *halfedge_; }
            const reference operator*() const { return *halfedge_; }
            operator bool() const { return end_ == nullptr; }
            // comparison
            friend bool operator==(const circulator& lhs, const circulator& rhs) { return lhs.end_ == rhs.end_; }
            friend bool operator!=(const circulator& lhs, const circulator& rhs) { return lhs.end_ != rhs.end_; }
           private:
            bool last_ = false;
            pointer halfedge_, end_;
        };
    };
    struct cell_t {
        cell_t() : h_(nullptr) { }
        cell_t(int id) : id_(id), h_(nullptr) { }
        cell_t(int id, halfedge_t* h) : id_(id), h_(h) { }
        std::list<cell_t>::iterator it_;   // iterator to the cell in the list
        // observers
        halfedge_t* halfedge() const { return h_; }
        int id() const { return id_; }
        std::list<cell_t>::iterator it() const { return it_; }
        // modifiers
        void set_halfedge(halfedge_t* h) { h_ = h; }
        void set_id(int id) {id_=id;}
        void set_it(std::list<cell_t>::iterator it) { it_ = it; }
          
         bool operator==(const cell_t& other) const {
           return id_ == other.id_;  
         }
 
         // code for conflict graph algorithm
         void add_conflict(node_t* point) { conflicting_points_.push_back(point); }
         std::vector<node_t*>& conflicting_points() const{ return conflicting_points_; }
         std::vector<node_t*>& conflicting_points() { return conflicting_points_; }
         void clear_conflicts() { conflicting_points_.clear(); }
      
 
        private:
         int id_;
         halfedge_t* h_;
         std::vector<node_t*> conflicting_points_;  // code needed for conflict graph algorithm
     };


    using halfedge_iterator = std::list<halfedge_t>::iterator;
    using node_iterator = std::list<node_t>::iterator;
    using cell_iterator = std::list<cell_t>::iterator;

    // constructors
    DCEL() : nodes_(), halfedges_(), n_nodes_(0), n_halfedges_(0), n_cells_(0) { }
    // constructs a closed loop structure linking nodes one after the other
    /*static DCEL<local_dim, embed_dim> make_polygon(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& nodes_entry) {
        fdapde_assert(nodes_entry.cols() == embed_dim);
        Eigen::Matrix<double, Dynamic, Dynamic> nodes=nodes_entry;
        if (!internals::are_2d_counterclockwise_sorted(nodes_entry)) {
            int n_nodes = nodes_entry.rows();
            for (int i = 0; i < n_nodes; ++i)
                nodes.row(i) = nodes_entry.row(n_nodes - 1 - i);
        }

        int n_nodes = nodes.rows();
        int n_halfedges = 2 * (n_nodes);
        DCEL<local_dim, embed_dim> dcel;
        // create polygon cell
        dcel.cells_.push_back(cell_t(0));
        cell_t* c = std::addressof(dcel.cells_.back());
        c->set_it(std::prev(dcel.cells_.end()));
        dcel.n_cells_ = 1;
        // push nodes
        for (int i = 0; i < n_nodes; ++i) {
            node_t* n = dcel.insert_node(node_t(i, true, nodes.row(i)));
            halfedge_t* h = dcel.emplace_halfedge_(n);
            n->set_halfedge(h);
	        h->set_cell(c);
        }
        c->set_halfedge(dcel.nodes_begin()->halfedge());
        // create twin edges
        for (auto it = dcel.nodes_begin(); it != dcel.nodes_end(); ++it) {
            node_t* n1 = std::addressof(*it);
            node_t* n2 = std::addressof(*((it->id() == n_nodes - 1) ? dcel.nodes_begin() : std::next(it, 1)));
            halfedge_t* h1 = n1->halfedge();
            halfedge_t* h2 = dcel.emplace_halfedge_(n2);   // push twin edge
            h2->set_cell(nullptr);
            h2->set_twin(h1);
            h1->set_twin(h2);
        }
        // finalize next-prev pointer pairs
        for (auto it = dcel.nodes_begin(); it != dcel.nodes_end(); ++it) {
            halfedge_t* h1 = it->halfedge();
            halfedge_t* h2 = ((it->id() == n_nodes - 1) ? dcel.nodes_begin() : std::next(it, 1))->halfedge();
            h1->set_next(h2);
            h2->set_prev(h1);
	        h1->twin()->set_prev(h2->twin());
            h2->twin()->set_next(h1->twin());}
        return dcel;
    }*/
    

    // overloading of make_polygon to also add holes
    static DCEL<local_dim, embed_dim> make_polygon(
        const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& boundary,
        const std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>>& holes={}) {
    
        fdapde_assert(boundary_entry.cols() == embed_dim);
        DCEL<local_dim, embed_dim> dcel;

        // check if the boundaries are counterclockwise sorted
        /*Eigen::Matrix<double, Eigen::Dynamic, embed_dim> boundary=boundary_entry;
        if (!internals::are_2d_counterclockwise_sorted(boundary_entry)) {
            int n_nodes = boundary_entry.rows();
            for (int i = 0; i < n_nodes; ++i)
                boundary.row(i) = boundary_entry.row(n_nodes - 1 - i);
        }
        std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>> holes;
        for (const auto& hole : holes_entry) {
            if (internals::are_2d_counterclockwise_sorted(hole)) {
                Eigen::Matrix<double, Eigen::Dynamic, embed_dim> reversed_hole(hole.rows(), embed_dim);
                for (int i = 0; i < hole.rows(); ++i)
                    reversed_hole.row(i) = hole.row(hole.rows() - 1 - i);
                holes.push_back(reversed_hole);
            } else {
                holes.push_back(hole);
            }
        }*/
    
        // external boundary
        int n_nodes = boundary.rows();
        dcel.cells_.push_back(cell_t(0)); 
        cell_t* c = std::addressof(dcel.cells_.back());
        c->set_it(std::prev(dcel.cells_.end()));
        dcel.n_cells_ = 1;
        // nodes and halfedges for external boundary
        for (int i = 0; i < n_nodes; ++i) {
            node_t* n = dcel.insert_node(node_t(i, true, boundary.row(i)));
            halfedge_t* h = dcel.emplace_halfedge_(n, true);
            n->set_halfedge(h);
            h->set_cell(c);
        }
        c->set_halfedge(dcel.nodes_begin()->halfedge());
        // twin halfedges for external boundary
        for (auto it = dcel.nodes_begin(); it != dcel.nodes_end(); ++it) {
            node_t* n1 = std::addressof(*it);
            node_t* n2 = std::addressof(*((it->id() == n_nodes - 1) ? dcel.nodes_begin() : std::next(it, 1)));
            halfedge_t* h1 = n1->halfedge();
            halfedge_t* h2 = dcel.emplace_halfedge_(n2, true); // Twin edge
            h2->set_cell(nullptr);
            h2->set_twin(h1);
            h1->set_twin(h2);
        }
        // next and prev for external boundary
        for (auto it = dcel.nodes_begin(); it != dcel.nodes_end(); ++it) {
            halfedge_t* h1 = it->halfedge();
            halfedge_t* h2 = ((it->id() == n_nodes - 1) ? dcel.nodes_begin() : std::next(it, 1))->halfedge();
            h1->set_next(h2);
            h2->set_prev(h1);
            h1->twin()->set_prev(h2->twin());
            h2->twin()->set_next(h1->twin());
        }

        
        /*int node_offset_bd= n_nodes;
        auto boundaries_vector = dcel.extract_internal_paths_with_connections_(boundaries);
        //adding internal segments (fixed regions)
        for(auto points_matrix: boundaries_vector){
            int n_nodes= points_matrix.rows();
            for(auto row: points_matrix.rowwise()){
                std::cout << "Boundary point: " << row.transpose() << std::endl;
            }
            for(int i=1; i< n_nodes-1; ++i){
                coords_t co= points_matrix.row(i).transpose();
                if(dcel.find_node(co)) continue;
                node_t* n1 = dcel.insert_node(node_t(dcel.n_nodes(), false, points_matrix.row(i)));
                halfedge_t* h1 = dcel.emplace_halfedge_(n1, true);
                n1->set_halfedge(h1);
                h1->set_cell(c);
                std::cout << "node: " << n1->id() << " coords: " << n1->coords().transpose() << std::endl;
            } 
            for (int i = 0; i < n_nodes-1; ++i) {
                coords_t co1= points_matrix.row(i).transpose();
                coords_t co2= points_matrix.row(i+1).transpose();
                node_t* n1= dcel.find_node(co1);
                node_t* n2= dcel.find_node(co2);
                std::cout << "n1: " << n1->id() << " n2: " << n2->id() << std::endl;
                if(!dcel.find_halfedge_between(co1,co2)){
                    halfedge_t* h1 = n1->halfedge();
                    halfedge_t* h2 = n2->halfedge();
                    dcel.insert_edge(h1, h2, true);
                    std::cout << "QUI" << std::endl;
                }
            }
            node_offset_bd += n_nodes -2;
        }*/

    
        // adding holes
        int node_offset = dcel.n_nodes(); 
        int hole_index = 1;
    
        for (const auto& hole : holes) {
            int hole_nodes = hole.rows();
            // nodes and halfedges for the hole
            for (int i = 0; i < hole_nodes; ++i) {
                node_t* n = dcel.insert_node(node_t(node_offset + i,  true, hole.row(i)));
                halfedge_t* h = dcel.emplace_halfedge_(n, true);
                n->set_halfedge(h);
                h->set_cell(c);
            }
            
            // twin halfedges for the hole
            for (int i = 0; i < hole_nodes; ++i) {
                node_t* n1 = std::addressof(*(std::next(dcel.nodes_begin(), node_offset + i)));
                node_t* n2 = std::addressof(*(std::next(dcel.nodes_begin(), node_offset + (i + 1) % hole_nodes)));
                halfedge_t* h1 = n1->halfedge();
                halfedge_t* h2 = dcel.emplace_halfedge_(n2, true); // Twin edge
                h2->set_twin(h1);
                h1->set_twin(h2);
                h2->set_cell(nullptr);
            }
            // connect hole's nodes
            for (int i = 0; i < hole_nodes; ++i) {
                halfedge_t* h1 = std::next(dcel.nodes_begin(), node_offset + i)->halfedge();
                halfedge_t* h2 = std::next(dcel.nodes_begin(), node_offset + (i + 1) % hole_nodes)->halfedge();
                h1->set_next(h2);
                h2->set_prev(h1);
                h1->twin()->set_prev(h2->twin());
                h2->twin()->set_next(h1->twin());
            }
            
            node_offset += hole_nodes; 
        }

        return dcel;
    }


    // observers
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> nodes() const {   // matrix of nodes coordinates
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coords(n_nodes_, embed_dim);
        for (int i = 0; i < n_nodes_; ++i) { coords.row(nodes[i].id()) = nodes[i].coords(); }
        return coords;
    }
    int n_nodes() const { return n_nodes_; }
    int n_halfedges() const { return halfedges_.size(); }   //USE size() ?? since we are not changing them
    int n_cells() const { return n_cells_; }
    int n_edges() const { return n_halfedges_ / 2; }


    // iterators
    // cyclic iteration over half-edge chain
    typename halfedge_t::circulator halfedge_circulator(halfedge_t* halfedge) {
        return typename halfedge_t::circulator(halfedge);
    }
    halfedge_iterator halfedges_begin() { return halfedges_.begin(); }
    halfedge_iterator halfedges_end() { return halfedges_.end(); }
    node_iterator nodes_begin() { return nodes_.begin(); }
    node_iterator nodes_end() { return nodes_.end(); }
    cell_iterator cells_begin() { return cells_.begin(); }
    cell_iterator cells_end() { return cells_.end(); }

   
    // to be removed, to export DCEL to file json
    void export_to_json(const std::string& filename) {
        json j;  
        j["nodes"] = json::array();
        for (auto it = nodes_begin(); it != nodes_end(); ++it) {
            json node;
            node["id"] = it->id();
            node["coords"] = {it->coords()(0), it->coords()(1)};
            node["boundary"] = it->on_boundary();
            j["nodes"].push_back(node);
        }  
        j["edges"] = json::array();
        for (auto it = halfedges_begin(); it != halfedges_end(); ++it) {
            json edge;
            edge["id"] = it->id();
            edge["from"] = it->node()->id();
            edge["to"] = it->next()->node()->id();
            edge["twin"] = it->twin() ? it->twin()->id() : -1; 
            j["edges"].push_back(edge);
        }
        j["cells"] = json::array();
        for (auto it = cells_begin(); it != cells_end(); ++it) {
            json cell;
            cell["id"] = it->id();
            cell["edges"] = json::array();
    
            auto h = it->halfedge();
            if (!h) { 
              std::cerr << "Error: cell with null halfedge: cell's ID: " << it->id() << std::endl;
            continue;
            }
            std::cout << "Exporting cell with ID: " << it->id() << std::endl;
            do {
                if (!h) { 
                    std::cerr << "ERROR: null halfedge in cell " << it->id() << std::endl;
                    break;
                }
                cell["edges"].push_back(h->id());
                h = h->next();
            } while (h && h != it->halfedge());
            j["cells"].push_back(cell);
        }
        std::ofstream file(filename);
        file << j.dump(4);
        file.close();
        std::cout << "Exporting DCEL to " << filename << std::endl;
    }


    // modifiers

    void set_n_cells_(int cont){n_cells_ = cont;}

    node_t* insert_node(const node_t& node) {
        nodes_.push_back(node);
	    n_nodes_++;
        return std::addressof(nodes_.back());
    }
    
    //DA RENDERE PRIVATA
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>> extract_internal_paths_with_connections_(const std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>>& boundaries) {

        const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& outer = boundaries[0];
        std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>> results(boundaries.size() - 1); 

        for (size_t r = 1; r < boundaries.size(); ++r) {
            const Eigen::MatrixXd& region = boundaries[r];
            std::list<Eigen::Matrix<double, 1, embed_dim>> internal_points;
            std::list<Eigen::Matrix<double, 1, embed_dim>> boundary_matches;
            bool first_non_match_found = false;
            bool post_match_after_nonmatch = false;
            for (int i = 0; i < region.rows(); ++i) {
                Eigen::Matrix<double, 1, embed_dim> pt = region.row(i);
                bool found = false;
                for (int j = 0; j < outer.rows(); ++j) {
                    if ((pt - outer.row(j)).norm() < 1e-10) {
                        boundary_matches.push_back(pt);
                        found = true;
                        break;
                    }
                }
               if (!first_non_match_found) {
                    if (found) {
                        boundary_matches.push_back(pt);  
                    } else {
                        first_non_match_found = true;
                        internal_points.push_back(boundary_matches.back());
                        internal_points.push_back(pt);
                    }
                } else {
                    if (found) {
                        internal_points.push_back(boundary_matches.back()); 
                        post_match_after_nonmatch = true;
                        break;  // interrompi il ciclo
                    } else {
                        internal_points.push_back(pt); 
                    }
                }
            }
            Eigen::Matrix<double, Eigen::Dynamic, embed_dim> internal_matrix(internal_points.size(), embed_dim);
            for(auto it = internal_points.begin(); it != internal_points.end(); ++it) {
                internal_matrix.row(std::distance(internal_points.begin(), it)) = *it;
            }

            results[r-1]=internal_matrix;
        }

        return results;
    }



 /*   halfedge_t* insert_edge(halfedge_t* v1, halfedge_t* v2) {

        if( (!v1 || !v2) ||
            (v1->cell() && v2-> cell() && v1->cell()!=v2->cell()) ||
            (v2->next() && v1->next() && ( v1->node() == v2->next()->node() || v2->node()==v1->next()->node())) ||
            (v1==v2 || v1->node()==v2->node()) ) {
            return v1;
        }
        // get exiting halfedges from n1 and n2
        node_t* n1 = v1->node();
        node_t* n2 = v2->node();
        // create a pair of twin half-edges
        halfedge_t* h1;
        halfedge_t* h2;
        if(v1->twin())  // if v1 is already structured
            h1 = emplace_halfedge_(n1);
        else  // if v1 was created just to call insert_edge
            h1 = v1;
        if(v2->twin())
            h2 = emplace_halfedge_(n2);
        else
            h2 = v2;

        h1->set_twin(h2);
        h2->set_twin(h1);

        h2->set_next(v1);
        if (v1->prev()){  // if v1 is already structured
            v1->prev()->set_next(h2->twin());
            h2->twin()->set_prev(v1->prev());
            v1->set_prev(h2);
        }
        else{
            h1->set_prev(h2);
            h2->set_next(h1);
        }
        h2->next()->set_prev(h2);
        h2->set_node(n2);

        h1->set_next(v2);
        if (v2->prev()) {  // if v2 is already structured
            v2->prev()->set_next(h1->twin());
            h1->twin()->set_prev(v2->prev());
            v2->set_prev(h1);
        }
        else{
            h2->set_prev(h1);
            h1->set_next(h2);
        }
        h1->next()->set_prev(h1);
        h1->set_node(n1);

        // the newly created halfedges are not one after the other
        if(h1->next() != h2 && h1->prev() != h2){ 
            h1->set_cell(h1->prev()->cell());
            h1->cell()->set_halfedge(h1);
            cells_.push_back(cell_t(n_cells_++));
            cell_t* c1 = std::addressof(cells_.back());
            c1->set_it(std::prev(cells_.end()));
            c1->set_halfedge(h2);
            halfedge_t* end = h2;
            do {         
              h2->set_cell(c1);   
              h2 = h2->next();
            } while (h2 != end );
        }
        // h1 and h2 are one after the other (e.g. an extremity of the mesh)
        else {
            if(h1->next()!=h2)
                h1->set_cell(h1->next()->cell());
            else
                h1->set_cell(h1->prev()->cell());
            h2->set_cell(h1->cell());
        }

        return h1;
    }*/
    halfedge_t* insert_edge(halfedge_t* v1, halfedge_t* v2, bool ignore_diff_cells=false) {

        if( (!v1 || !v2) ||
            (v1->cell() && v2-> cell() && v1->cell()!=v2->cell() && !ignore_diff_cells) ||
            (v2->next() && v1->next() && ( v1->node() == v2->next()->node() || v2->node()==v1->next()->node())) ||
            (v1==v2 || v1->node()==v2->node()) ) {
            return v1;
        }
        
        // get exiting halfedges from n1 and n2
        node_t* n1 = v1->node();
        node_t* n2 = v2->node();
        // create a pair of twin half-edges
        halfedge_t* h1;
        halfedge_t* h2;
        if(v1->twin()!=nullptr)  // if v1 is already structured
            h1 = emplace_halfedge_(n1);
        else  // if v1 was created just to call insert_edge
            h1 = v1;
        if(v2->twin()!=nullptr)
            h2 = emplace_halfedge_(n2);
        else
            h2 = v2;

        h1->set_twin(h2);
        h2->set_twin(h1);

        h2->set_next(v1);
        if (v1->prev()){  // if v1 is already structured
            v1->prev()->set_next(h2->twin());
            h2->twin()->set_prev(v1->prev());
            v1->set_prev(h2);
        }
        else{
            h1->set_prev(h2);
            h2->set_next(h1);
        }
        h2->next()->set_prev(h2);
        h2->set_node(n2);

        h1->set_next(v2);
        if (v2->prev()) {  // if v2 is already structured
            v2->prev()->set_next(h1->twin());
            h1->twin()->set_prev(v2->prev());
            v2->set_prev(h1);
        }
        else{
            h2->set_prev(h1);
            h1->set_next(h2);
        }
        h1->next()->set_prev(h1);
        h1->set_node(n1);

        // the newly created halfedges don't belong to the same cell
        //if(h1->next() != h2 && h1->prev() != h2){ 
        halfedge_t* l=h1;
        do{
            l=l->next();
        }while(l!=h1 && l!=h2);
        if(l==h1){
            h1->set_cell(h1->prev()->cell());
            if(h1->cell())  // h1 is not a boundary edge
                h1->cell()->set_halfedge(h1);
            cells_.push_back(cell_t(n_cells_++));
            cell_t* c1 = std::addressof(cells_.back());
            c1->set_it(std::prev(cells_.end()));
            c1->set_halfedge(h2);
            halfedge_t* end = h2;
            do {         
              h2->set_cell(c1);   
              h2 = h2->next();
            } while (h2 != end );
        }
        // h1 and h2 are one after the other (e.g. an extremity of the mesh)
        else {
            if(h1->next()!=h2)
                h1->set_cell(h1->next()->cell());
            else
                h1->set_cell(h1->prev()->cell());
            h2->set_cell(h1->cell());
        }

        return h1;
    }

    halfedge_t* add_polygon(halfedge_t* v, const std::vector<node_t*>& nodes, bool building_dcel=false){
        int nodes_polygon= nodes.size();                                                                         
        cell_t* c= v->cell();                                                                                    
    
        std::vector<halfedge_t*> halfedges_to_call(nodes_polygon +2 ); 
        halfedges_to_call[0] = v;
        halfedges_to_call[1] = v->next();
    
        // add nodes and create ghost halfedges
        for (int i = 0; i < nodes_polygon; ++i) {
            halfedge_t* h = nullptr;
            if(find_halfedge(nodes[i],c, building_dcel)) // if the halfedge already exists
                    h = find_halfedge(nodes[i],c, building_dcel);
            else{
                    h = emplace_halfedge_(nodes[i]);
                    h->set_cell(c);
            }
            std::cout << h->id() << std::endl;
            halfedges_to_call[i+2] = h;
        }
        // add edges
        for (int i = 0; i < nodes_polygon+2 ; ++i) {
            halfedge_t* h1 = halfedges_to_call[i];
            halfedge_t* h2 = halfedges_to_call[(i + 1) % (nodes_polygon + 2)];
            halfedges_to_call[(i + 1) % (nodes_polygon + 2)] = insert_edge(h1, h2)->next();
        }
        c->set_halfedge(v);
        return v;
    }


 
    halfedge_t* remove_edge(halfedge_t* v1){        
        if(!v1) return nullptr;

        if(v1->cell()==nullptr) // v1 external halfedge on boundary (null cell)
            v1=v1->twin();
        halfedge_t* v2 = v1->twin();
        //std::cout << "v1: " << v1->id() << std::endl;
        //std::cout << "v2: " << v2->id() << std::endl;
        //std::cout << "cell: " << v1->cell()->id() << std::endl;
        halfedge_t* end;
        halfedge_t* begin;
        cell_t* c1 = v1->cell();

        // if v1 is not a boundary edge, insert halfedges of v2's cell to v1's cell
        if(!v1->on_boundary()){
            begin = v2->next();
            end = v2;
            c1->set_halfedge(begin); 
            do {
                begin->set_cell(c1);
                begin = begin->next();
            } while (begin != end);
        }
        // if v1 is a boundary edge, set the halfedges in v1's cell to boundary halfedges,
        // cell=null and remove v1's cell
        else{
            begin= v1->next();
            end = v1;
            do {
                begin->set_cell(nullptr);
                begin->node()->set_boundary(true);
                begin->twin()->node()->set_boundary(true);
                begin = begin->next();
            } while (begin != end);
            // remove v1's cell
            cells_.erase(c1->it());  
        }

        // set next of v1_prev to v2_next
        v1->prev()->set_next(v2->next());      
        v2->next()->set_prev(v1->prev());
        // set next of v2_prev to v1_next
        v2->prev()->set_next(v1->next());
        v1->next()->set_prev(v2->prev());

        // remove v2's cell
        cell_t* c2=v2->cell();
        if (c2) {
            cells_.erase(c2->it());
        }
        // remove v1 and v2
        halfedge_t* next= v2->next();
        halfedges_.erase(v1->it());  
        halfedges_.erase(v2->it());
        return next; 
    }


    // remove polygon inside dcel by calling its id 
    void remove_polygon(int cell_id) {

        auto it = cells_.begin();
        for (; it != cells_.end(); ++it) 
            if (it->id() == cell_id) 
                break;
        cell_t* cell= &(*it);
        if (!cell) return;
    
        // removing cell's halfedges
        halfedge_t* h1 = cell->halfedge();
        halfedge_t* ending = h1->twin() ? h1->twin()->next() : nullptr; 
        do {
            halfedge_t* next = remove_edge(h1->twin());
            h1 = next; 
        } while (h1 && h1 != ending);
    
    }
        

    // return the node of the halfedge previous to h 
    node_t* adjacent(halfedge_t* h) const {return (h->twin()) ? h->twin()->prev()->node() : nullptr;  }

    // find the halfedge given its node and cell
    halfedge_t* find_halfedge(node_t* n, cell_t* cell, bool building_dcel=false) {
        if(!building_dcel){
            halfedge_t* h= cell->halfedge();
            halfedge_t* end=h;
            do{
                if(h->node()==n)
                    return h;
                h = h->next();
            }while(h!=end);
        }
        else{
            for(auto it = halfedges_begin(); it != halfedges_end(); ++it) {
                halfedge_t* h = &(*it);
                if (h->node() == n && h->cell() == cell) {
                    return h;
                }
            }
        }
        return nullptr;
    }

    halfedge_t* find_halfedge_between(coords_t from, coords_t to) {
        for (auto it = halfedges_begin(); it != halfedges_end(); ++it) {
            halfedge_t* h = &(*it);
            if (h->node()->coords() == from && h->twin() && h->twin()->node()->coords() == to) {
                return h;
            }
        }
        return nullptr;
    }

    node_t* find_node(const coords_t& coords) {
        for (auto it = nodes_begin(); it != nodes_end(); ++it) {
            if (it->coords() == coords) {
                return std::addressof(*it);
            }
        }
        return nullptr;
    }
    
  

    // internal utils
    template <typename... Args> halfedge_t* emplace_halfedge_(Args&&... args) {
        halfedges_.emplace_back(n_halfedges_++, std::forward<Args>(args)...);
        auto it = std::prev(halfedges_.end()); // oppure std::next(halfedges_.begin(), n_halfedges_ - 1);
        it->set_it(it); // = it;
        return std::addressof(halfedges_.back());
    }

    template <typename... Args> node_t* emplace_node_(Args&&... args) {
        nodes_.emplace_back(n_nodes_++, std::forward<Args>(args)...);
        return std::addressof(nodes_.back());
    }  

    template <typename TriangulationType>
    TriangulationType to_triangulation() const {
        Eigen::Matrix<double, Eigen::Dynamic, embed_dim> nodes_mat(n_nodes(), embed_dim);
        Eigen::Matrix<int, Eigen::Dynamic, 3> cells_mat(n_cells(), 3);
        Eigen::Matrix<int, Eigen::Dynamic, 1> boundary_markers(n_nodes());

        int idx = 0;
        for (auto it = nodes_.begin(); it != nodes_.end(); ++it, ++idx) {
            nodes_mat.row(it->id()) = it->coords().transpose();
            boundary_markers(it->id()) = it->on_boundary() ? 1 : 0;
        }

        idx = 0;
        for (auto it = cells_.begin(); it != cells_.end(); ++it, ++idx) {
            halfedge_t* h = it->halfedge();
            cells_mat(idx, 0) = h->node()->id();
            cells_mat(idx, 1) = h->next()->node()->id();
            cells_mat(idx, 2) = h->prev()->node()->id();
        }
        return TriangulationType(nodes_mat, cells_mat, boundary_markers);
    }

    template <typename TriangulationType>
    void from_triangulation(const TriangulationType& triangulation, const std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>>& holes) {
        int n_hole_nodes = 0;
        for (const auto& hole : holes)
            n_hole_nodes += hole.rows();
        int n_boundary_external = triangulation.n_boundary_nodes() - n_hole_nodes;
        Eigen::Matrix<double, Eigen::Dynamic, embed_dim> boundary_nodes(n_boundary_external, embed_dim);
        const auto& coords = triangulation.nodes();
        const auto& markers = triangulation.boundary_nodes();
        int n = boundary_nodes.rows();
        
        int idx = 0;
        for (int i = 0; i < n_boundary_external; ++i) {
            if (markers(i, 0) == 1) {
                boundary_nodes.row(idx++) = coords.row(i);
            }
        }
        if(this->nodes_.empty() && this->halfedges_.empty() && this->cells_.empty()) { 
            *this = DCEL::make_polygon(boundary_nodes, holes);
        }
        else{
            std::cout << "Warning: DCEL is not empty, cannot overwrite with triangulation data." << std::endl;
        }
        for (int i = 0; i < coords.rows(); ++i) {
            if (markers(i, 0) == 0) {
                insert_node(typename DCEL::node_t(n_nodes(), false, coords.row(i)));
            }
        }

        auto cells= triangulation.cells();

        // creating a map to associate nodes to holes
        std::unordered_map<int, int> node_to_hole;
        std::vector<std::set<std::pair<int, int>>> hole_edges(holes.size());
        for (int k = 0; k < holes.size(); ++k) {
            const auto& hole = holes[k];
            std::vector<int> hole_node_ids;  
            for (int i = 0; i < hole.rows(); ++i) {
                const auto& pt = hole.row(i);
                for (int j = 0; j < triangulation.nodes().rows(); ++j) {
                    if ((triangulation.nodes().row(j) - pt).norm() < 1e-10) {
                        node_to_hole[j] = k;
                        hole_node_ids.push_back(j);
                        break;
                    }
                }
            }
            int m = hole_node_ids.size();
            for (int i = 0; i < m; ++i) {
                int a = hole_node_ids[i];
                int b = hole_node_ids[(i + 1) % m];  
                hole_edges[k].insert({a, b});  
            }
        }

        // vector of bools to understand if an edge of a hole has been connected to the rest of the dcel being created
        std::vector<bool> is_connected(holes.size(),false);
        int cont=0;
        cell_t* longest_cell = nullptr;
        std::list< Eigen::Matrix<int, Eigen::Dynamic, 3> > cells_list;
        for(int i = 0; i < triangulation.n_cells(); ++i) {
            cells_list.push_back(Eigen::Matrix<int, Eigen::Dynamic, 3>(triangulation.cells().row(i)));
        }
        
        while(!cells_list.empty()) {

            auto cell = cells_list.front();
            cells_list.pop_front();
            
            int id0 = cell(0, 0);
            int id1 = cell(0, 1);
            int id2 = cell(0, 2);

            // first, connect for each hole a triangle that has an edge on the hole, to better deal with the new cells being created 
            bool all_connected = std::all_of(is_connected.begin(), is_connected.end(), [](bool v) { return v; });
            if (!all_connected){
                // map node -> hole (-1 if external)
                int h0 = node_to_hole.count(id0) ? node_to_hole.at(id0) : -1;
                int h1 = node_to_hole.count(id1) ? node_to_hole.at(id1) : -1;
                int h2 = node_to_hole.count(id2) ? node_to_hole.at(id2) : -1;

                // verify if it's a priority triangle (one that connects the hole to the boundary)
                bool is_priority = false;
                int hole_to_connect = -1;

                // to check if 2 nodes in same hole are actually consecutive
                auto is_consecutive_in_hole = [&](int a, int b, int hole_id) {
                    return hole_edges[hole_id].count({a, b}) || hole_edges[hole_id].count({b, a});
                };

                if (h0 == h1 && h0 != -1 && h2 == -1 && !is_connected[h0] && is_consecutive_in_hole(id0, id1, h0)) {
                    is_priority = true;
                    hole_to_connect = h0;
                }
                else if (h1 == h2 && h1 != -1 && h0 == -1 && !is_connected[h1] && is_consecutive_in_hole(id1, id2, h1)) {
                    is_priority = true;
                    hole_to_connect = h1;
                }
                else if (h2 == h0 && h2 != -1 && h1 == -1 && !is_connected[h2] && is_consecutive_in_hole(id2, id0, h2)) {
                    is_priority = true;
                    hole_to_connect = h2;
                }
                
                if (!is_priority) {
                    // go to the next cell
                    cells_list.push_back(cell);
                    continue;
                }
                is_connected[hole_to_connect] = true;
            }

            /*auto n0 = std::find_if(nodes_.begin(), nodes_.end(),
                [&](const node_t& n) { return n.id() == id0; });
            auto n1 = std::find_if(nodes_.begin(), nodes_.end(),
                [&](const node_t& n) { return n.id() == id1; });
            auto n2 = std::find_if(nodes_.begin(), nodes_.end(),
                [&](const node_t& n) { return n.id() == id2; });

            node_t* p0 = std::addressof(*n0);
            node_t* p1 = std::addressof(*n1);
            node_t* p2 = std::addressof(*n2);*/

            int row0 = cell(0, 0);
            int row1 = cell(0, 1);
            int row2 = cell(0, 2);
            node_t* p0 = find_node(coords.row(row0));
            node_t* p1 = find_node(coords.row(row1));
            node_t* p2 = find_node(coords.row(row2));

            if (!fdapde::internals::are_2d_counterclockwise_sorted(p0->coords(), p1->coords(), p2->coords())) {
                std::swap(p1, p2);
                std::cout<< "swapping p1 and p2" << std::endl;
            }
            std::cout << "p0: " << p0->id() << " p1: " << p1->id() << " p2: " << p2->id() << std::endl;
            std::cout << " p0: " << p0->coords().transpose() << " p1: " << p1->coords().transpose() << " p2: " << p2->coords().transpose() << std::endl;
            halfedge_t* h = find_halfedge_between(p0->coords(), p1->coords());
            bool building_dcel = true; 
            if (h) {
                std::cout << "qui 1"    << std::endl;
                add_polygon(h, {p2}, building_dcel);
            } else if ((h = find_halfedge_between(p1->coords(), p2->coords()))) {
                std::cout << "qui 2"    << std::endl;
                add_polygon(h, {p0}, building_dcel);
            } else if ((h = find_halfedge_between(p2->coords(), p0->coords()))) {
                std::cout << "qui 3"    << std::endl;
                add_polygon(h, {p1}, building_dcel);
            } else {
                // there might be 2 or more halfedges departing from same node and belonging to same cell
                // e.g. if edges belong to boundary, hole and another hole
                std::cout << "qui 4"    << std::endl;
                cells_list.push_back(cell);
                continue;
            }
              
            if(!all_connected){
                // find the cell that connects all the halfedges that don't belong to triangles yet 
                for(auto it = cells_begin(); it!= cells_end(); ++it){
                    int cont=0;
                    halfedge_t* h= it->halfedge();
                    do{
                        cont++;
                        h=h->next();
                    }while(h!=it->halfedge());
                    if(cont>3){  // cell is not a triangle
                        longest_cell= &(*it);
                        break;
                    }

                }

                auto it = halfedges_.begin();
                // go to the first half-edge of the internal boundary (hole)
                for (int i = 0; i < n_boundary_external*2; ++i)
                    ++it;
                for (int k = 0; k < holes.size(); ++k) {
                    // if the hole is not connected yet, change the cell of the halfedges to longest_cell
                    if (!is_connected[k]) { 
                        for (int h = 0; h < holes[k].rows(); ++h) {
                            it->set_cell(longest_cell);  
                            ++it;
                        }
                        // disregard the twins since they must have null cell
                        for (int h = 0; h < holes[k].rows(); ++h) {
                            it++;
                        }
                    } else {
                        // skip halfedges of the hole
                        for (int h = 0; h < holes[k].rows()*2; ++h)
                            ++it;
                    }
                }
            }
            cont++;
            //if(cont==1) break;
            
        }
    }

private:
    // internal storage (use list to avoid reallocations)
    std::list<node_t> nodes_;
    std::list<halfedge_t> halfedges_;
    std::list<cell_t> cells_;
    int n_nodes_, n_halfedges_, n_cells_;
};


}   // namespace fdapde



#endif // __DCEL_H__
