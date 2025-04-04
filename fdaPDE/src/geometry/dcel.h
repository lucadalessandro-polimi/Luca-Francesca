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
        public:
         halfedge_t() : node_(nullptr), prev_(nullptr), next_(nullptr), twin_(nullptr) { }
         halfedge_t(int id, halfedge_t* prev, halfedge_t* next, halfedge_t* twin, node_t* node) :
             id_(id), prev_(prev), next_(next), twin_(twin), node_(node) { }
         // no twin constructors
         halfedge_t(int id, halfedge_t* prev, halfedge_t* next, node_t* node) :
             halfedge_t(id, prev, next, nullptr, node) { }
         // minimal constructor
         halfedge_t(int id, node_t* node) : halfedge_t(id, nullptr, nullptr, nullptr, node) { }
 
         // observers
         halfedge_t* prev() const { return prev_; }
         halfedge_t* next() const { return next_; }
         halfedge_t* twin() const { return twin_; }
         node_t* node() const { return node_; }
         cell_t* cell() const { return cell_; }
         int id() const { return id_; }
         bool on_boundary() const { return (node_->on_boundary() && twin_->node()->on_boundary() && (cell()==nullptr || twin()->cell()==nullptr)); }
         std::list<halfedge_t>::iterator it() const { return it_; }
         // modifiers
         void set_prev(halfedge_t* prev) { prev_ = prev; }
         void set_next(halfedge_t* next) { next_ = next; }
         void set_twin(halfedge_t* twin) { twin_ = twin; }
         void set_node(node_t* node) { node_ = node; }
         void set_cell(cell_t* cell) { cell_ = cell; }
         void set_id(int id) {id_=id;}
         void set_it(std::list<halfedge_t>::iterator it) { it_ = it; }
 
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
    static DCEL<local_dim, embed_dim> make_polygon(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& nodes) {
        fdapde_assert(nodes.cols() == embed_dim);
        int n_nodes = nodes.rows();
        int n_halfedges = 2 * (n_nodes);
        DCEL<local_dim, embed_dim> dcel;
        // create polygon cell
        dcel.cells_.push_back(cell_t(0));
        cell_t* c = std::addressof(dcel.cells_.back());
        dcel.n_cells_ = 1;
        // push nodes
        for (int i = 0; i < n_nodes; ++i) {
            node_t* n = dcel.insert_node(node_t(i, /* boundary = */ true, nodes.row(i)));
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
    }
    

    // overloading of make_polygon to also add holes
    static DCEL<local_dim, embed_dim> make_polygon(
        const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& boundary,
        const std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>>& holes) {
    
        fdapde_assert(boundary.cols() == embed_dim);
        DCEL<local_dim, embed_dim> dcel;
    
        // external boundary
        int n_nodes = boundary.rows();
        dcel.cells_.push_back(cell_t(0)); // Cella principale
        cell_t* c = std::addressof(dcel.cells_.back());
        dcel.n_cells_ = 1;
        // nodes and halfedges for external boundary
        for (int i = 0; i < n_nodes; ++i) {
            node_t* n = dcel.insert_node(node_t(i, true, boundary.row(i)));
            halfedge_t* h = dcel.emplace_halfedge_(n);
            n->set_halfedge(h);
            h->set_cell(c);
        }
        c->set_halfedge(dcel.nodes_begin()->halfedge());
        // twin halfedges for external boundary
        for (auto it = dcel.nodes_begin(); it != dcel.nodes_end(); ++it) {
            node_t* n1 = std::addressof(*it);
            node_t* n2 = std::addressof(*((it->id() == n_nodes - 1) ? dcel.nodes_begin() : std::next(it, 1)));
            halfedge_t* h1 = n1->halfedge();
            halfedge_t* h2 = dcel.emplace_halfedge_(n2); // Twin edge
    
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
    
        // adding holes
        int node_offset = n_nodes; 
        int hole_index = 1;
    
        for (const auto& hole : holes) {
            int hole_nodes = hole.rows();
            // dcel.cells_.push_back(cell_t(hole_index++)); // Nuova cella per il buco
            // cell_t* hole_cell = std::addressof(dcel.cells_.back());
            // nodes and halfedges for the hole
            for (int i = 0; i < hole_nodes; ++i) {
                node_t* n = dcel.insert_node(node_t(node_offset + i,  true, hole.row(i)));
                halfedge_t* h = dcel.emplace_halfedge_(n);
                n->set_halfedge(h);
                h->set_cell(c);
               // std::cout << "HALF CHE STO CREANDO: " << h->id() << " ASSEGNATO A CELLA: " << h->cell()->id() << std::endl;
            }
            // hole_cell->set_halfedge(std::next(dcel.nodes_begin(), node_offset)->halfedge());
            // twin halfedges for the hole
            for (int i = 0; i < hole_nodes; ++i) {
                node_t* n1 = std::addressof(*(std::next(dcel.nodes_begin(), node_offset + i)));
                node_t* n2 = std::addressof(*(std::next(dcel.nodes_begin(), node_offset + (i + 1) % hole_nodes)));
                halfedge_t* h1 = n1->halfedge();
                halfedge_t* h2 = dcel.emplace_halfedge_(n2); // Twin edge
                h2->set_twin(h1);
                h1->set_twin(h2);
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
    int n_halfedges() const { return n_halfedges_; }   //USE size() ?? since we are not changing them
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
            int count = 0;
            do {
                if (!h) { 
                    std::cerr << "ERROR: null halfedge in cell " << it->id() << std::endl;
                    break;
                }
                if (count > 20) {
                    std::cerr << "Infinite loop in cell " << it->id() << std::endl;
                    break;
                }
                cell["edges"].push_back(h->id());
                h = h->next();
                count++;
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


    halfedge_t* insert_edge(halfedge_t* v1, halfedge_t* v2) {

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
    halfedge_t* add_polygon(halfedge_t* v, const std::vector<node_t*>& nodes){
      int nodes_polygon= nodes.size();                                                                         
        cell_t* c= v->cell();                                                                                    
    
        std::vector<halfedge_t*> halfedges_to_call(nodes_polygon +2 ); 
        halfedges_to_call[0] = v;
        halfedges_to_call[1] = v->next();
    
        // add nodes and create ghost halfedges
       for (int i = 0; i < nodes_polygon; ++i) {
            halfedge_t* h = nullptr;
            if(find_halfedge(nodes[i],c))
                    h = find_halfedge(nodes[i],c);
            else{
                    h = emplace_halfedge_(nodes[i]);
                    h->set_cell(c);
            } 
            halfedges_to_call[i+2] = h;
        }
        // add edges
        for (int i = 0; i < nodes_polygon+2 ; ++i) {
            halfedge_t* h1 = halfedges_to_call[i];
            halfedge_t* h2 = halfedges_to_call[(i + 1) % (nodes_polygon + 2)];
            halfedges_to_call[(i + 1) % (nodes_polygon + 2)] = insert_edge(h1, h2)->next();
        }

        c->set_halfedge(v);
        return c->halfedge();
    }


 
    halfedge_t* remove_edge(halfedge_t* v1){        
        if(!v1) return nullptr;

        if(!v1->cell()) // v1 external halfedge on boundary (null cell)
            v1=v1->twin();

        halfedge_t* v2 = v1->twin();

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
            cell_iterator it = cells_.begin();
            while (it != cells_.end()) {
                if (it->id() == c1->id()) 
                    break;  
                ++it;
            }
            if (it != cells_.end()) 
                cells_.erase(it);
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
            cell_iterator it = cells_.begin();
            while (it != cells_.end()) {
                if (it->id() == c2->id()) 
                    break;  
                ++it;
            }
            if (it != cells_.end()) 
                cells_.erase(it);
        }
        // remove v1 and v2
        halfedge_t* next= v2->next();
        if (v1 && v2) {
            halfedge_iterator it1 = halfedges_.begin();
            while (it1 != halfedges_.end()) {
                if (it1->id() == v1->id()) {
                    break;
                }
                ++it1;
            }
            halfedge_iterator it2 = halfedges_.begin();
            while (it2 != halfedges_.end()) {
                if (it2->id() == v2->id()) {
                    break;
                }
                ++it2;
            }
            if (it1 != halfedges_.end() && it2 != halfedges_.end()) {
                halfedges_.erase(it1);  
                halfedges_.erase(it2);
            }
        }
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
    halfedge_t* find_halfedge(node_t* n, cell_t* cell) {
        halfedge_t* h= cell->halfedge();
        halfedge_t* end=h;
        do{
            if(h->node()==n)
                return h;
            h = h->next();
        }while(h!=end);
        return nullptr;
    }
    halfedge_t* find_halfedge_between(node_t* from, node_t* to) {
        for (auto it = halfedges_begin(); it != halfedges_end(); ++it) {
            halfedge_t* h = &(*it);
            if (h->node() == from && h->twin() && h->twin()->node() == to) {
                return h;
            }
        }
        return nullptr;
    }
    
  

    // internal utils
    template <typename... Args> halfedge_t* emplace_halfedge_(Args&&... args) {
        halfedges_.emplace_back(n_halfedges_++, std::forward<Args>(args)...);
        return std::addressof(halfedges_.back());
    }

    template <typename... Args> node_t* emplace_node_(Args&&... args) {
        nodes_.emplace_back(n_nodes_++, std::forward<Args>(args)...);
        return std::addressof(nodes_.back());
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
