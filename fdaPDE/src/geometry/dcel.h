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
        //using coords_t = Eigen::Matrix<double, embed_dim, 1>;
        int id_;                  // global node index
        halfedge_t* halfedge_;    // any edge having this node as its origin
        bool boundary_;           // asserted true if node is on boundary
        coords_t coords_;
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


        // observers
        const Eigen::Matrix<double, embed_dim, 1>& coords() const { return coords_; }
        halfedge_t* halfedge() const { return halfedge_; }
        void set_halfedge(halfedge_t* halfedge) { halfedge_ = halfedge; }
        int id() const { return id_; }
        bool on_boundary() const { return boundary_; }
        node_t* next() const { return halfedge_->next()->node(); }
        node_t* prev() const { return halfedge_->prev()->node(); }

    };
    struct halfedge_t {
       private:
        int id_;   // global halfedge index
        halfedge_t *prev_, *next_, *twin_;
        node_t* node_;
        cell_t* cell_;   // cell to which this halfedge belongs to
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
        // modifiers
        void set_prev(halfedge_t* prev) { prev_ = prev; }
        void set_next(halfedge_t* next) { next_ = next; }
        void set_twin(halfedge_t* twin) { twin_ = twin; }
        void set_node(node_t* node) { node_ = node; }
        void set_cell(cell_t* cell) { cell_ = cell; }
        void set_id(int id) {id_=id;}

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
        // observers
        halfedge_t* halfedge() const { return h_; }
        int id() const { return id_; }
        void set_id(int id) {id_=id;}
        // modifiers
        void set_halfedge(halfedge_t* h) { h_ = h; }
           // Operatore di uguaglianza per confrontare due celle
        bool operator==(const cell_t& other) const {
          return id_ == other.id_;  // Confrontiamo solo l'ID, che dovrebbe essere univoco
        }

       private:
        int id_;
        halfedge_t* h_;
    };


    using halfedge_iterator = std::list<halfedge_t>::iterator;
    using node_iterator = std::list<node_t>::iterator;
    using cell_iterator = std::list<cell_t>::iterator;
    //per fare le funzioni costanti lavorando con gli iteratori
    using const_cell_iterator = std::list<cell_t>::const_iterator;

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

    
    // modifiers
    node_t* insert_node(const node_t& node) {
        nodes_.push_back(node);
	    n_nodes_++;
        return std::addressof(nodes_.back());
    }

    // observers
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> nodes() const {   // matrix of nodes coordinates
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coords(n_nodes_, embed_dim);
        for (int i = 0; i < n_nodes_; ++i) { coords.row(nodes[i].id()) = nodes[i].coords(); }
        return coords;
    }
    int n_nodes() const { return n_nodes_; }
    int n_halfedges() const { return n_halfedges_; }
    int n_cells() const { return n_cells_; }
    int n_edges() const { return n_halfedges_ / 2; }

    void set_n_cells_(int num){n_cells_=num;}

    //non li facciamo funzionare quindi accedere con size !!!


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


   //////////////////////////////////// FUNZIONI NUOVE DI DCEL /////////////////////////////////////////
   const_cell_iterator cells_begin() const { return cells_.cbegin(); }
   const_cell_iterator cells_end() const { return cells_.cend(); }
   
    void export_to_json(const std::string& filename) {
        json j;  
        // Salva i nodi
        j["nodes"] = json::array();
        for (auto it = nodes_begin(); it != nodes_end(); ++it) {
            json node;
            node["id"] = it->id();
            node["coords"] = {it->coords()(0), it->coords()(1)};
            node["boundary"] = it->on_boundary();
            j["nodes"].push_back(node);
        }  
        // Salva gli archi
        j["edges"] = json::array();
        for (auto it = halfedges_begin(); it != halfedges_end(); ++it) {
            json edge;
            edge["id"] = it->id();
            edge["from"] = it->node()->id();
            edge["to"] = it->next()->node()->id();
            edge["twin"] = it->twin() ? it->twin()->id() : -1;  // -1 se non ha twin
            j["edges"].push_back(edge);
        }
        // Salva le celle
        j["cells"] = json::array();
        for (auto it = cells_begin(); it != cells_end(); ++it) {
            json cell;
            cell["id"] = it->id();
            cell["edges"] = json::array();
    
            auto h = it->halfedge();
            if (!h) { // Se h è nullptr, saltiamo questa cella per evitare crash
              std::cerr << "Errore: cella con halfedge nullo!ID Cella: " << it->id() << std::endl;
            continue;
            }
            std::cout << "Esportando cella ID: " << it->id() << std::endl;
            int count = 0;
            /*do {
              if (h) {  // Controllo extra per sicurezza
                cell["edges"].push_back(h->id());}
              h = h->next();
            } while (h && h != it->halfedge());*/
            do {
                if (!h) { 
                    std::cerr << "ERRORE: Half-edge nullo in cella " << it->id() << std::endl;
                    break;
                }
                if (count > 20) {
                    std::cerr << "Loop infinito rilevato nella cella " << it->id() << std::endl;
                    break;
                }
                cell["edges"].push_back(h->id());
                h = h->next();
                count++;
            } while (h && h != it->halfedge());
            j["cells"].push_back(cell);
        }
        std::ofstream file(filename);
        file << j.dump(4); // Indentazione di 4 spazi per leggibilità
        file.close();
        std::cout << "DCEL esportata in " << filename << std::endl;
    }

    
    halfedge_t* remove_edge(halfedge_t* v1){        
        if(!v1) return nullptr;

        //if(v1->on_boundary() || v1->twin()->on_boundary()){
        //    std::cout << "Edge on boundary. Removing is not allowed" << std::endl;
        //    return nullptr;
        //}
        if(!v1->cell()) //v1 external halfedge on boundary (null cell)
            v1=v1->twin();
        halfedge_t* v2 = v1->twin();
        // insert halfedges in v1's cell
        halfedge_t* end = v2;
        halfedge_t* begin = v2->next();
        cell_t* c1 = v1->cell();
        c1->set_halfedge(begin); // either this or v1->prev()
        do {
            begin->set_cell(c1);
            begin = begin->next();
        } while (begin != end);
        // set next of v1_prev to v2_next
        v1->prev()->set_next(v2->next());      
        v2->next()->set_prev(v1->prev());
        // set next of v2_prev to v1_next
        v2->prev()->set_next(v1->next());
        v1->next()->set_prev(v2->prev());
        // remove v2's cell
        //n_cells_--;
        cell_t* c2=v2->cell();
        if (c2) {
            //auto it = std::find(cells_begin(), cells_end(), v2->cell());
            //auto it = std::find_if(cells_.begin(), cells_.end(), [=](const cell_t& c) { return c.id() == c2->id(); });
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
        halfedge_t* next= v1->next();
        //n_halfedges_ -= 2;
        if (v1 && v2) {
            //auto it1 = std::find(halfedges_begin(), halfedges_end(), v1);
            //auto it2 = std::find(halfedges_begin(), halfedges_end(), v2);
            //auto it1 = std::find_if(halfedges_.begin(), halfedges_.end(), [=](const halfedge_t& h) { return h.id() == some_halfedge_ptr1->id();});
            //auto it2 = std::find_if(halfedges_.begin(), halfedges_.end(), [=](const halfedge_t& h) { return h.id() == some_halfedge_ptr2->id();});
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

    void remove_polygon(const cell_t* cell) {  
        if (!cell) return;
        
        halfedge_t* temp = cell->halfedge();
        int safety_counter = 0; 
        do {
            std::cout << temp->id() << " ";
            temp = temp->next();
            if (++safety_counter > 500) {    //USEFUL?
                std::cerr << "ERROR: infinite loop in the initial scan of cell " << cell->id() << std::endl;
                return;
            }
        } while (temp != cell->halfedge());
        std::cout << std::endl;
    
        // Rimozione degli half-edge
        halfedge_t* h1 = cell->halfedge();
        halfedge_t* ending = h1->twin() ? h1->twin()->next() : nullptr;
    
        safety_counter = 0;
        do {
            std::cout << "Removing half-edge ID: " << h1->id() << " (next will be "
                      << (h1->next() ? std::to_string(h1->next()->id()) : "nullptr") << ")" << std::endl;
            halfedge_t* next = remove_edge(h1);
            if (!next || next == h1) {        //USEFUL?
                std::cout << "Interrupt the removing: next halfedge is null or identical to the current one.\n";
                break;
            }
            h1 = next; 
            if (++safety_counter > 500) {
                std::cerr << "ERROR: infinite loop while removing cell " << cell->id() << std::endl;
                return;
            }
        } while (h1 && h1 != ending);
    
        auto it = std::find_if(cells_.begin(), cells_.end(), [&](const cell_t& c) {
            return c.id() == cell->id();
        });
    
        if (it != cells_.end()) 
            cells_.erase(it);
    }


    static DCEL<local_dim, embed_dim> make_polygon(
        const Eigen::Matrix<double, Eigen::Dynamic, embed_dim>& boundary,
        const std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>>& holes) {
    
        fdapde_assert(boundary.cols() == embed_dim);
        DCEL<local_dim, embed_dim> dcel;
    
        // Creazione del bordo esterno
        int n_nodes = boundary.rows();
        dcel.cells_.push_back(cell_t(0)); // Cella principale
        cell_t* c = std::addressof(dcel.cells_.back());
        dcel.n_cells_ = 1;
    
        // Aggiunta dei nodi e creazione degli half-edge per il bordo esterno
        for (int i = 0; i < n_nodes; ++i) {
            node_t* n = dcel.insert_node(node_t(i, true, boundary.row(i)));
            halfedge_t* h = dcel.emplace_halfedge_(n);
            n->set_halfedge(h);
            h->set_cell(c);
        }

        c->set_halfedge(dcel.nodes_begin()->halfedge());
    
        // Creazione delle twin edges per il bordo esterno
        for (auto it = dcel.nodes_begin(); it != dcel.nodes_end(); ++it) {
            node_t* n1 = std::addressof(*it);
            node_t* n2 = std::addressof(*((it->id() == n_nodes - 1) ? dcel.nodes_begin() : std::next(it, 1)));
            halfedge_t* h1 = n1->halfedge();
            halfedge_t* h2 = dcel.emplace_halfedge_(n2); // Twin edge
    
            h2->set_twin(h1);
            h1->set_twin(h2);
        }
    
        // Impostiamo next e prev per il bordo esterno
        for (auto it = dcel.nodes_begin(); it != dcel.nodes_end(); ++it) {
            halfedge_t* h1 = it->halfedge();
            halfedge_t* h2 = ((it->id() == n_nodes - 1) ? dcel.nodes_begin() : std::next(it, 1))->halfedge();
            h1->set_next(h2);
            h2->set_prev(h1);
            h1->twin()->set_prev(h2->twin());
            h2->twin()->set_next(h1->twin());
        }
    
        // Aggiunta dei buchi
        int node_offset = n_nodes; // Indice iniziale per i nodi dei buchi
        int hole_index = 1;
    
        for (const auto& hole : holes) {
            int hole_nodes = hole.rows();
        //    dcel.cells_.push_back(cell_t(hole_index++)); // Nuova cella per il buco
        //    cell_t* hole_cell = std::addressof(dcel.cells_.back());
    
            // Creazione dei nodi e degli half-edge per il buco
            for (int i = 0; i < hole_nodes; ++i) {
                node_t* n = dcel.insert_node(node_t(node_offset + i,  true, hole.row(i)));
                halfedge_t* h = dcel.emplace_halfedge_(n);
                n->set_halfedge(h);
                h->set_cell(c);
               // std::cout << "HALF CHE STO CREANDO: " << h->id() << " ASSEGNATO A CELLA: " << h->cell()->id() << std::endl;
            }

           // hole_cell->set_halfedge(std::next(dcel.nodes_begin(), node_offset)->halfedge());
    
            // Creazione delle twin edges per il buco
            for (int i = 0; i < hole_nodes; ++i) {
                node_t* n1 = std::addressof(*(std::next(dcel.nodes_begin(), node_offset + i)));
                node_t* n2 = std::addressof(*(std::next(dcel.nodes_begin(), node_offset + (i + 1) % hole_nodes)));
                halfedge_t* h1 = n1->halfedge();
                halfedge_t* h2 = dcel.emplace_halfedge_(n2); // Twin edge
    
                h2->set_twin(h1);
                h1->set_twin(h2);
            }
    
            // Collegare i nodi del buco
            for (int i = 0; i < hole_nodes; ++i) {
                halfedge_t* h1 = std::next(dcel.nodes_begin(), node_offset + i)->halfedge();
                halfedge_t* h2 = std::next(dcel.nodes_begin(), node_offset + (i + 1) % hole_nodes)->halfedge();
                h1->set_next(h2);
                h2->set_prev(h1);
                h1->twin()->set_prev(h2->twin());
                h2->twin()->set_next(h1->twin());
            }

            

            node_offset += hole_nodes; // Aggiorna l'offset per il prossimo buco
        }
    
        return dcel;
    }
        

    node_t* adjacent(halfedge_t* edge) const {return (edge->twin()) ? edge->twin()->prev()->node() : nullptr;  }

    halfedge_t* insert_edge(halfedge_t* v1, halfedge_t* v2) {
        if (v1->cell() && v2-> cell() && v1->cell()!=v2->cell()){
            std::cout<< "cella " << v1->next()->id() << std::endl;
            std::cout<< "cella " << v1->next()->next()->id() << std::endl;
            std::cout<< "cella " << v1->next()->next()->next()->id() << std::endl;
            std::cout<< "cella " << v1->next()->next()->next()->next()->id() << std::endl;
            std::cout<< "cella " << v1->next()->next()->next()->next()->next()->id() << std::endl;
            std::cout << v1->node()->id() << std::endl;
            std::cout << v2->node()->id() << std::endl;
            std::cout << "halfedge v1 : "<<v1->id() << std::endl;
            std::cout << "halfedge v2 : "<<v2->id() << std::endl;
            std::cout << "cella1 " << v1->cell()->id() << std::endl;
            std::cout << "cella2 " << v2->cell()->id() << std::endl;
            std::cerr << "Error: the two halfedges belong to different cells. It's not allowed to insert an edge in between them" << std::endl;
            return v1;
        }
        if (v1 && v2 && v2->next() && v1->next() && ( v1->node() == v2->next()->node() || v2->node()==v1->next()->node()) ) {
           // std::cout << "Consecutive halfedges" << std::endl;
            return v1;   
        }
        // get exiting halfedges from n1 and n2
        node_t* n1 = v1->node();
        node_t* n2 = v2->node();
        // create a pair of twin half-edges
        halfedge_t* h1;
        halfedge_t* h2;
        if(v1->twin()){  //if the halfedge is already structured
            h1 = emplace_halfedge_(n1);
        }
        else{  //if the halfedge was created just to call insert_edge
            h1 = v1;
        }
        if(v2->twin()){
            h2 = emplace_halfedge_(n2);
        }
        else{
            h2 = v2;
        }
        h1->set_twin(h2);
        h2->set_twin(h1);
        h2->set_next(v1);
        if (v1->prev()){
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
        if (v2->prev()) {
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
        else {
            if(h1->next()!=h2)
                h1->set_cell(h1->next()->cell());
            else
                h1->set_cell(h1->prev()->cell());
            h2->set_cell(h1->cell());
        }
        return h1;
    }


    halfedge_t* add_polygon(halfedge_t* v, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& nodes){  //DARE IN INGRESSO ANCHE I BOUNDARY EDGES...
        int nodes_polygon= nodes.rows();                                                                         //passare boundary edges se ci sono....
        cell_t* c= v->cell();                                                                                    //ma per i poligoni va bene??
    
        std::vector<halfedge_t*> ghost_halfedges(nodes_polygon +2 ); //O(n)
        ghost_halfedges[0] = v;
        ghost_halfedges[1] = v->next();
    
        // add nodes and create ghost halfedges
       for (int i = 0; i < nodes_polygon; ++i) {
            node_t* n;
            if(find_node(nodes.row(i))==nullptr){  
                n = insert_node(node_t(n_nodes_, /* boundary = */ false, nodes.row(i)));
            }
            else {
                n = find_node(nodes.row(i));  //if node is already a vertex of the mesh
            }
            halfedge_t* h = nullptr;
            if(find_halfedge(n,c)){
                    h = find_halfedge(n,c);}
            else{
                    h = emplace_halfedge_(n);
                    //halfedges_.emplace_back(n_halfedges_++, n);
                    //h = std::addressof(halfedges_.back());  
            }
            //NON STIAMO ASSEGNANDO L'HALFEDGE AL NODO 
            ghost_halfedges[i+2] = h;
        }
        // add edges
        std::vector<std::pair<coords_t, coords_t>> boundary_edges = get_boundary_edges();
        for (int i = 0; i < nodes_polygon+2 ; ++i) {
            halfedge_t* h1 = ghost_halfedges[i];
            halfedge_t* h2 = ghost_halfedges[(i + 1) % (nodes_polygon + 2)];

            // Check for intersection with boundary edges
            coords_t A = h1->node()->coords();
            coords_t B = h2->node()->coords();

            bool flag=false;
            if(!(h1->on_boundary() && h2->on_boundary())){  //if the edges are not both on the boundary
               int num = 0;
                for (const auto& edge : boundary_edges) {
                    if(num==boundary_edges.size()/2) break;   //DA RIVEDERE IL CHECK DOPO AVER AGGIUNTO COME ATTRIBUTO I BOUNDARY EDGES
                    num++;
                    if(A!=edge.first && A!=edge.second && B!=edge.first && B!=edge.second){
                      if (do_segments_intersect(A, B, edge.first, edge.second)) {
                        std::cout << "Error:intersection with boundary!" << std::endl;
                        flag=true;
                      }
                    }
                }
            }
            if(!flag){
                ghost_halfedges[(i + 1) % (nodes_polygon + 2)] = insert_edge(h1, h2)->next();
            } 
            else{ 
                ghost_halfedges[(i + 1) % (nodes_polygon + 2)] = insert_edge(h1->prev(), h1->next())->next();
            }
        }

        c->set_halfedge(v);

        return c->halfedge();
    }

    node_t* find_node(const Eigen::Matrix<double, Eigen::Dynamic, 1>& position) {
        for (auto& n : nodes_) 
            if (n.coords().isApprox(position, 1e-6)) 
                return &n;
        return nullptr;
    }
 
    halfedge_t* find_halfedge(node_t* n, cell_t* cell) {
        for (auto& h : halfedges_) {
            if (h.cell() == cell && h.node() == n) 
                return &h;
        }
        return nullptr;
    }


    bool do_segments_intersect(const coords_t& A, const coords_t& B, const coords_t& C, const coords_t& D) const {
            // Funzione di orientazione: restituisce il segno dell'area del parallelogramma formato dai tre punti
            auto orientation = [](const coords_t& P, const coords_t& Q, const coords_t& R) -> double {
                return (Q.x() - P.x()) * (R.y() - P.y()) - (Q.y() - P.y()) * (R.x() - P.x());
            };

            double o1 = orientation(A, B, C);
            double o2 = orientation(A, B, D);
            double o3 = orientation(C, D, A);
            double o4 = orientation(C, D, B);
        
            //general case: the 2 segments intersect properly
            if ((o1 * o2 < 0) && (o3 * o4 < 0)) {
                return true;
            }
            //degenerate cases: points are colinear and one is inside the other segment
            auto on_segment = [](const coords_t& P, const coords_t& Q, const coords_t& R) -> bool {
                return (std::min(P.x(), Q.x()) <= R.x() && R.x() <= std::max(P.x(), Q.x()) &&
                        std::min(P.y(), Q.y()) <= R.y() && R.y() <= std::max(P.y(), Q.y()));
            };
            if (o1 == 0 && on_segment(A, B, C)) return true;
            if (o2 == 0 && on_segment(A, B, D)) return true;
            if (o3 == 0 && on_segment(C, D, A)) return true;
            if (o4 == 0 && on_segment(C, D, B)) return true;
        
            return false;
    }
        
      
    //Boundary edges   // DA TOGLIERE
    std::vector<std::pair<coords_t, coords_t>> get_boundary_edges() const {
            std::vector<std::pair<coords_t, coords_t>> boundary_edges;
            for (const auto& h : halfedges_) {
                if (h.on_boundary() && h.twin()->on_boundary()) {
                    boundary_edges.emplace_back(h.node()->coords(), h.twin()->node()->coords());
                    
                }
            }
            return boundary_edges;
    }


        

 
 //////////////////////////////   FINE FUNZIONI NUOVE DCEL ////////////////////////////////////
   private:
    // internal utils
    template <typename... Args> halfedge_t* emplace_halfedge_(Args&&... args) {
        halfedges_.emplace_back(n_halfedges_++, std::forward<Args>(args)...);
        return std::addressof(halfedges_.back());
    }

    template <typename... Args> node_t* emplace_node_(Args&&... args) {
        nodes_.emplace_back(n_nodes_++, std::forward<Args>(args)...);
        return std::addressof(nodes_.back());
    }  

    // internal storage (use list to avoid reallocations)
    std::list<node_t> nodes_;
    std::list<halfedge_t> halfedges_;
    std::list<cell_t> cells_;
    int n_nodes_, n_halfedges_, n_cells_;
};


}   // namespace fdapde



#endif // __DCEL_H__
