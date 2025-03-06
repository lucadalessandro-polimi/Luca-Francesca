      
        void set_on_boundary(bool b) {
            boundary_ = b;
        }

        halfedge_t* cell_on_boundary() const { 
            if (!h_) return nullptr;  
        
            halfedge_t* h1 = h_;
            do {
                if (h1->on_boundary()) return h1;  
                h1 = h1->next(); 
            } while (h1 != h_);  
        
            return nullptr;  
        }


           //////////////////CODICE DI PROVA /////////////////////////////
    //(deve rimuovere tutti i lati interni e creare una cavità)
    void remove_polygon(const cell_t* cell){    //DA CORREGGERE QUEI DO-WHILE
       
        halfedge_t* b = cell->on_boundary();
        if (b){ // only remove edges on boundary    DOVREBBE ARRIVARE FINO A CASO b. IN CUI RIMANE SOLO 1 LATO SUL BORDO
        do{
         b = remove_edge(b);
        }while(b->on_boundary())
        return; 
        }

        halfedge_t* h1 = cell->halfedge();
        halfedge_t* ending = h1->twin()->next();
        do{
        h1 = remove_edge(h1);
        }while(h1!=ending)

    // remove cells --> already done by remove_edge  IN TEORIA, RICONTROLLARE
    }
    /////////////////FINE CODICE DI PROVA //////////////////

        ////////////////////CODICE DI PROVA ////////////////////////
        bool remove_node(const node_t& node) {
            auto it = std::find(nodes_begin(), nodes_end(), node);
            if (it == nodes_end()) return false; 
            
            nodes_.remove(node);
            n_nodes_--;
            return true;
        }
        ///////////////////FINE CODICE DI PROVA ////////////////////

            ///////////////////////////// CODICE DI PROVA///////////////////////////////////////7
    halfedge_t* remove_edge(halfedge_t* v1){
        if (v1->on_boundary()) { // v1 is on the boundary
        //2 cases: 
        //a. v1 and its next/prev are on boundary --> elongate v1->next or prev (SE NO NON PERMETTERE)
        if(v1->next->on_boundary() || v1->prev()->on_boundary()){
            v1->next()->set_prev(v1->prev());
            v1->prev()->set_next(v1->next());
            v1->next()->set_node(v1->node());     
            halfedge_t* v2= v1->twin();
            v2->next()->set_prev(v2->prev());
            v2->prev()->set_next(v2->next());
            v2->next()->set_node(v2->node());
            // remove edges, cell is the same so it doesn't need to be removed
            halfedge_t* next= v1->next();
            n_halfedges_ -= 2;
            auto it1 = std::find(halfedges_begin(), halfedges_end(), v1);
            auto it2 = std::find(halfedges_begin(), halfedges_end(), v2);
            halfedges_.erase(it1);
            halfedges_.erase(it2); 
            return next;
        }
        //b. only v1 is on boundary and its next/prev isn't
        else{
            halfedge_t* v2= v1->twin();
            v1->next()->set_prev(v2->prev());
            v2->prev()->set_next(v1->next());    
            v1->prev()->set_next(v2->next());
            v2->next()->set_prev(v1_prev());
            // remove edges and cell
            n_cells_--;
            auto it = std::find(cells_begin(), cells_end(), v1->cell());
            halfedge_t* end = v1, begin = v1->next();
            do {
               begin->set_cell(nullptr);
               begin->node()->set_on_boundary(true);
               begin = begin->next();
            } while (begin!=end)
            cells_.erase(it);
            halfedge_t* next= v1->next();
            n_halfedges_ -= 2;
            auto it1 = std::find(halfedges_begin(), halfedges_end(), v1);
            auto it2 = std::find(halfedges_begin(), halfedges_end(), v2);
            halfedges_.erase(it1);
            halfedges_.erase(it2); 
            return next; 
        }
    } 

    halfedge_t* v2 = v1->twin();
    // set next of v1_prev to v2_next
    v1->prev()->set_next(v2->next());      
    v2->next()->set_prev(v1->prev());
    // set next of v2_prev to v1_next
    v2->prev()->set_next(v1->next());
    v1->next()->set_prev(v2->prev());
    // remove v2's cell
    n_cells_--;
    auto it = std::find(cells_begin(), cells_end(), v2->cell());
    cells_.erase(it);    
    // insert halfedges in v1's cell
    halfedge_t* end = v2, begin = v2->next();
    cell_t* c1 = v1->cell();
    c1->set_halfedge(begin); // either this or v1->prev()
    do {
        begin->set_cell(c1);
        begin = begin->next();
    } while (begin != end);
    // remove v1 and v2
    halfedge_t* next= v1->next();
    n_halfedges_ -= 2;
    auto it1 = std::find(halfedges_begin(), halfedges_end(), v1);
    auto it2 = std::find(halfedges_begin(), halfedges_end(), v2);
    halfedges_.erase(it1);
    halfedges_.erase(it2);  
    return next;
    }

    /////////////////////////////// FINE CODICE DI PROVA//////////////////////////


        ////////////////////////////// CODICE DI PROVA //////////////////////////////
        node_t* adjacent(halfedge_t* edge) const {
            return (edge->twin()) ? edge->twin()->prev()->node() : nullptr;  //da controllare
        }
        ////////////////////////////// FINE CODICE DI PROVA //////////////////////////

/////////////GIA IMPLEMENTATATA DA LORO IN primitives.h///////////////////////////////////7777
        bool is_point_inside_triangle(const node_t& node_P, const node_t& node_A, const node_t& node_B, const node_t& node_C) const{
            // baricenter test
            coords_t P= node_P.coords();
            coords_t A= node_A.coords();
            coords_t B= node_B.coords();
            coords_t C= node_C.coords();
            Eigen::Matrix<double, LocalDim, LocalDim> M;
            M << (B - A), (C - A); 
            coords_t lambda = M.inverse() * (P - A);
            double lambda1 = lambda.x();
            double lambda2 = lambda.y();
            double lambda3 = 1 - lambda1 - lambda2;
            return (lambda1 >= 0 && lambda1<=1 && lambda2 >= 0 && lambda2<=1 && lambda3 >= 0 && lambda3<=1);  
        }

//////////////////////////////// fine ////////////////////////////////////////////

/////////////////////////////// riadattata in primiteves.h //////////////////////
bool in_circle(const node_t& node_A, const node_t& node_B, const node_t& node_C, const node_t& node_D) const{
        
    coords_t A= node_A.coords();
    coords_t B= node_B.coords();
    coords_t C= node_C.coords();
    coords_t D= node_D.coords();
    
    Eigen::Matrix<double, LocalDim+1, LocalDim+1> M;
    M << (A-D).x(), (A-D).y(), (A-D).squaredNorm(), (B-D).x(), (B-D).y(), (B-D).squaredNorm(), (C-D).x(), (C-D).y(), (C-D).squaredNorm();
    double det = M.determinant();
    // if det is positive, D is inside A-B-C circumcircle
    return det > 0; 
} 
////////////////////////////// fine //////////////////////////////////////////77