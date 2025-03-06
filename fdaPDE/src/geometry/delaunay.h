#ifndef __FDAPDE_DELAUNAY_H__
#define __FDAPDE_DELAUNAY_H__

#include "header_check.h"
#include "dcel.h"

namespace fdapde {
  
template <int LocalDim = 2, int EmbedDim = 2>
class Delaunay {
   public:
    static constexpr int local_dim = LocalDim;
    static constexpr int embed_dim = EmbedDim;

    using coords_t = Eigen::Matrix<double, embed_dim, 1>;
    using node_t = typename DCEL<local_dim, embed_dim>::node_t;
    using halfedge_t = typename DCEL<local_dim, embed_dim>::halfedge_t;
    using cell_t = typename DCEL<local_dim, embed_dim>::cell_t;

    // Costruttori

    //

    cell_t* find_triangle(const node_t& P) const {
        for (auto it = dcel_.cells_begin(); it != dcel_.cells_end(); ++it) {
            cell_t* cell = &(*it);

            const node_t& A = cell->halfedge()->node();
            const node_t& B = cell->halfedge()->next()->node();
            const node_t& C = cell->halfedge()->prev()->node();

            if (is_point_inside_triangle(P, A, B, C)) {
                return cell;
            }
        }
        return nullptr; // Il punto non è contenuto in nessun triangolo
    }

    void add_triangle(node_t* u, node_t* v, node_t* w) {

        Eigen::Matrix<double, 3, embed_dim> nodes;
        nodes.row(0) = u->coords();
        nodes.row(1) = v->coords();
        nodes.row(2) = w->coords();

        DCEL<local_dim, embed_dim> triangle = DCEL<local_dim, embed_dim>::make_polygon(nodes);
      //MANCA TUTTA LA PARTE SUL COME INSERIRE NELLA STRUTTURA DATI GIA ESISTENTE 
    }
    

    // Funzione per eliminare un triangolo dalla DCEL
    void delete_triangle(cell_t* triangle) {
        if (!triangle) return;  
        dcel_.remove_polygon(triangle);
    }
   

    void dig_cavity(node_t* u, halfedge_t* vw) { 
        
        node_t* x = dcel_.adjacent(vw);
        if (!x) return;

        if (in_circle(*u, *(vw->node()), *(vw->next()->node()), *x)) {  
        //itero ragionamento fatto nella insert_vertex 
        // si possono evitare i tre sotto e scrivere tutto esplicito 
          halfedge_t* wv = vw->twin();
          halfedge_t* vx = vw->twin()->next();
          halfedge_t* xw = vw->twin()->prev(); 
          
          bool wv_is_boundary = wv->on_boundary();
          bool vx_is_boundary = vx->on_boundary();
          bool xw_is_boundary = xw->on_boundary();

          delete_triangle(vw->twin()->cell());
          
          if (!vx_is_boundary) dig_cavity(u, vx);
          if (!xw_is_boundary) dig_cavity(u, xw);
        } 
        else {
          add_triangle(u, vw->node(), vw->next()->node());
        } 
    }
//il find_triangle mi dice che u sta in tringle ---> concetto da esprimere in funzione che fa la triangolazione 
    void insert_vertex(node_t* u,cell_t* triangle) {
     // si possono evitare i tre sotto e scrivere tutto esplicito
        halfedge_t* vw = triangle->halfedge();
        halfedge_t* wx = vw->next();
        halfedge_t* xv = vw->prev();
    // da ricontrollare se il bordo è gestito correttamente: idea è se un lato è sul bordo poi il delete_triangle lo cancella 
    // quindi non posso fare la dig_cavity li 
        bool vw_is_boundary = vw->on_boundary();
        bool wx_is_boundary = wx->on_boundary();
        bool xv_is_boundary = xv->on_boundary();

        delete_triangle(triangle);
        
        if (!vw_is_boundary) dig_cavity(u, vw);
        if (!wx_is_boundary) dig_cavity(u, wx);
        if (!xv_is_boundary) dig_cavity(u, xv);
    }
    
    // Getter
    const DCEL<local_dim, embed_dim>& dcel() const { return dcel_; }
   
   private:
    DCEL<local_dim, embed_dim> dcel_;  
};
  
}   // namespace fdapde

#endif // __DELAUNAY_H__