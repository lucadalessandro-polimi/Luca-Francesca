#ifndef __FDAPDE_DCEL_H__
#define __FDAPDE_DCEL_H__

#include "header_check.h"
#include "dcel.h"

namespace fdapde {

using LocalDim = 2;
using EmbedDim = 2; 
using coords_t = Eigen::Matrix<double, LocalDim, 1>;
using node_t = DCEL<LocalDim, EmbedDim>::node_t;
using halfedge_t = DCEL<LocalDim, EmbedDim>::halfedge_t;
using cell_t = DCEL<LocalDim, EmbedDim>::cell_t;

bool is_point_inside_triangle(const node_t& node_P, const node_t& node_A, const node_t& node_B, const node_t& node_C)}{
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

bool InCircle(const node_t& node_A, const node_t& node_B, const node_t& node_C, const node_t& node_D){
    
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



// Funzione per trovare il triangolo adiacente dall'altro lato di un edge
DCEL<2,2>::cell_t* Adjacent(DCEL<2,2>::halfedge_t* edge) {
// Se l'half-edge ha un twin, il triangolo adiacente è quello del twin
if (edge->twin()) {
return edge->twin()->cell();
}
return nullptr; // Se non ha twin, significa che il lato è sul bordo
}


//Da inserire in DCEL come duale di insert_edge (non specifico per triangoli ma per poligoni --> serve anche make_triangle?)
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
            auto it1 = std::find(halfedges_.begin(), halfedges_.end(), v1);
            auto it2 = std::find(halfedges_.begin(), halfedges_.end(), v2);
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
            auto it = std::find(cells_.begin(), cells_.end(), v1->cell());
            cells_.erase(it);
            halfedge_t* next= v1->next();
            n_halfedges_ -= 2;
            auto it1 = std::find(halfedges_.begin(), halfedges_.end(), v1);
            auto it2 = std::find(halfedges_.begin(), halfedges_.end(), v2);
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
    auto it = std::find(cells_.begin(), cells_.end(), v2->cell());
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
    auto it1 = std::find(halfedges_.begin(), halfedges_.end(), v1);
    auto it2 = std::find(halfedges_.begin(), halfedges_.end(), v2);
    halfedges_.erase(it1);
    halfedges_.erase(it2);  
    return next;
}

halfedge_t* cell_on_boundary(cell_t* cell){  //metodo di cell_t
    halfedge_t* h1 = cell->halfedge();
    do{
        if(h1->on_boundary()) return h1;
        h1 = h1->next();
    }while(h1!=cell->halfedge())
    return nullptr;
}

// va inserita in DCEL come duale di make_polygon (deve rimuovere tutti i lati interni e creare una cavità)
void remove_polygon(const cell_t* cell){ 

    
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



// Funzione per eliminare un triangolo dalla DCEL
void delete_triangle(DCEL<LocalDim,EmbedDim>& dcel, cell_t* triangle) {
if (!triangle) return; // 

dcel.remove_polygon(triangle);

}


// Walk Search ITERANDO SU TUTTI I TRIANGOLI FINCHE NON TROVI QUELLO CHE CONTIENE IL PUNTO 
DCEL<2,2>::cell_t* FindTriangle(DCEL<2,2>& dcel, const Eigen::Vector2d& P) {
// Partiamo da un triangolo iniziale (scegliamone uno qualsiasi, ad esempio il primo)
DCEL<2,2>::cell_t* current_cell = &(*dcel.cells_begin());

while (true) {
// Otteniamo i vertici del triangolo corrente
Eigen::Vector2d A = current_cell->halfedge()->node()->coords();
Eigen::Vector2d B = current_cell->halfedge()->next()->node()->coords();
Eigen::Vector2d C = current_cell->halfedge()->prev()->node()->coords();

// Se il punto è dentro, lo abbiamo trovato!
if (is_point_inside_triangle(P, A, B, C)) {
return current_cell;
}

// Se non è dentro, troviamo il lato più vicino e ci spostiamo nel triangolo adiacente
DCEL<2,2>::halfedge_t* outgoing_edge = nullptr;

// Controlliamo il lato AB usando il prodotto scalare
// Se il vettore (P - A) punta "all'indietro" rispetto ad AB, significa che P è oltre AB
if ((B - A).dot(P - A) < 0) {
outgoing_edge = current_cell->halfedge(); // Lato AB
}
// Controlliamo il lato BC
// Se il vettore (P - B) punta "all'indietro" rispetto a BC, significa che P è oltre BC
else if ((C - B).dot(P - B) < 0) {
outgoing_edge = current_cell->halfedge()->next(); // Lato BC
}
// Controlliamo il lato CA
// Se il vettore (P - C) punta "all'indietro" rispetto a CA, significa che P è oltre CA
else if ((A - C).dot(P - C) < 0) {
outgoing_edge = current_cell->halfedge()->prev(); // Lato CA
}

// Se non troviamo un lato oltre cui andare, restituiamo il triangolo corrente (errore nella struttura dati?)
if (!outgoing_edge || !outgoing_edge->twin()) {
return current_cell;
}

// Passiamo al triangolo adiacente
current_cell = outgoing_edge->twin()->cell();
}
}


// Funzione per identificare la cavità e rimuovere i triangoli non Delaunay, quindi aggiungere nuovi triangoli
void DigCavity(DCEL<2,2>& dcel, DCEL<2,2>::node_t* u, DCEL<2,2>::halfedge_t* vw) {
DCEL<2,2>::cell_t* wvx = Adjacent(vw);
if (!wvx) return;

DCEL<2,2>::node_t* x = wvx->halfedge()->next()->node();

if (InCircle(u->coords(), vw->node()->coords(), vw->next()->node()->coords(), x->coords())) {
DeleteTriangle(dcel, wvx);
DigCavity(dcel, u, vw->next());
DigCavity(dcel, u, vw->prev());
} else {
// Aggiungiamo il triangolo manualmente utilizzando la logica di make_polygon
DCEL<2,2>::halfedge_t* h1 = dcel.insert_edge(u->halfedge(), vw->node()->halfedge());
DCEL<2,2>::halfedge_t* h2 = dcel.insert_edge(u->halfedge(), vw->next()->node()->halfedge());
DCEL<2,2>::halfedge_t* h3 = dcel.insert_edge(vw->node()->halfedge(), vw->next()->node()->halfedge());

dcel.insert_cell(DCEL<2,2>::cell_t(dcel.n_cells(), h1));
dcel.insert_cell(DCEL<2,2>::cell_t(dcel.n_cells(), h2));
dcel.insert_cell(DCEL<2,2>::cell_t(dcel.n_cells(), h3));
}    //DA VERIFICARE PERCHE DI SICURO BISOGNA USARE MAKE_POLYGON (DA CAPIRE SE è IN GRADO DI COSTRUIRE LA MESH IN MODO INCREMENTALE)
}

// Funzione per inserire un nuovo vertice nella triangolazione
void InsertVertex(DCEL<2,2>& dcel, DCEL<2,2>::node_t* u, DCEL<2,2>::cell_t* triangle) {
if (!triangle) return;

// Otteniamo gli half-edges del triangolo
DCEL<2,2>::halfedge_t* vw = triangle->halfedge();
DCEL<2,2>::halfedge_t* wx = vw->next();
DCEL<2,2>::halfedge_t* xv = wx->next();

// Rimuoviamo il triangolo esistente
DeleteTriangle(dcel, triangle);

// Espandiamo la cavità e creiamo nuovi triangoli
DigCavity(dcel, u, vw);
DigCavity(dcel, u, wx);
DigCavity(dcel, u, xv);
}

// Funzione per costruire una triangolazione di Delaunay da un insieme di punti
DCEL<2,2> DelaunayTriangulation(const std::vector<Eigen::Vector2d>& points) {
DCEL<2,2> dcel;

// Creiamo un maxi-triangolo abbastanza grande da contenere tutti i punti
double min_x = std::numeric_limits<double>::max();
double max_x = std::numeric_limits<double>::lowest();
double min_y = std::numeric_limits<double>::max();
double max_y = std::numeric_limits<double>::lowest();

for (const auto& p : points) {
min_x = std::min(min_x, p.x());
max_x = std::max(max_x, p.x());
min_y = std::min(min_y, p.y());
max_y = std::max(max_y, p.y());
}

double dx = max_x - min_x;
double dy = max_y - min_y;
double delta = std::max(dx, dy) * 10; // Un margine sufficiente

Eigen::Vector2d p1(min_x - delta, min_y - delta);
Eigen::Vector2d p2(max_x + delta, min_y - delta);
Eigen::Vector2d p3((min_x + max_x) / 2, max_y + delta);

// Inseriamo il maxi-triangolo nella DCEL
std::vector<Eigen::Vector2d> super_triangle = {p1, p2, p3};
dcel = DCEL<2,2>::make_polygon(super_triangle);

// Inseriamo i punti uno alla volta
for (const auto& p : points) {
DCEL<2,2>::cell_t* containing_triangle = FindTriangle(dcel, p);
if (!containing_triangle) continue; // Se non troviamo un triangolo valido, saltiamo il punto

DCEL<2,2>::node_t* new_node = dcel.insert_node(DCEL<2,2>::node_t(dcel.n_nodes(), false, p));
InsertVertex(dcel, new_node, containing_triangle);
}

// Rimuoviamo i vertici del maxi-triangolo
for (auto it = dcel.nodes_begin(); it != dcel.nodes_end();) {
if (it->coords() == p1 || it->coords() == p2 || it->coords() == p3) {
it = dcel.nodes_.erase(it);
} else {
++it;
}
}

return dcel;
}




}   // namespace fdapde

#endif // __DCEL_H__