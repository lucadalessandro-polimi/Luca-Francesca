#ifndef __FDAPDE_DCEL_H__
#define __FDAPDE_DCEL_H__

#include "header_check.h"
#include "dcel.h"

namespace fdapde {





    //DA GENERALIZZARE A N-DIM (TUTTO)
    bool is_point_inside_triangle(const Eigen::Vector2d& P, const Eigen::Vector2d& A, 
        const Eigen::Vector2d& B, const Eigen::Vector2d& C) {
// Calcoliamo le coordinate baricentriche
Eigen::Matrix2d M;
M << (B - A), (C - A); //CAPIRE IL RIEMPIEMENTO
Eigen::Vector2d lambda = M.inverse() * (P - A); //COSTO COMPUTAZIONE CON INVERSO RISPETTO AL CILCO FOR 
double lambda1 = lambda.x();
double lambda2 = lambda.y();
double lambda3 = 1 - lambda1 - lambda2;

return (lambda1 >= 0 && lambda2 >= 0 && lambda3 >= 0);  //DEVE STARE TRA ZERO E UNO I COEFFICIENTI
}

// GLI PASSI UN FLAG PER SCARTARE LE DIMENSIONI MAGGIORI DI 3?
// Funzione per verificare se un punto D è dentro il circumcerchio del triangolo ABC
bool InCircle(const Eigen::Vector2d& A, const Eigen::Vector2d& B, const Eigen::Vector2d& C, const Eigen::Vector2d& D) {
Eigen::Matrix4d M;
M << A.x(), A.y(), A.squaredNorm(), 1,
B.x(), B.y(), B.squaredNorm(), 1,
C.x(), C.y(), C.squaredNorm(), 1,
D.x(), D.y(), D.squaredNorm(), 1;

// Calcoliamo il determinante della matrice
double det = M.determinant();

return det > 0; // Se positivo, il punto D è dentro il circumcerchio
}

// Funzione per trovare il triangolo adiacente dall'altro lato di un edge
DCEL<2,2>::cell_t* Adjacent(DCEL<2,2>::halfedge_t* edge) {
// Se l'half-edge ha un twin, il triangolo adiacente è quello del twin
if (edge->twin()) {
return edge->twin()->cell();
}
return nullptr; // Se non ha twin, significa che il lato è sul bordo
}

// Funzione per eliminare un triangolo dalla DCEL
void DeleteTriangle(DCEL<2,2>& dcel, DCEL<2,2>::cell_t* triangle) {
if (!triangle) return; // Se il triangolo è nullo, non facciamo nulla

// Otteniamo gli half-edge del triangolo
DCEL<2,2>::halfedge_t* h1 = triangle->halfedge();
DCEL<2,2>::halfedge_t* h2 = h1->next();
DCEL<2,2>::halfedge_t* h3 = h2->next();

// DA IMPLEMENTARE 
dcel.remove_cell(triangle);

// DA IMPLEMENTARE REMOVE_HALFEDGE
if (h1->twin() == nullptr) dcel.remove_edge(h1);
if (h2->twin() == nullptr) dcel.remove_edge(h2);
if (h3->twin() == nullptr) dcel.remove_edge(h3);
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