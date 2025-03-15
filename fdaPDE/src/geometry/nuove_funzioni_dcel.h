      
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
            return (edge->twin()) ? edge->twin()->prev()->node() : nullptr;  
        }
        ////////////////////////////// FINE CODICE DI PROVA //////////////////////////







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



    // dcel.remove_polygon(find_halfedge_from_id(20)->cell());
    /*
    DCEL<2, 2>::node_t node1(dcel.n_nodes(),nullptr,false,0.5,1.), node2(dcel.n_nodes()+1,nullptr,false,1.5,1.);
    dcel.insert_node(node1);
    dcel.insert_node(node2);
    DCEL<2, 2>::halfedge_t h_i1(dcel.n_halfedges()+1000,&node1);
    DCEL<2, 2>::halfedge_t h_i2(dcel.n_halfedges()+1001,&node2);
    h_i1.set_cell(nullptr);
    h_i2.set_cell(nullptr);  //SE NO BISOGNA MODIFICARE IL COSTRUTTORE
    node1.set_halfedge(&h_i1);
    node2.set_halfedge(&h_i2);
    DCEL<2, 2>::halfedge_t* he1 = find_halfedge_from_node(1);
    DCEL<2, 2>::halfedge_t* hi1=&h_i1;
    DCEL<2, 2>::halfedge_t* hi2=&h_i2;

    //TEST insert e remove edge
    DCEL<2, 2>::halfedge_t* h10= dcel.insert_edge(he1, hi1);
    DCEL<2, 2>::halfedge_t* h12= dcel.insert_edge(hi2, he1);
    DCEL<2, 2>::halfedge_t* h3=dcel.insert_edge(h10->twin(), h12);
    dcel.insert_edge(find_halfedge_from_id(3), find_halfedge_from_id(12));
    dcel.insert_edge(find_halfedge_from_id(14), find_halfedge_from_id(4));
    dcel.insert_edge(find_halfedge_from_id(14), find_halfedge_from_id(3));
    dcel.insert_edge(find_halfedge_from_id(0), find_halfedge_from_id(18));   
    dcel.remove_edge(find_halfedge_from_id(10));
    dcel.remove_edge(find_halfedge_from_id(13));
    dcel.remove_edge(find_halfedge_from_id(18));
    dcel.remove_edge(find_halfedge_from_id(21));
    dcel.remove_edge(find_halfedge_from_id(16));
    dcel.remove_edge(find_halfedge_from_id(14));
    */







    halfedge_t* remove_edge(halfedge_t* v1){
        if(!v1) return nullptr;
        /*
        if (v1->on_boundary()) { // v1 is on the boundary
        //2 cases: 
        //a. v1 and its next/prev are on boundary --> elongate v1->next or prev (SE NO NON PERMETTERE)
        if(v1->next()->on_boundary() || v1->prev()->on_boundary()){
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
            v2->next()->set_prev(v1->prev());
            // remove edges and cell
            n_cells_--;
            auto it = std::find(cells_begin(), cells_end(), v1->cell());
            halfedge_t* end = v1;
            halfedge_t* begin = v1->next();
            do {
               begin->set_cell(nullptr);
               begin->node()->set_on_boundary(true);
               begin = begin->next();
            } while (begin!=end);
            cells_.erase(it);
            halfedge_t* next= v1->next();
            n_halfedges_ -= 2;
            auto it1 = std::find(halfedges_begin(), halfedges_end(), v1);
            auto it2 = std::find(halfedges_begin(), halfedges_end(), v2);
            halfedges_.erase(it1);
            halfedges_.erase(it2); 
            return next; 
        }
        } */
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
        n_cells_--;
        cell_t* c2=v2->cell();
        if (c2) {
            //auto it = std::find(cells_begin(), cells_end(), v2->cell());
            //auto it = std::find_if(cells_.begin(), cells_.end(), [=](const cell_t& c) { return c.id() == c2->id(); });
            cell_iterator it = cells_.begin();
            while (it != cells_.end()) {
                if (it->id() == c2->id()) 
                    break;  // Fermiamo il loop se troviamo la cella
                ++it;
            }
            if (it != cells_.end()) 
                cells_.erase(it);
        }
        // remove v1 and v2
        halfedge_t* next= v1->next();
        n_halfedges_ -= 2;
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

    halfedge_t* add_polygon(halfedge_t* v, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& nodes){
        // update n_halfedges_ and n_cells_
        // n_nodes_ is already updated by insert_node
        int nodes_polygon= nodes.rows();
        cell_t* c= v->cell();

        std::vector<halfedge_t*> ghost_halfedges(nodes_polygon +2 ); //O(n)
        ghost_halfedges[0] = v;
        ghost_halfedges[1] = v->next();

        // add nodes and create ghost halfedges
        for (int i = 0; i < nodes_polygon; ++i) {
            node_t* n = insert_node(node_t(n_nodes_, /* boundary = */ false, nodes.row(i)));
            halfedges_.emplace_back(n_halfedges_+ 1000 + i, n);
            halfedge_t* h = std::addressof(halfedges_.back());    //NON USARE emplace_halfedge PERCHè INCASINA GLI INDICI, così li posso controllare io!
            ghost_halfedges[i+2] = h;
        }
        // add edges
        int count=0;
        for (int i = 0; i < nodes_polygon+2 ; ++i) {
            halfedge_t* h1 = ghost_halfedges[i];
            halfedge_t* h2 = ghost_halfedges[(i + 1) % (nodes_polygon+2)];  
            ghost_halfedges[(i + 1) % (nodes_polygon+2)]= insert_edge(h1, h2)->next();  
            if(!h2->next()) {
                auto it = std::find_if(halfedges_.begin(), halfedges_.end(), [=](const halfedge_t& h) { return h.id() == h2->id();});
                if (it != halfedges_.end()){
                    halfedges_.erase(it);
                    count++;
                } 
            } 
        }

        return c->halfedge();
    }




    //   void remove_polygon(const cell_t* cell){    //PRENDE IN INGRESSO CELLA O HALFEDGE?
       
        /*halfedge_t* b = cell->on_boundary();
        if (b){ // only remove edges on boundary    DOVREBBE ARRIVARE FINO A CASO b. IN CUI RIMANE SOLO 1 LATO SUL BORDO
        do{
         b = remove_edge(b);
        }while(b->on_boundary())
        return; 
        }*/
/*
        halfedge_t* h1 = cell->halfedge();
        halfedge_t* ending = h1->twin()->next();
        do{
            h1 = remove_edge(h1);
        }while(h1!=ending);
    }
*/





void add_triangle(const coords_t& A, const coords_t& B, const coords_t& C) {
    std::cout << "🔺 Aggiungo triangolo con vertici: " << A.transpose() << ", " 
              << B.transpose() << ", " << C.transpose() << std::endl;
    // ⚠️ Controllo per evitare triangoli degeneri
    if (A == B || B == C || C == A) {
        std::cerr << "❌ Errore: triangolo degenere con vertici coincidenti!" << std::endl;
        return;
    }


    // 1️⃣ Cerchiamo i nodi esistenti
    node_t* nA = find_node(A);
    node_t* nB = find_node(B);
    node_t* nC = find_node(C);

    if (!nA) nA = insert_node(node_t(n_nodes_, false, A));
    if (!nB) nB = insert_node(node_t(n_nodes_, false, B));
    if (!nC) nC = insert_node(node_t(n_nodes_, false, C));

    std::cout << "Nodi: A=" << nA->id() << " B=" << nB->id() << " C=" << nC->id() << std::endl;

    // 2️⃣ Cerchiamo se gli half-edge esistono già
    halfedge_t* hAB = find_halfedge(nA, nB);
    halfedge_t* hBC = find_halfedge(nB, nC);
    halfedge_t* hCA = find_halfedge(nC, nA);

    std::cout << "Edge trovati: AB=" << (hAB ? hAB->id() : -1) 
              << " BC=" << (hBC ? hBC->id() : -1) 
              << " CA=" << (hCA ? hCA->id() : -1) << std::endl;

    // 3️⃣ Se non esistono, li creiamo
    if (!hAB) {
        hAB = emplace_halfedge_(nA);
        std::cout << "Creato half-edge AB con ID: " << hAB->id() << std::endl;
    }
    if (!hBC) {
        hBC = emplace_halfedge_(nB);
        std::cout << "Creato half-edge BC con ID: " << hBC->id() << std::endl;
    }
    if (!hCA) {
        hCA = emplace_halfedge_(nC);
        std::cout << "Creato half-edge CA con ID: " << hCA->id() << std::endl;
    }

    // **CHECK: Tutti gli half-edge devono esistere**
    if (!hAB || !hBC || !hCA) {
        std::cerr << "❌ Errore: uno degli half-edge è NULL dopo l'inserimento!" << std::endl;
        return;
    }

    // 4️⃣ Creiamo e colleghiamo i twin
    halfedge_t* hBA = find_halfedge(nB, nA);
    halfedge_t* hCB = find_halfedge(nC, nB);
    halfedge_t* hAC = find_halfedge(nA, nC);

    if (!hBA) hBA = emplace_halfedge_(nB);
    if (!hCB) hCB = emplace_halfedge_(nC);
    if (!hAC) hAC = emplace_halfedge_(nA);

    // 5️⃣ Impostiamo i twin
    hAB->set_twin(hBA); hBA->set_twin(hAB);
    hBC->set_twin(hCB); hCB->set_twin(hBC);
    hCA->set_twin(hAC); hAC->set_twin(hCA);

    // 6️⃣ Colleghiamo next e prev
    hAB->set_next(hBC); 
    hBC->set_next(hCA); 
    hCA->set_next(hAB);

    hBA->set_prev(hCB); 
    hCB->set_prev(hAC); 
    hAC->set_prev(hBA);

    hAB->set_prev(hCA); 
    hBC->set_prev(hAB); 
    hCA->set_prev(hBC);

    hBA->set_next(hAC); 
    hCB->set_next(hBA); 
    hAC->set_next(hCB);

    // 🔍 **Check extra per assicurarci che tutti i next siano assegnati correttamente**
    std::cout << "🔍 Debug connessioni next dopo assegnazione: " << std::endl;
    std::cout << "AB: " << hAB->id() << " -> " << (hAB->next() ? std::to_string(hAB->next()->id()) : "nullptr") << std::endl;
    std::cout << "BC: " << hBC->id() << " -> " << (hBC->next() ? std::to_string(hBC->next()->id()) : "nullptr") << std::endl;
    std::cout << "CA: " << hCA->id() << " -> " << (hCA->next() ? std::to_string(hCA->next()->id()) : "nullptr") << std::endl;

    if (!hAB->next() || !hBC->next() || !hCA->next()) {
        std::cerr << "❌ Errore: almeno un half-edge non ha next() assegnato correttamente!" << std::endl;
        return;
    }
    

    // 7️⃣ Creiamo la cella
    cells_.push_back(cell_t(n_cells_++));
    cell_t* triangle = std::addressof(cells_.back());

    // Assegniamo la cella agli edge
    hAB->set_cell(triangle);
    hBC->set_cell(triangle);
    hCA->set_cell(triangle);

    triangle->set_halfedge(hAB);

    std::cout << "✅ Triangolo aggiunto con successo: " << triangle->id() << std::endl;
}

halfedge_t* find_halfedge(node_t* from, node_t* to) {
    if (!from || !to) return nullptr; // Protezione extra

    for (auto it = halfedges_.begin(); it != halfedges_.end(); ++it) {
        if (it->node() == from && it->next() && it->next()->node() == to) {
            std::cout << "✅ Trovato half-edge esistente da " << from->id() << " a " << to->id() << ": " << it->id() << std::endl;
            return std::addressof(*it);
        }
    }

    std::cout << "❌ Nessun half-edge trovato da " << from->id() << " a " << to->id() << std::endl;
    return nullptr;
}



node_t* find_node(const coords_t& coords) {
    for (auto it = nodes_begin(); it != nodes_end(); ++it) {
        if ((it->coords() - coords).norm() < 1e-8) { // Tolleranza per errore numerico
            return std::addressof(*it);
        }
    }
    return nullptr; // Se non trova nulla, restituisce nullptr
}




void add_triangle(const coords_t& p1, const coords_t& p2, const coords_t& p3) {
    // Funzione per trovare un nodo esistente
    auto find_node = [&](const coords_t& coords) -> node_t* {
        for (auto it = nodes_begin(); it != nodes_end(); ++it) {
            if ((it->coords() - coords).norm() < 1e-10) { // Tolleranza per errore numerico
                return std::addressof(*it);
            }
        }
        return nullptr;
    };

    // Creazione o ricerca dei nodi
    node_t* n1 = find_node(p1);
    if (!n1) n1 = insert_node(node_t(n_nodes_, /* boundary = */ false, p1));
    
    node_t* n2 = find_node(p2);
    if (!n2) n2 = insert_node(node_t(n_nodes_, /* boundary = */ false, p2));
    
    node_t* n3 = find_node(p3);
    if (!n3) n3 = insert_node(node_t(n_nodes_, /* boundary = */ false, p3));

    // Funzione per trovare un half-edge esistente
    auto find_halfedge = [&](node_t* from, node_t* to) -> halfedge_t* {
        for (auto it = halfedges_begin(); it != halfedges_end(); ++it) {
            if (it->node() == from && it->next() && it->next()->node() == to) {
                return std::addressof(*it);
            }
        }
        return nullptr;
    };

    // Creazione o ricerca dei half-edge e delle twin
    halfedge_t* h1 = find_halfedge(n1, n2);
    if (!h1) h1 = emplace_halfedge_(n1);
    if (!h1) { std::cerr << "Errore: Creazione h1 fallita!\n"; return; }
    
    halfedge_t* h2 = find_halfedge(n2, n3);
    if (!h2) h2 = emplace_halfedge_(n2);
    if (!h2) { std::cerr << "Errore: Creazione h2 fallita!\n"; return; }
    
    halfedge_t* h3 = find_halfedge(n3, n1);
    if (!h3) h3 = emplace_halfedge_(n3);
    if (!h3) { std::cerr << "Errore: Creazione h3 fallita!\n"; return; }

    halfedge_t* h1_twin = find_halfedge(n2, n1);
    if (!h1_twin) h1_twin = emplace_halfedge_(n2);
    
    halfedge_t* h2_twin = find_halfedge(n3, n2);
    if (!h2_twin) h2_twin = emplace_halfedge_(n3);
    
    halfedge_t* h3_twin = find_halfedge(n1, n3);
    if (!h3_twin) h3_twin = emplace_halfedge_(n1);

    // Controlli sui twin
    if (h1 && h1_twin) { h1->set_twin(h1_twin); h1_twin->set_twin(h1); }
    if (h2 && h2_twin) { h2->set_twin(h2_twin); h2_twin->set_twin(h2); }
    if (h3 && h3_twin) { h3->set_twin(h3_twin); h3_twin->set_twin(h3); }

    // Collegamento delle strutture next-prev
    if (h1 && h2) { h1->set_next(h2); h2->set_prev(h1); }
    if (h2 && h3) { h2->set_next(h3); h3->set_prev(h2); }
    if (h3 && h1) { h3->set_next(h1); h1->set_prev(h3); }

    if (h1_twin && h3_twin) { h1_twin->set_next(h3_twin); h3_twin->set_prev(h1_twin); }
    if (h3_twin && h2_twin) { h3_twin->set_next(h2_twin); h2_twin->set_prev(h3_twin); }
    if (h2_twin && h1_twin) { h2_twin->set_next(h1_twin); h1_twin->set_prev(h2_twin); }

    // Creazione della cella per il triangolo
    cells_.push_back(cell_t(n_cells_++));
    cell_t* c = std::addressof(cells_.back());
    c->set_halfedge(h1);

    // Associazione cella agli half-edge
    if (h1) h1->set_cell(c);
    if (h2) h2->set_cell(c);
    if (h3) h3->set_cell(c);

    // Associa i nodi ai loro half-edge
    if (n1) n1->set_halfedge(h1);
    if (n2) n2->set_halfedge(h2);
    if (n3) n3->set_halfedge(h3);
}




void remove_triangle(const cell_t* triangle) {
    dcel_.remove_polygon(triangle);  // Riutilizziamo remove_polygon
}


halfedge_t* add_polygon(halfedge_t* v, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& nodes){
    // update n_halfedges_ and n_cells_
    // n_nodes_ is already updated by insert_node
    int nodes_polygon= nodes.rows();
    cell_t* c= v->cell();

    std::vector<halfedge_t*> ghost_halfedges(nodes_polygon +2 ); //O(n)
    ghost_halfedges[0] = v;
    ghost_halfedges[1] = v->next();

    // add nodes and create ghost halfedges
   for (int i = 0; i < nodes_polygon; ++i) {
       node_t* n = insert_node(node_t(n_nodes_, /* boundary = */ false, nodes.row(i)));
   // add nodes and create ghost halfedges

        halfedges_.emplace_back(n_halfedges_+ 1000 + i, n);
        halfedge_t* h = std::addressof(halfedges_.back());    //NON USARE emplace_halfedge PERCHè INCASINA GLI INDICI, così li posso controllare io!
        ghost_halfedges[i+2] = h;
    }
    // add edges
    int count=0;
    for (int i = 0; i < nodes_polygon+2 ; ++i) {
        halfedge_t* h1 = ghost_halfedges[i];
        halfedge_t* h2 = ghost_halfedges[(i + 1) % (nodes_polygon+2)];  
        ghost_halfedges[(i + 1) % (nodes_polygon+2)]= insert_edge(h1, h2)->next();  
        if(!h2->next()) {
            auto it = std::find_if(halfedges_.begin(), halfedges_.end(), [=](const halfedge_t& h) { return h.id() == h2->id();});
            if (it != halfedges_.end()){
                halfedges_.erase(it);
                count++;
            } 
        } 
    }

    return c->halfedge();
}

 /*   halfedge_t* insert_edge(halfedge_t* v1, halfedge_t* v2) {
        // AGGIUNGO CONTROLLO ULTERIORE  DENTRO ALL'IF per essere sicura che v1,v2 e v2->next() non siano null
        if (v1->cell() && v2-> cell() && v1->cell()!=v2->cell()){
            std::cout << "Errore: i due half-edge appartengono a celle diverse. Non è possibile inserire un edge tra i due" << std::endl;
            return nullptr;
        }
        if (v1 && v2 && v2->next() && v1->next() && ( v1->node() == v2->next()->node() || v2->node()==v1->next()->node()) ) {
            std::cout << "Halfedges consecutivi" << std::endl;
            return v1;   // v1 and v2 are next halfedges
        }
        // get exiting halfedges from n1 and n2
        node_t* n1 = v1->node();
        node_t* n2 = v2->node();
        // create a pair of twin half-edges
        halfedge_t* h1 = emplace_halfedge_(n1);
        halfedge_t* h2 = emplace_halfedge_(n2);
        h1->set_twin(h2);
        h2->set_twin(h1);
	    // insert halfedge h1 between v1 and v1->prev
        h2->set_next(v1);
        // AGGIUNGO CONTROLLO SU v1->prev()
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
	    // insert halfedge h2 between v2 and v2->prev
        h1->set_next(v2);
        // AGGIUNGO CONTROLLO SU v2->prev()
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
	    // set cell pointers
        // create new cell
        // SOLO SE L'EDGE NUOVO HA PREV E NEXT DIVERSI DA Sè STESSO
        if(h1->next() != h2 && h1->prev() != h2){
            //h1->set_cell(v1->cell());  //POTENZIALI PROBLEMI
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
    */


    /*
   // delaunay.print_dcel();
    delaunay.dcel().add_polygon(find_halfedge_from_id(0),internal.row(0));
    
    delaunay.dcel().add_polygon(find_halfedge_from_id(1),internal.row(0));
    //delaunay.print_dcel();
    delaunay.dcel().add_polygon(find_halfedge_from_id(2),internal.row(0));
    //delaunay.print_dcel();
    delaunay.dcel().add_polygon(find_halfedge_from_id(3),internal.row(0));

    Eigen::Matrix<double, 2, 2> prova;
    prova << 1.5, 0.5,
             0.5, 0.5;
    delaunay.dcel().add_polygon(find_halfedge_from_id(3),prova);
  
    delaunay.print_dcel();*/

        /*DCEL<2,2>::coords_t u = internal.row(1);
    std::cout << "🔍 Inserimento del punto interno: " << u.transpose()<< std::endl;
    
    const DCEL<2,2>::cell_t* triangle = delaunay.find_triangle(u);
    
    delaunay.insert_vertex(u,triangle);

    DCEL<2,2>::coords_t u_ = internal.row(2);
    std::cout << "🔍 Inserimento del punto interno: " << u_.transpose()<< std::endl;
    
    const DCEL<2,2>::cell_t* triangle_ = delaunay.find_triangle(u_);
    
    delaunay.insert_vertex(u_,triangle_);*/