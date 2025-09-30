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

#ifndef __FDAPDE_ADAPTIVITY_COSTS_H__
#define __FDAPDE_ADAPTIVITY_COSTS_H__

#include "header_check.h"

namespace fdapde {

// cost object using the Curiously Recurring Template Pattern (CRTP)
/*template<int LocalDim, int EmbedDim, typename ConcreteCost>
struct CostObjBase{
	double max_ = std::numeric_limits<double>::lowest();
	double min_ = std::numeric_limits<double>::max();
	double threshold_ = 1.5;

	void setup_(){}

	template<typename... CostArgs>
	double operator()(CostArgs&&... cost_args)
	{  	
		double cost = static_cast<ConcreteCost&>(*this).get_cost(std::forward<CostArgs>(cost_args)...);
		if(cost < min_) {min_ = cost;}
		double ret = cost/max_;
		return ret;
	}

	template<typename... CostArgs>
	void update_min(CostArgs&&... cost_args)
	{
		double cost = static_cast<ConcreteCost&>(*this).get_cost(std::forward<CostArgs>(cost_args)...);
		if(cost < min_) {min_ = cost;}
	}

	void update_max()
	{
		if(min_ > max_) {max_ = min_;}
		min_ = std::numeric_limits<double>::max();
	}

	template<typename... UpdateArgs>
	void update(UpdateArgs&&...) {}

	bool check_update()
	{
		double current_min = min_;
		min_ = std::numeric_limits<double>::max();
		return (current_min > threshold_*max_);
	}

	void set_threshold(double new_threshold) {threshold_ = new_threshold;}

};*/


template<int LocalDim, int EmbedDim>
struct CostObjBase {
    using dcel_t = DCEL<LocalDim, EmbedDim>;
    using cell_t = typename dcel_t::cell_t;
    virtual ~CostObjBase() = default;

    virtual double get_cost(const std::vector<std::unordered_set<cell_t*>>& elems, const std::unordered_map<int, std::pair<bool,int>>& cavity_info = {}) const = 0;
    virtual void update(const std::unordered_set<cell_t*>& elems) = 0;
};



template<int LocalDim, int EmbedDim>
class DataEquiCost : public CostObjBase<LocalDim, EmbedDim> {
public:
    using dcel_t  = DCEL<LocalDim, EmbedDim>;
    using node_t = typename dcel_t::node_t;
    using halfedge_t = typename dcel_t::halfedge_t;
    using cell_t  = typename dcel_t::cell_t;
    using coords_t = Eigen::Matrix<double, 1, EmbedDim>;

    enum class PointTriangleRelation {
          Outside,
          Inside,
          OnAB,OnBC,OnCA,
          OnA,OnB,OnC
    };


    DataEquiCost(dcel_t& dcel, const Eigen::Matrix<double, Eigen::Dynamic, EmbedDim>& data_points): dcel_(dcel), data_points_(data_points) {
        setup_();
    }

    // compute cost on modified cells
    double get_cost(const std::vector<std::unordered_set<cell_t*>>& elems, const std::unordered_map<int, std::pair<bool,int>>& cavity_info={}) const override {
        double disp_cost = 0.0;
        auto elems_modified = elems[0]; // after
        for (auto* cell : elems_modified) {
            double Nt = compute_qoi_(cell, cavity_info);
            disp_cost += (mean_qoi_ - Nt) * (mean_qoi_ - Nt);
        }
        return elems_modified.empty() ? 0.0 : disp_cost / elems_modified.size();
    }

    void update(const std::unordered_set<cell_t*>& elems_modified) {
          double sum_qoi = mean_qoi_ * dcel_.n_cells();
          for (auto* cell : elems_modified) {
                    auto it = cell_to_qoi_it_.find(cell);  //pair<cell, it>, it->second = iterator in multimap
                    if (it != cell_to_qoi_it_.end()) {
                              sum_qoi -= it->second->first;  // subtract old qoi from sum
                              qoi_.erase(it->second);       // remove entry from multimap
                              cell_to_qoi_it_.erase(it);    
                    }
                    double new_qoi = compute_qoi_(cell);
                    sum_qoi += new_qoi;
                    // insert new qoi in multimap and update mapping
                    auto new_it = qoi_.emplace(new_qoi, cell);
                    cell_to_qoi_it_[cell] = new_it;
          }
          mean_qoi_ = sum_qoi / dcel_.n_cells();
    }

private:
    dcel_t& dcel_;
    const Eigen::Matrix<double, Eigen::Dynamic, EmbedDim> data_points_;
    std::multimap<double, cell_t*> qoi_;  // ordered for increasing qoi_ for adaptive cycle
    std::unordered_map<cell_t*, typename std::multimap<double, cell_t*>::iterator> cell_to_qoi_it_;
    double mean_qoi_;

    void setup_() {
        double sum = 0.0;
        for (auto it = dcel_.cells_begin(); it != dcel_.cells_end(); ++it) {
            cell_t* c = &(*it);
            double q = compute_qoi_(c);
            auto it_new = qoi_.emplace(q, c);
            cell_to_qoi_it_[c] = it_new;
            sum += q;
        }
        mean_qoi_ = sum / dcel_.n_cells();
    }

    double compute_qoi_(cell_t* c, const std::unordered_map<int, std::pair<bool,int>>& cavity_info = {}) const {
        
        // calculate the bounding box of c
        auto* h  = c->halfedge();
        double xmin = std::min({h->node()->coords()(0), h->next()->node()->coords()(0), h->prev()->node()->coords()(0)});
        double xmax = std::max({h->node()->coords()(0), h->next()->node()->coords()(0), h->prev()->node()->coords()(0)});
        double ymin = std::min({h->node()->coords()(1), h->next()->node()->coords()(1), h->prev()->node()->coords()(1)});
        double ymax = std::max({h->node()->coords()(1), h->next()->node()->coords()(1), h->prev()->node()->coords()(1)});

        double Nt = 0.0;
        for (int i=0; i<data_points_.rows(); ++i) {
            coords_t p = data_points_.row(i);
            if (p(0) < xmin || p(0) > xmax || p(1) < ymin || p(1) > ymax) continue;  // outside the bounding box
            auto rel= point_in_triangle_(c->halfedge()->node()->coords(), c->halfedge()->next()->node()->coords(), c->halfedge()->prev()->node()->coords(), p);
            int patch_size=0;
            if (rel != PointTriangleRelation::Outside) {
                if(rel== PointTriangleRelation::Inside) {
                    patch_size = 1;
                }
                else if (rel == PointTriangleRelation::OnAB || rel == PointTriangleRelation::OnBC || rel == PointTriangleRelation::OnCA) {
                    halfedge_t* edge = nullptr;
                    if (rel == PointTriangleRelation::OnAB) edge = c->halfedge();
                    else if (rel == PointTriangleRelation::OnBC) edge = c->halfedge()->next();
                    else edge = c->halfedge()->prev();
                    if (cavity_info.empty()) {
                              // real mode: use true DCEL info
                              patch_size = edge->on_boundary() ? 1 : 2;
                    } else {
                              // simulated mode: use cavity information
                              int id = edge->id();
                              auto it = cavity_info.find(id);
                              if(it == cavity_info.end())
                                        patch_size= 2;  // edge is internal in the cavity
                              else
                                        patch_size = it->second.first ? 1 : 2;
                    }
                }
                else{
                    node_t* v = nullptr;
                    if (rel == PointTriangleRelation::OnA) v = c->halfedge()->node();
                    else if (rel == PointTriangleRelation::OnB) v = c->halfedge()->next()->node();
                    else if (rel == PointTriangleRelation::OnC) v = c->halfedge()->prev()->node();
                    if (cavity_info.empty()) {
                              // real mode: use true DCEL info
                              patch_size = count_triangles_from_vertex_(v);
                    } else {
                              // simulated mode: use cavity information
                              int vid = v->id(); 
                              patch_size = cavity_info.find(vid)->second.second; 
                    }
                }
                Nt +=  1.0/patch_size;
            }
        }
        return Nt;
    }


    PointTriangleRelation point_in_triangle_(const coords_t& a,const coords_t& b,const coords_t& c, const coords_t& p) const {
          // check if p is on a vertex
          if(p.isApprox(a)) return PointTriangleRelation::OnA;
          if(p.isApprox(b)) return PointTriangleRelation::OnB;
          if(p.isApprox(c)) return PointTriangleRelation::OnC;
          // check if p is on an edge
          if (fdapde::internals::contains(p, a, b)) return PointTriangleRelation::OnAB;
          if (fdapde::internals::contains(p, b, c)) return PointTriangleRelation::OnBC;
          if (fdapde::internals::contains(p, c, a)) return PointTriangleRelation::OnCA;          
          // check if p is inside ABC
          Eigen::Matrix<double,3,2> tri;
          tri << a[0], a[1], b[0], b[1], c[0], c[1];
          return fdapde::internals::point_in_polygon(tri, p)? PointTriangleRelation::Inside : PointTriangleRelation::Outside;
    }

    int count_triangles_from_vertex_(const node_t* v) const {
        int count=0;
        halfedge_t* h_start = v->halfedge();
        halfedge_t* h = h_start;
        do {
            count++;
            h = h->prev()->twin();
            if(!h->cell()) // on boundary
               h = h->prev()->twin();  // skip the boundary (no cell in between halfedges)
        } while (h != h_start);
        return count;
    }
};




// ============================================================================
// SharpElemsCost
// ============================================================================
// Questo cost object misura la "qualità angolare" degli elementi attraverso il massimo coseno di un triangolo (≈ angolo più acuto o più ottuso).
// La metrica confronta direttamente lo stato *prima* e *dopo* il collasso:
//     csi  = ((1 - max_cos_before)^2) / ((1 - max_cos_after)^2)
//     cost = 1/tanh(csi) - 1/csi

// In questo modo il costo è già espresso in forma *relativa*:
// - se gli angoli peggiorano, csi cresce → il costo aumenta;
// - se gli angoli migliorano, csi si riduce → il costo si abbassa.

// Differenza con DataEquiCost:
// In DataEquiCost la quantità Nt è un valore "assoluto" che dipende dal numero di punti di dato nella cavità, quindi richiede una normalizzazione
// rispetto al massimo globale per poter essere confrontata tra step diversi. Per questo DataEquiCost eredita la logica di normalizzazione da CostObjBase.
// In SharpElemsCost invece la formula del costo contiene già un confronto esplicito tra lo stato precedente e quello simulato. Non è necessario
// riscalare i valori con un massimo globale, perché il costo è intrinsecamente "autosufficiente" come misura relativa.



template<int LocalDim, int EmbedDim>
class SharpElemsCost : public CostObjBase<LocalDim, EmbedDim> {
public:
    using dcel_t     = DCEL<LocalDim, EmbedDim>;
    using cell_t     = typename dcel_t::cell_t;
    using coords_t   = Eigen::Matrix<double,1,EmbedDim>;

    SharpElemsCost(dcel_t& dcel) : dcel_(dcel) {}

    // no normalization, since costs are already comparable between elements
    //template<typename... CostArgs>
    //double operator()(CostArgs&&... args) {
    //    return get_cost(std::forward<CostArgs>(args)...);
    //}
    // since normalization is not needed, check_update always returns false    
    //void   update_max() {}
    //bool   check_update() { return false; }

    // compute cost on modified cells versus the previous ones     
    double get_cost(const std::vector<std::unordered_set<cell_t*>>& elems, const std::unordered_map<int, std::pair<bool,int>>& cavity_info={}) const override{
          if(elems.size()<2) return 0.0;
          auto elems_before = elems[0]; // before
          auto elems_after  = elems[1]; // after
          if (elems_after.empty()) return 0.0;

          // lambda to compute a triangle cell's max cosine
          auto max_cos_in_cell = [](const cell_t* c) -> double {
              if (!c || !c->halfedge()) return 1.0; // degenerate triangle
              auto h = c->halfedge();
              const coords_t A = h ->node()->coords();
              const coords_t B = h->next()->node()->coords();
              const coords_t C = h->prev()->node()->coords();

              auto cos_at = [](const coords_t& a, const coords_t& b, const coords_t& c) {  // absolute value
                    Eigen::Matrix<double,1,EmbedDim> u = a - b;
                    Eigen::Matrix<double,1,EmbedDim> v = c - b;
                    const double nu = u.norm(), nv = v.norm();
                    if (nu==0.0 || nv==0.0) return 1.0; // degenere
                    double cc = u.dot(v)/(nu*nv);
                    return std::clamp(std::abs(cc), 0.0, 1.0);
              };

              double c1 = cos_at(C,A,B);
              double c2 = cos_at(A,B,C);
              double c3 = cos_at(B,C,A);
              return std::max({c1,c2,c3});
          };

          // before: max cavity cosine before collapse
          double max_cos_before = 0.0;
          for (auto* c : elems_before) {
             if (!c || !c->halfedge()) continue;
             double mc = max_cos_in_cell(c);
             if (mc > max_cos_before) max_cos_before = mc;
          }

          // after: max cavity cosine after collapse
          double max_cos_after = 0.0;
          for (auto* c : elems_after) {
             if (!c || !c->halfedge()) continue;
             double mc = max_cos_in_cell(c);
             if (mc > max_cos_after) max_cos_after = mc;
          }

          auto sq = [](double x){ return x*x; };
          const double num = sq(1.0 - max_cos_before); 
          const double den = sq(1.0 - max_cos_after);
          if (den <= std::numeric_limits<double>::epsilon()) 
              return std::numeric_limits<double>::infinity();   // to avoid numerical instabilities
          
          const double csi = num / den;
          if (csi <= std::numeric_limits<double>::epsilon())    // to avoid numerical instabilities
              return 0.0;

          // high cost: worse angles after collapse; low cost: improvement 
          return 1.0/std::tanh(csi) - 1.0/csi;   
    }


private:
    dcel_t& dcel_;

};






} //namespace fdapde

#endif // __ADAPTIVITY_COSTS_H__