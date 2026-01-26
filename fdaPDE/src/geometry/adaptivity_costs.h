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


enum AdaptiveStrategy {
    NodeDensity,        
    GradientMagnitude   
};
// forward declaration of the main Metric class template
template <int LocalDim, int EmbedDim, AdaptiveStrategy Strategy>
class Metric;
template <int EmbedDim, AdaptiveStrategy S>
struct StrategyData{ };
template <int EmbedDim>
struct StrategyData<EmbedDim, AdaptiveStrategy::NodeDensity> { using storage_t = Eigen::Matrix<double, Eigen::Dynamic, EmbedDim>; };
template <int EmbedDim>
struct StrategyData<EmbedDim, AdaptiveStrategy::GradientMagnitude> {
    struct DataPoint {   
        Eigen::Matrix<double, 1, EmbedDim> coords;
        double grad_norm;
        double dx; // dz/dx 
        double dy; // dz/dy 
    };
    using storage_t = std::vector<DataPoint>;
};

Eigen::Matrix<double, Eigen::Dynamic, 2> convex_hull(const Eigen::Matrix<double, Eigen::Dynamic, 2>& boundary, bool keep_collinear);






template<int LocalDim, int EmbedDim, AdaptiveStrategy Strategy>
class NodeMetric {
public:
    using dcel_t  = DCEL<LocalDim, EmbedDim>;
    using node_t = typename dcel_t::node_t;
    using halfedge_t = typename dcel_t::halfedge_t;
    using cell_t  = typename dcel_t::cell_t;
    using coords_t = Eigen::Matrix<double, 1, EmbedDim>;

    using storage_t = typename StrategyData<EmbedDim, Strategy>::storage_t;

    // constructor
    NodeMetric(dcel_t& dcel, const storage_t& data_points, double max_area): dcel_(dcel), data_points_(data_points)  {
        if constexpr (Strategy == AdaptiveStrategy::NodeDensity) {
            kdtree_ = KDTree<EmbedDim>(data_points_);
        }
        else if constexpr (Strategy == AdaptiveStrategy::GradientMagnitude) {
            Eigen::Matrix<double, Eigen::Dynamic, 2> locations(data_points_.size(), 2);
            for(int i = 0; i < data_points_.size(); ++i)
                locations.row(i) = data_points_[i].coords;
            kdtree_ = KDTree<EmbedDim>(locations);
        }

        h_target_ = std::sqrt(4.0 * max_area / std::sqrt(3.0));
        
        build_metric_();     
        auto stats  = length_stats(0.1);
        double mean_length = stats.mean;
        rescale_c(mean_length);        // find best c and rescale node_metric_ accordingly
        build_metric_();  
  
    }
    
    struct MetricParams {
        // to avoid extreme triangles
        double rho_min = 1e-6;  // lower cap for density
        double rho_max = 1e+6;  // upper cap for density

        double Lmin = 0.7;   //   lM < Lmin  -> collapse candidate
        double Lmax = 1.3;   //   lM > Lmax  -> split candidate
        // global scaling factor for the metric: M_i = (c * rho_i) * I
        double c = 1.0;
    };

    struct LenStats {
        double mean = 1.0;
        double tol = 0.1;
        double frac_in_band = 0.0;   // percentage of edges e with |l(e)-1| <= tol
    };

    enum class PointTriangleRelation {
          Outside,
          Inside,
          OnAB,OnBC,OnCA,
          OnA,OnB,OnC
    };


    void update_metric(const std::unordered_set<cell_t*>& cells_modified, const std::unordered_set<node_t*>& deleted_nodes) {
        std::unordered_set<node_t*> nodes_modified;

        for (auto* c : cells_modified) {
            if (!c || !c->halfedge()) continue;
            nodes_modified.insert(c->halfedge()->node());
            nodes_modified.insert(c->halfedge()->next()->node());
            nodes_modified.insert(c->halfedge()->prev()->node());
        }
        update_metric(nodes_modified, deleted_nodes);
    }

    // update the metric only at the given nodes
    void update_metric(const std::unordered_set<node_t*>& nodes_modified, const std::unordered_set<node_t*>& deleted_nodes) {
        auto clamp_rho = [&](double x) {
            return std::max(metPar_.rho_min, std::min(metPar_.rho_max, x));
        };

        for (node_t* v: deleted_nodes){
            node_metric_.erase(v);
        }

        for (node_t* v : nodes_modified) {
            //double rho_mean = 0.0;
            //for (auto &kv : node_resolution_metric_all_()) {rho_mean += kv.second;}
            //rho_mean /= std::max<size_t>(1, node_resolution_metric_all_().size());
            //auto clamp = [](double x, double a, double b){ return std::max(a, std::min(b,x)); };

            double rho_v; 
            if constexpr(Strategy==AdaptiveStrategy::NodeDensity) 
                rho_v = node_density_(v);
            else if constexpr(Strategy==AdaptiveStrategy::GradientMagnitude){
                rho_v = node_resolution_metric_single_(v);
            }
            //double rho = clamp(rho_v, metPar_.rho_min, metPar_.rho_max);
            //double h_i = h_target_ / std::sqrt(rho / rho_mean);
            //h_i = std::min(h_i, h_target_);
            //double lambda = 1.0 / (h_i * h_i);
            //node_metric_[v] = lambda;
            node_metric_[v] = rho_v;
        }
    }

    // rescale the global factor c so that the mean metric edge length is ~ 1.
    void rescale_c(double mean_length, double h_target=1.0) {
        metPar_.c *= (h_target / mean_length);
    }

    double metric_edge_length(halfedge_t* e) const {
        node_t* A = e->node();
        node_t* B = e->twin()->node();
        auto ia = node_metric_.find(A);
        auto ib = node_metric_.find(B);

        const double lam = 0.5 * (ia->second + ib->second);
        coords_t AB = B->coords() - A->coords();
        return std::sqrt(std::max(0.0, lam)) * AB.norm();
    }

    LenStats length_stats(double tol) const {
        double sum = 0.0; 
        int count = 0;
        int in_band = 0;
        for (auto it = dcel_.halfedges_begin(); it != dcel_.halfedges_end(); ++it) {
            halfedge_t* e = &(*it);
            if (e->id() >= e->twin()->id()) continue;
            double len = metric_edge_length(e);
            sum += len;
            ++count;
            if (std::abs(len - 1.0) <= tol) {
                ++in_band;
            }
        }
        LenStats S;   
        if (count > 0) {  // or dcel_.n_edges()
            S.mean = sum / count;
            S.frac_in_band = static_cast<double>(in_band) / count;
        }
        return S;
    }

    std::unordered_map<node_t*, Eigen::Matrix<double,EmbedDim, EmbedDim>> node_metric() const{
        std::unordered_map<node_t*, Eigen::Matrix<double,EmbedDim, EmbedDim>> nodal_metric;
        for (auto it = dcel_.nodes_begin(); it != dcel_.nodes_end(); ++it) {
            node_t* v = &(*it);
            Eigen::Matrix<double,EmbedDim, EmbedDim> M = Eigen::Matrix<double,EmbedDim, EmbedDim>::Identity() * node_metric_.at(v);  
            nodal_metric.emplace(v, M);
        }
        return nodal_metric;
    }

    MetricParams metPar() const { return metPar_; }
    void set_Lmin(double Lmin) { metPar_.Lmin = Lmin; }
    void set_Lmax(double Lmax) { metPar_.Lmax = Lmax; }

    std::unordered_map<node_t*, Eigen::Matrix<double,EmbedDim, EmbedDim>> node_density_knn() const 
    {
        std::unordered_map<node_t*, Eigen::Matrix<double,EmbedDim, EmbedDim>> rho;
        //std::unordered_map<node_t*, double> rho;
        rho.reserve(dcel_.n_nodes());

        auto it = dcel_.nodes_begin();
        coords_t minBB = it->coords();
        coords_t maxBB = it->coords();
        ++it;
        for (; it != dcel_.nodes_end(); ++it) {
            const coords_t& p = it->coords();
            for (int d = 0; d < EmbedDim; ++d) {
                if (p(d) < minBB(d)) minBB(d) = p(d);
                if (p(d) > maxBB(d)) maxBB(d) = p(d);
            }
        }
        int N = data_points_.rows();
        Eigen::Vector2d diag = maxBB - minBB;
        double L = diag.norm();
        double h = L / std::sqrt(N); 
        double R = 2.0 * h; 
        int k = std::sqrt(N);  

        for (auto it = dcel_.nodes_begin(); it != dcel_.nodes_end(); ++it)
        {
            node_t* v = &(*it);
            coords_t x = v->coords();

            std::unordered_set<int> ids;

            // expand box until it contains at least k points
            while (true)
            {
                typename KDTree<EmbedDim>::RangeType box;
                box.ll = x.array() - R;
                box.ur = x.array() + R;

                ids = kdtree_.range_search(box);

                if ((int)ids.size() >= k)
                    break;

                R *= 2;
                if (R > 1e6) break;
            }

            // calculate exact distance to the found points
            std::vector<double> dists;
            dists.reserve(ids.size());
            for (int idx : ids)
            {
                coords_t p = data_points_.row(idx);
                dists.push_back((p - x).norm());  //euclidean distance
            }

            std::nth_element(dists.begin(), dists.begin() + (k-1), dists.end());
            double Rk = dists[k-1];   // distance from k-th nearest
            // KNN density
            double density = double(k) / (M_PI * Rk * Rk);

            rho[v] = 1/density * Eigen::Matrix<double,EmbedDim,EmbedDim>::Identity();
        }

        return rho;
    }

    // returns a map associating each mesh node with an isotropic metric tensor derived from a k-nearest-neighbor density estimate of the data locations,
    // where distances are computed geodesically along the mesh topology
    std::unordered_map<node_t*, Eigen::Matrix<double,EmbedDim, EmbedDim>> node_density_knn_geodesic() {

        std::vector<node_t*> nodes_;
        nodes_.reserve(dcel_.n_nodes());
        for (auto it = dcel_.nodes_begin(); it != dcel_.nodes_end(); ++it)
            nodes_.push_back(&(*it));
        node_coords_.resize(nodes_.size(), EmbedDim);
        for (size_t i = 0; i < nodes_.size(); ++i)
            node_coords_.row(i) = nodes_[i]->coords();

        auto kdtree_nodes_ = KDTree<EmbedDim>(node_coords_);   // KDTREE on nodes coords, not data

        auto nearest_mesh_node = [&](const coords_t& p) -> node_t* {
            auto it = kdtree_nodes_.nn_search(p.transpose());
            if (!it) return nullptr;
            int idx = *it;
            return nodes_[idx];   // node at position idx, not its actual id
        };

        // associate data to nearest node          
        std::unordered_map<node_t*, std::vector<int>> data_on_node;
        for (int i = 0; i < data_points_.rows(); ++i) {
            coords_t p = data_points_.row(i);
            node_t* v = nearest_mesh_node(p);
            if (!v) continue;
            data_on_node[v].push_back(i);
        }

        // node density via geodesic KNN
        std::unordered_map<node_t*, Eigen::Matrix<double,EmbedDim, EmbedDim>> rho;
        int N = data_points_.rows();
        int k = std::sqrt(N);   

        struct QueueItem {
            node_t* v;
            double dist;
        };
        auto comparison = [](const QueueItem& a, const QueueItem& b) {
            return a.dist > b.dist;
        };

        for (auto it = dcel_.nodes_begin(); it != dcel_.nodes_end(); ++it)
        {
            node_t* v0 = &(*it);
            std::priority_queue<QueueItem, std::vector<QueueItem>, decltype(comparison)> pq(comparison);
            std::unordered_map<node_t*, double> best;        // best geodesic distance from v0
            std::unordered_map<node_t*, node_t*> parent;     // predecessor to reconstruct shortest path

            best.reserve(dcel_.n_nodes());
            parent.reserve(dcel_.n_nodes());

            best[v0] = 0.0;
            parent[v0] = nullptr;
            pq.push({v0, 0.0});

            // KNN stopping criterion 
            int found = 0;
            double Rk = 0.0;
            node_t* v_k = nullptr; // node that "closes" found>=k (determines Rk)

            // Dijkstra loop 
            while (!pq.empty() && found < k)
            {
                auto [v, dist] = pq.top();
                pq.pop();

                // discard outdated queue entries
                auto itb = best.find(v);
                if (itb == best.end() || dist != itb->second) continue;

                // data associated to this node
                auto it_map = data_on_node.find(v);
                int n_here = (it_map != data_on_node.end()) ? int(it_map->second.size()) : 0;

                if (n_here > 0)
                {
                    found += n_here;
                    Rk = dist;
                    v_k = v; 
                }

                // expand to neighboring nodes
                for (halfedge_t* h : halfedges_from_node_(v))
                {
                    node_t* u = h->twin()->node();
                    if (!u) continue;

                    double w = (u->coords() - v->coords()).norm();
                    if (!std::isfinite(w) || w < 0.0) continue;

                    double candidate = dist + w;

                    auto itu = best.find(u);
                    if (itu == best.end() || candidate < itu->second)
                    {
                        best[u] = candidate;
                        parent[u] = v;           // store predecessor
                        pq.push({u, candidate});
                    }
                }
            }

            if (!v_k || found == 0)  // fallback metric if no data is reached
            {
                rho[v0] = Eigen::Matrix<double,EmbedDim,EmbedDim>::Identity();
                continue;
            }

            // change Rk=0 (all k points mapped onto v0) to something reasonable (Rk goes at denominator)
            if (Rk <= 0.0) 
                Rk = local_mean_edge_length_(v0)*0.3;

            // compute density + isotropic metric tensor 
            double density = double(k) / (double(N) * M_PI * Rk * Rk);
            if (!std::isfinite(density) || density <= 0.0) {
                // fallback if something went weird
                rho[v0] = Eigen::Matrix<double,EmbedDim,EmbedDim>::Identity();
                continue;
            }
            rho[v0] = (1.0 / density) * Eigen::Matrix<double,EmbedDim,EmbedDim>::Identity();
        }

        return rho;
    }

    std::list<halfedge_t*> halfedges_from_node_(node_t* n) const{
        std::list<halfedge_t*> star;
        if(n->halfedge()->cell()==nullptr) n->set_halfedge(n->halfedge()->twin()->next()); // in case the boundary node starts from an external halfedge
        halfedge_t* h_start = n->halfedge();
        halfedge_t* h = h_start;
        do {
            star.push_back(h);
            h = h->prev()->twin();
            if(!h->cell() && h->on_boundary()) {h=h->prev()->twin(); } // in case of boundary node don't count the external halfedge
        } while (h != h_start);
        return star;
    }

    double local_mean_edge_length_(node_t* v0) const{
        double sum = 0.0;
        int cnt = 0;

        for (halfedge_t* h : halfedges_from_node_(v0)) {
            node_t* u = h->twin()->node();
            if (!u) continue;

            double l = (u->coords() - v0->coords()).norm();
            if (std::isfinite(l) && l > 0.0) {
                sum += l;
                cnt++;
            }
        }

        return sum / cnt;
    }



private:
    dcel_t& dcel_;
    const storage_t data_points_;
    Eigen::Matrix<double, Eigen::Dynamic, EmbedDim> node_coords_;
    MetricParams metPar_;
    std::unordered_map<node_t*, double> node_metric_; // node* -> lambda
    KDTree<EmbedDim> kdtree_;
    double h_target_;

    std::multimap<double, cell_t*> qoi_;  // ordered for increasing qoi_ for adaptive cycle   // SERVE LA MULTIMAP ?????
    std::unordered_map<cell_t*, typename std::multimap<double, cell_t*>::iterator> cell_to_qoi_it_;


    void build_metric_() {

        std::unordered_map<node_t*, double> rho_i;
        if constexpr(Strategy==AdaptiveStrategy::NodeDensity) 
            rho_i = node_density_();
        else if constexpr(Strategy==AdaptiveStrategy::GradientMagnitude){
            rho_i = node_resolution_metric_all_();
        }
        auto clamp = [](double x, double a, double b){ return std::max(a, std::min(b,x)); };

        double rho_mean = 0.0;
        for (auto &kv : rho_i) {rho_mean += kv.second;}
        rho_mean /= std::max<size_t>(1, rho_i.size());

        node_metric_.clear();
        for (auto it = dcel_.nodes_begin(); it != dcel_.nodes_end(); ++it) {
            node_t* v = &(*it);
            double rho = clamp(rho_i.at(v), metPar_.rho_min, metPar_.rho_max);
            double h_i = h_target_ / std::sqrt(rho / rho_mean);
            h_i = std::min(h_i, h_target_);
            double lambda = 1.0 / (h_i * h_i);
            lambda = 1/lambda;
            node_metric_.emplace(v, lambda);
        }
    }



    // ----------------------------------- Metric based on density of data points -----------------------------------

    // nodal density: rho(node) = (weighted point count at node) / (barycentric dual area)
    std::unordered_map<node_t*, double> node_density_() const {
        std::unordered_map<node_t*, double> nodal_count; // accumulates weights
        std::unordered_map<node_t*, double> nodal_area;  // barycentric dual area

        // init maps with all nodes present (0.0)
        for (auto it = dcel_.nodes_begin(); it != dcel_.nodes_end(); ++it) {
            node_t* v = &(*it);
            nodal_count.emplace(v, 0.0);
            nodal_area.emplace(v, 0.0);
        }

        // precompute barycentric dual area: A*(i) = sum_T( area(T)/3 ) over incident triangles
        for (auto it = dcel_.cells_begin(); it != dcel_.cells_end(); ++it) {
            const cell_t* c = &(*it);
            node_t* A = c->halfedge()->node();
            node_t* B = c->halfedge()->next()->node();
            node_t* C = c->halfedge()->prev()->node();
            const double area = fdapde::internals::measure_2d_tri(A->coords(), B->coords(), C->coords());
            if (area <= 0.0) continue;
            nodal_area[A] += area/3.0;
            nodal_area[B] += area/3.0;
            nodal_area[C] += area/3.0;
        }

        auto add = [&](node_t* v, double w){
            auto it = nodal_count.find(v);
            if (it != nodal_count.end()) it->second += w;
            else nodal_count.emplace(v, w); // safety (in case mesh changed)
        };

        cell_t* c_start = &(*dcel_.cells_begin());
        for (int ip = 0; ip < data_points_.rows(); ++ip) {
            coords_t P = data_points_.row(ip);
            cell_t*  c = dcel_.find_cell(P, c_start );
            if (!c) continue;
            auto* h = c->halfedge();
            node_t* A = h->node();
            node_t* B = h->next()->node();
            node_t* C = h->prev()->node();

            const auto rel = point_in_triangle_(A->coords(),B->coords(),C->coords(),P);

            switch (rel) {
                case PointTriangleRelation::Inside:
                    // equally to the 3 vertices (simple, uses your helpers only)
                    add(A, 1.0/3.0); add(B, 1.0/3.0); add(C, 1.0/3.0);
                    break;

                case PointTriangleRelation::OnAB:
                    add(A, 0.5); add(B, 0.5);
                    break;
                case PointTriangleRelation::OnBC:
                    add(B, 0.5); add(C, 0.5);
                    break;
                case PointTriangleRelation::OnCA:
                    add(C, 0.5); add(A, 0.5);
                    break;

                case PointTriangleRelation::OnA:
                    add(A, 1.0);
                    break;
                case PointTriangleRelation::OnB:
                    add(B, 1.0);
                    break;
                case PointTriangleRelation::OnC:
                    add(C, 1.0);
                    break;

                default: // Outside
                    break;
            }
            c_start = c; // next search starts from here

        }

        // convert counts to densities
        for (auto& kv : nodal_count) {
            const node_t* v = kv.first;
            //std::cout << "Node " << v->id() << ": count = " << kv.second << ", area = " << nodal_area[v] << std::endl;
            const double  area = 1; //nodal_area[v];
            kv.second = (area > 0.0) ? (kv.second / area) : 0.0;
        }

        return nodal_count; 
    }

    // overloaded version to calculate the density in only one node (faster if it's not all nodes)
    double node_density_(node_t* v) const{
        if (!v) return 0.0;

        double area = 0.0;
        std::list<const cell_t*> star_around;

        double xmin =  std::numeric_limits<double>::infinity();
        double ymin =  std::numeric_limits<double>::infinity();
        double xmax = -std::numeric_limits<double>::infinity();
        double ymax = -std::numeric_limits<double>::infinity();

        if (halfedge_t* h0 = v->halfedge()) {
            halfedge_t* h = h0;
            do {
                if (cell_t* c = h->cell()) {
                    star_around.push_back(c);
                    auto* he = c->halfedge();
                    const coords_t& A = he->node()->coords();
                    const coords_t& B = he->next()->node()->coords();
                    const coords_t& C = he->prev()->node()->coords();

                    const double aT = fdapde::internals::measure_2d_tri(A, B, C);
                    if (aT > 0.0) area += aT / 3.0;

                    xmin = std::min({xmin, A(0), B(0), C(0)});
                    xmax = std::max({xmax, A(0), B(0), C(0)});
                    ymin = std::min({ymin, A(1), B(1), C(1)});
                    ymax = std::max({ymax, A(1), B(1), C(1)});
                }
                
                h = h->prev()->twin();
                if(!h->cell()){ // on boundary
                    if(h==h0 || h->twin()==h0) break; 
                    h = h->prev()->twin();  // skip the boundary (no cell in between halfedges)
                }
            } while (h != h0);
        }

        if (area <= 0.0 || star_around.empty()) 
            return 0.0;

        typename KDTree<EmbedDim>::RangeType q;
        q.ll << xmin, ymin; q.ur << xmax, ymax;
        const auto idxs = kdtree_.range_search(q);

        double count = 0.0;
        for (int i : idxs) {
            const coords_t P = data_points_.row(i);
            for (const cell_t* c : star_around) {
                auto* he = c->halfedge();
                const node_t* A = he->node();
                const node_t* B = he->next()->node();
                const node_t* C = he->prev()->node();
                const auto rel = point_in_triangle_(A->coords(), B->coords(), C->coords(), P);
                if (rel == PointTriangleRelation::Outside) continue;

                double w = 0.0;
                switch (rel) {
                    case PointTriangleRelation::Inside: w = (v==A||v==B||v==C) ? 1.0/3.0 : 0.0; break;
                    case PointTriangleRelation::OnA:    w = (v==A) ? 1.0 : 0.0; break;
                    case PointTriangleRelation::OnB:    w = (v==B) ? 1.0 : 0.0; break;
                    case PointTriangleRelation::OnC:    w = (v==C) ? 1.0 : 0.0; break;
                    case PointTriangleRelation::OnAB:   w = (v==A||v==B) ? 0.5 : 0.0; break;
                    case PointTriangleRelation::OnBC:   w = (v==B||v==C) ? 0.5 : 0.0; break;
                    case PointTriangleRelation::OnCA:   w = (v==C||v==A) ? 0.5 : 0.0; break;
                    default: break;
                }
                count += w;
                break;
            }
        }

        count = count / area;
        return count;
    }


    // cell density: rho(T) = Nt(T) / area(T)
    std::unordered_map<const cell_t*, double> cell_density_() const {
        std::unordered_map<const cell_t*, double> rho;
        rho.reserve(dcel_.n_cells());

        for (auto it = dcel_.cells_begin(); it != dcel_.cells_end(); ++it) {
            cell_t* c = &(*it);
            const double  A = fdapde::internals::measure_2d_tri(c->halfedge()->node()->coords(),c->halfedge()->next()->node()->coords(),c->halfedge()->prev()->node()->coords());
            const double  Nt = compute_qoi_(c); 
            rho.emplace(c, (A > 0.0) ? (Nt / A) : 0.0);
        }
        return rho; 
    }

    double compute_qoi_(cell_t* c, const std::unordered_map<int, std::pair<bool,int>>& cavity_info = {}) const {
        
        // calculate the bounding box of c
        auto* h  = c->halfedge();
        auto A = h->node()->coords();
        auto B = h->next()->node()->coords();
        auto C = h->prev()->node()->coords();
        double xmin = std::min({A(0), B(0), C(0)});
        double xmax = std::max({A(0), B(0), C(0)});
        double ymin = std::min({A(1), B(1), C(1)});
        double ymax = std::max({A(1), B(1), C(1)});

        double Nt = 0.0;
        typename KDTree<EmbedDim>::RangeType q;
        q.ll << xmin, ymin; q.ur << xmax, ymax;
        const auto idxs = kdtree_.range_search(q);

        //for (int i=0; i<data_points_.rows(); ++i) {
        for(int i : idxs) {
            coords_t p = data_points_.row(i);
            
            //cell_t* cell = dcel_.find_cell(p, c);
            //if(cell!=c) continue;  // point outside the cell c

            auto rel= point_in_triangle_(A, B, C, p);
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
               if(h==h_start || h->twin()==h_start) break; 
               h = h->prev()->twin();  // skip the boundary (no cell in between halfedges)
            }
        } while (h != h_start);
        return count;
    }

    // ----------------------------------- Metric based on gradient magnitude of data points -----------------------------------
    
    // node density based on gradient magnitude of data_points_
    std::unordered_map<const node_t*, double> node_resolution_metric_all_(double alpha = 2.0, double eps = 0.05) const {
        using DataPoint = typename StrategyData<EmbedDim, AdaptiveStrategy::GradientMagnitude>::DataPoint;
        std::unordered_map<const node_t*, double> nodal_rho;
        double r_initial = r_initial_();
        auto minmax_it = std::minmax_element(data_points_.begin(),data_points_.end(),[](const DataPoint& a, const DataPoint& b) {return a.grad_norm < b.grad_norm;});
        double gmin = minmax_it.first->grad_norm;
        double gmax = minmax_it.second->grad_norm;
        const double RHO_MAX_CALC = rho_from_grad_(gmax, gmin, gmax, alpha, eps); // norm=1 -> max rho
        const double RHO_MIN_CALC = rho_from_grad_(gmin, gmin, gmax, alpha, eps); // norm=0 -> min rho
        double rho_range_calc = RHO_MAX_CALC - RHO_MIN_CALC;
        if (rho_range_calc < 1e-12) { // if all gradients are equal, assign uniform density
            for (auto it = dcel_.nodes_begin(); it != dcel_.nodes_end(); ++it)
                nodal_rho.emplace(&(*it), 1.0);
            return nodal_rho;
        }
        constexpr double RHO_MIN_TARGET = 0.0; 
        constexpr double RHO_MAX_TARGET = 1.0; 
        double rho_range_target = RHO_MAX_TARGET - RHO_MIN_TARGET;
        for (auto it = dcel_.nodes_begin(); it != dcel_.nodes_end(); ++it) {
            const node_t* v = &(*it);
            const coords_t& P = v->coords();
            const double g_interpolated = shepard_interpolation_kd_combined_(P, r_initial);
            const double rho_calc = rho_from_grad_(g_interpolated, gmin, gmax, alpha, eps);
            double rho_normalized = ((rho_calc - RHO_MIN_CALC) / rho_range_calc) * rho_range_target + RHO_MIN_TARGET;
            nodal_rho.emplace(v, rho_normalized);
        }
        return nodal_rho;
        //for (auto it = dcel_.nodes_begin(); it != dcel_.nodes_end(); ++it) {
        //    const node_t* v = &(*it);
        //    const coords_t& P = v->coords();
        //    const double g_interpolated = shepard_interpolation_kd_combined_(P, r_initial);
        //    const double rho = rho_from_grad_(g_interpolated, gmin, gmax, alpha, eps);
        //    nodal_rho.emplace(v, rho);
        //}
        //return nodal_rho;
    }

    double node_resolution_metric_single_(const node_t* v, double alpha = 2.0, double eps = 0.05) const {
        using DataPoint = typename StrategyData<EmbedDim, AdaptiveStrategy::GradientMagnitude>::DataPoint;
        if (!v) return 0.0;
        const double g_interpolated = shepard_interpolation_kd_combined_(v->coords(), r_initial_());
        auto minmax_it = std::minmax_element(data_points_.begin(),data_points_.end(),[](const DataPoint& a, const DataPoint& b) {return a.grad_norm < b.grad_norm;});
        double gmin = minmax_it.first->grad_norm;
        double gmax = minmax_it.second->grad_norm;
        const double RHO_MAX_CALC = rho_from_grad_(gmax, gmin, gmax, alpha, eps); // norm=1 -> max rho
        const double RHO_MIN_CALC = rho_from_grad_(gmin, gmin, gmax, alpha, eps); // norm=0 -> min rho
        double rho_range_calc = RHO_MAX_CALC - RHO_MIN_CALC;
        constexpr double RHO_MIN_TARGET = 0.0; 
        constexpr double RHO_MAX_TARGET = 1.0; 
        double rho_range_target = RHO_MAX_TARGET - RHO_MIN_TARGET;
        const double rho_calc = rho_from_grad_(g_interpolated, gmin, gmax, alpha, eps);
        double rho_normalized = ((rho_calc - RHO_MIN_CALC) / rho_range_calc) * rho_range_target + RHO_MIN_TARGET;       
        return rho_normalized;
    }

    double rho_from_grad_(double g, double gmin, double gmax, double alpha=6.0, double eps=0.1) const{
        double norm = (g - gmin) / (gmax - gmin + 1e-12);
        norm = std::clamp(norm, 0.0, 1.0);
        return std::pow(norm + eps, alpha);
    }

    // function that implements Shepard interpolation with slope correction (A two-dimensional interpolation function for irregularly-spaced data, D. Shepard)
    // given the coordinates of point P, based on distribution of data points, returns the interpolated value f(P)
    double shepard_interpolation_kd_combined_(const coords_t& P, double r_initial) const{

        using DataPoint = typename StrategyData<EmbedDim, AdaptiveStrategy::GradientMagnitude>::DataPoint;
        const size_t C_MIN = 4;
        const size_t C_MAX = 10;

        double r_sq = r_initial * r_initial;
        typename KDTree<EmbedDim>::RangeType q;
        q.ll[0] = P(0) - r_initial;
        q.ll[1] = P(1) - r_initial;
        q.ur[0] = P(0) + r_initial;
        q.ur[1] = P(1) + r_initial;
        const auto idxs = kdtree_.range_search(q);

        std::vector<const DataPoint*> C_P;      
        C_P.reserve(idxs.size());
        for (int i : idxs) {
            const DataPoint& D_i = data_points_[i];
            if (distance_sq_(P, D_i.coords) < r_sq) {
                C_P.push_back(&D_i);
            }
        }
        size_t n_CP = C_P.size();
        
        std::vector<const DataPoint*> C_prime;
        double r_prime_sq; 

        if (n_CP <= C_MIN) {  // 4 nearest points 
            NearestResults results = find_k_neighbors_(P, C_MIN);
            C_prime = results.k_neighbors;
            r_prime_sq = results.r_prime_sq;
        } else if (n_CP <= C_MAX) {
            C_prime = C_P;
            r_prime_sq = r_sq; 
        } else { // 10 nearest points
            NearestResults results = find_k_neighbors_(P, C_MAX);
            C_prime = results.k_neighbors;
            r_prime_sq = results.r_prime_sq;
        }
        
        double sum_weighted_value = 0.0;
        double sum_weights = 0.0;
        const double EPSILON_SQ = 1e-12; 

        for (const auto* D_i_ptr : C_prime) {
            const auto& D_i = *D_i_ptr;
            double d_sq = distance_sq_(P, D_i.coords);
            
            if (d_sq < EPSILON_SQ) { // too close to be seen different from D_i
                return D_i.grad_norm;
            }

            double d_i = std::sqrt(d_sq);
            double r_prime = std::sqrt(r_prime_sq);
            double W_i = 0.0;
            if (d_sq < r_prime_sq && r_prime_sq > EPSILON_SQ) {
                double num = std::sqrt(r_prime_sq) - std::sqrt(d_sq);
                W_i = num*num / (r_prime_sq * d_sq);
            }
            
            // slope term
            double dx = P(0) - D_i.coords(0);
            double dy = P(1) - D_i.coords(1);
            double delta_z_i = D_i.dx * dx + D_i.dy * dy;

            sum_weighted_value += W_i * (D_i.grad_norm + delta_z_i);
            sum_weights += W_i;
        }

        if (sum_weights > 0.0) {
            return sum_weighted_value / sum_weights;
        } else {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }

    double r_initial_() const {
        if (data_points_.size()==0)  return 0.0;
        Eigen::Matrix<double, Eigen::Dynamic, 2> locations(data_points_.size(), 2);
            for(int i = 0; i < data_points_.size(); ++i)
                locations.row(i) = data_points_[i].coords;
        auto c_hull = fdapde::convex_hull(locations, false);
        double area = std::abs(fdapde::internals::signed_measure_2d_polygon(c_hull));
        size_t N = data_points_.size();
        // 7 is the target number of points contained in a circle of radius r_initial
        return std::sqrt((7.0 * area) / (M_PI * N));
    }

    double distance_sq_(const coords_t& p1, const coords_t& p2) const {
        double dx = p1(0) - p2(0);
        double dy = p1(1) - p2(1);
        return dx * dx + dy * dy;
    }

    struct NearestResults {
        using DataPoint = typename StrategyData<EmbedDim, AdaptiveStrategy::GradientMagnitude>::DataPoint;
        std::vector<const DataPoint*> k_neighbors;
        double r_prime_sq; // Distanza al quadrato del (K+1)-esimo vicino, usata come r_prime^2
    };

    NearestResults find_k_neighbors_(const coords_t& P, int K_target) const {
        using DataPoint = typename StrategyData<EmbedDim, AdaptiveStrategy::GradientMagnitude>::DataPoint;
        auto compare_dist = [](const std::pair<double, const DataPoint*>& a, const std::pair<double, const DataPoint*>& b) {
            return a.first < b.first; 
        };
        std::priority_queue<std::pair<double, const DataPoint*>, std::vector<std::pair<double, const DataPoint*>>, decltype(compare_dist)> max_heap(compare_dist);

        double d_sq;
        for (const auto& D_i : data_points_) {
            d_sq = distance_sq_(P, D_i.coords);
            
            if (max_heap.size() < K_target) {
                max_heap.push({d_sq, &D_i});
            } else if (d_sq < max_heap.top().first) {
                max_heap.pop(); 
                max_heap.push({d_sq, &D_i});
            }
        }

        NearestResults results;
        results.k_neighbors.reserve(max_heap.size());
        // r_prime is the distance from the farthest neighbor in the heap
        results.r_prime_sq = max_heap.empty() ? 0.0 : max_heap.top().first;

        while (!max_heap.empty()) {
            results.k_neighbors.push_back(max_heap.top().second);
            max_heap.pop();
        }
        return results;
    }
    

};

























































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

    DataEquiCost(dcel_t& dcel, const Eigen::Matrix<double, Eigen::Dynamic, EmbedDim>& data_points): dcel_(dcel), data_points_(data_points), kdtree_(data_points_) {
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
    KDTree<EmbedDim> kdtree_;

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
        auto A = h->node()->coords();
        auto B = h->next()->node()->coords();
        auto C = h->prev()->node()->coords();
        double xmin = std::min({A(0), B(0), C(0)});
        double xmax = std::max({A(0), B(0), C(0)});
        double ymin = std::min({A(1), B(1), C(1)});
        double ymax = std::max({A(1), B(1), C(1)});

        double Nt = 0.0;
        typename KDTree<EmbedDim>::RangeType q;
        q.ll << xmin, ymin; q.ur << xmax, ymax;
        const auto idxs = kdtree_.range_search(q);

        //for (int i=0; i<data_points_.rows(); ++i) {
        for(int i : idxs) {
            coords_t p = data_points_.row(i);
            
            //cell_t* cell = dcel_.find_cell(p, c);
            //if(cell!=c) continue;  // point outside the cell c

            auto rel= point_in_triangle_(A, B, C, p);
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