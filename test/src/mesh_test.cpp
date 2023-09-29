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

#include <fdaPDE/utils.h>
#include <fdaPDE/mesh.h>
#include <gtest/gtest.h>   // testing framework

#include <memory>
#include <random>
#include <unordered_set>
#include <vector>
using fdapde::core::Element;
using fdapde::core::Mesh;

#include "utils/mesh_loader.h"
#include "utils/utils.h"
using fdapde::testing::MESH_TYPE_LIST;
using fdapde::testing::MeshLoader;
using fdapde::testing::almost_equal;

// test suite for testing both non-manifold meshes (2D/3D) and manifold mesh (2.5D/1.5D)
template <typename E> struct mesh_test : public ::testing::Test {
    MeshLoader<E> meshLoader {};   // use default mesh
    static constexpr unsigned int M = MeshLoader<E>::M;
    static constexpr unsigned int N = MeshLoader<E>::N;
};
TYPED_TEST_SUITE(mesh_test, MESH_TYPE_LIST);

// check points' coordinate embedded in an element are loaded correctly
TYPED_TEST(mesh_test, element_construction) {
    for (std::size_t i = 0; i < this->meshLoader.mesh.n_elements(); ++i) {
        // request element with ID i
        auto e = this->meshLoader.mesh.element(i);
        // check coordinates stored in element built from Mesh object match raw informations
        int j = 0;
	auto raw_elements = this->meshLoader.elements_.row(i);
	auto raw_points = this->meshLoader.points_;
        for (int k = 0; k < raw_elements.size(); ++k) {
            int nodeID = raw_elements[k];
	    auto p = raw_points.row(nodeID);
            SVector<TestFixture::N> ePoint = e.coords()[j];
            for (std::size_t idx = 0; idx < TestFixture::N; ++idx) {
                EXPECT_TRUE(almost_equal(p[idx], ePoint[idx]));
            }
            j++;
        }
    }
}

// performs some checks on the mesh topology, e.g. checks that stated neighbors shares exactly M points
TYPED_TEST(mesh_test, boundary_checks) {
    // cycle over all mesh elements
    for (std::size_t i = 0; i < this->meshLoader.mesh.n_elements(); ++i) {
        auto e = this->meshLoader.mesh.element(i);
        // check that neighboing elements have always M points in common
        for (int neighID : e.neighbors()) {
            if (!e.is_on_boundary()) {
                // request neighboring element from mesh
                auto n = this->meshLoader.mesh.element(neighID);
                // take nodes of both elements
                std::array<SVector<TestFixture::N>, TestFixture::M + 1> eList, nList;
                for (std::size_t j = 0; j < TestFixture::M + 1; ++j) {
                    eList[j] = e.coords()[j];
                    nList[j] = n.coords()[j];
                }
                // check that the points in common between the two are exactly M
                std::size_t matches = 0;
                for (SVector<TestFixture::N> p : eList) {
                    if (std::find(nList.begin(), nList.end(), p) != nList.end()) matches++;
                }
                EXPECT_TRUE(matches == TestFixture::M);
            } else {
                // check that at least one vertex of e is detected as boundary point
                bool element_on_boundary = false;
                auto node_ids = e.node_ids();
                for (std::size_t n : node_ids) {
                    if (this->meshLoader.mesh.is_on_boundary(n)) {   // mesh detects this point as boundary point
                        element_on_boundary = true;
                    }
                }
                EXPECT_TRUE(element_on_boundary);
            }
        }
    }
}

// check the range for loop scans the whole mesh element by element
TYPED_TEST(mesh_test, range_for) {
    // prepare set with all indexes of IDs to touch
    std::unordered_set<int> meshIDs {};
    for (int i = 0; i < this->meshLoader.mesh.n_elements(); ++i) meshIDs.insert(i);

    // range-for over all elements removing the element's ID from the above set when the element is visited
    for (const auto& e : this->meshLoader.mesh) {
        // check element ID still present in the IDs set (ID not visisted by means of a different element)
        EXPECT_TRUE(meshIDs.find(e.ID()) != meshIDs.end());
        meshIDs.erase(e.ID());
    }
    // check that no ID is left in the initial set
    EXPECT_TRUE(meshIDs.empty());
}
