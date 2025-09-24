// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef TRIWILD_FEATURE_H
#define TRIWILD_FEATURE_H

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_AABB.h>
#include <igl/writeOBJ.h>

#include "FeatureElements.h"
#include "TrimeshElements.h"

#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace triwild {
    namespace feature {
        /*主特征 (primary feature curves)
        这些曲线是“重要”的。
        要在网格中以 高阶曲边三角形 (curved triangle) 来精确表示。
        网格的目标边长由用户指定，例如 --target-edge-length-r 参数。
        这些特征必须尽可能保留，因为它们对几何和 PDE 边界条件可能至关重要。*/
        extern std::vector<std::shared_ptr<FeatureElement>> features;

        /*次特征 (secondary feature curves)
        这些曲线“不那么重要”或者太细小，没必要完全保留。
        可以用 线性分片近似 (piecewise linear mesh) 来表示。
        允许在设定的公差 ε 内近似，从而避免在这些细节处生成大量小而质量差的三角形。*/
        extern std::vector<std::shared_ptr<FeatureElement>> secondary_features;
        
        extern double feature_eps;
        extern double feature_eps_2;

        bool init(const std::string& feature_file);
        bool init(json& feature_info);
        void map_feature2mesh(MeshData& mesh);
        bool is_on_segment(const Point_2& p, const Point_2& p1, const Point_2& p2);

        void snap_vertices(MeshData& mesh);

        void merge_inflection(MeshData& mesh);
        void curving(MeshData& mesh, GEO::MeshFacetsAABB &b_tree);
        void add_nodes(MeshData& mesh);
        void subdivide_into_2(MeshData& mesh);
        void subdivide_into_3(MeshData& mesh);
        void get_new_nodes(const std::array<int, 3>& tag_feature_e, const std::array<TriVertex, 3>& vs,
                std::vector<Point_2f>& new_nodes);
        bool is_valid_inversion(const std::array<Point_2f, 3>& ps, const std::vector<Point_2f>& ns);
        void fix_inversion(MeshData& mesh);

        void check_inversion(MeshData& mesh, bool is_output_objs = false);

        void visualize_features(MeshData& mesh);
        Point_2f json2point(const json& a);
        std::array<double, 2> json2d2stdarray(const json& a);
        std::vector<double> json2d2stdvector(const json& a);
        std::vector<double> json1d2stdvector(const json& a);

        ControlVector json2d2ctrlvector(const json& a);
        ControlVector json1d2ctrlvector(const json& a);

        void output_input_features(const std::vector<std::shared_ptr<FeatureElement>>& features,
                Eigen::MatrixXd& V, std::vector<std::array<int, 2>>& edges, const std::string& postfix);
        void output_stats(MeshData& mesh, std::ofstream& f);

        inline void reset(){
            features.clear();
            secondary_features.clear();
        }
    }
}


#endif //TRIWILD_FEATURE_H
