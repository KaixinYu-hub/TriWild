// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef TRIWILD_FEATUREELEMENTS_H
#define TRIWILD_FEATUREELEMENTS_H

#include "Curves.h"

#include <array>
#include <string>
#include <vector>
#include <iostream>
#include <memory>

#include <Eigen/Dense>

#include "Point_2f.h"

namespace triwild {
    namespace feature {
        class FeatureElement 
        {
        public:
            std::vector<int> v_ids;//这条特征经过的顶点序列（用全局顶点索引存）。
            std::vector<double> paras;//对应的参数值（通常是曲线参数 t），用于把离散顶点和连续曲线对齐。
            std::string type;//曲线类型字符串（如 "Line", "RationalBezierCurve", "BezierCurve"）。
            int degree = 1;//曲线阶数（基类中初始化；具体值对不同曲线含义不同）。
            int curve_id;//一条更高层“CAD 曲线”或“语义曲线”的 id，用于分组/着色等。
            std::array<bool, 2> is_inflection = {{false, false}};//端点是否为拐点的标记（某些操作可能要断开/切分）。

            virtual void print_info() const  = 0;

            virtual Point_2f eval(double t)  const = 0; //曲线上参数 t 处的三维坐标。
            virtual Point_2f eval_first_derivative(double t)  const = 0;//一阶导
            virtual Point_2f eval_second_derivative(double t)  const = 0;//二阶导

            //反求参数，给定三维空间点 p，在区间 [t0,t1] 内找最接近 p 的参数 t（或做投影）。
            virtual double inv_eval(Point_2f& p, double t, double t0, double t1) const = 0;

            virtual std::shared_ptr<FeatureElement> simplify(const Point_2f& start, const Point_2f& end) 
            const{ return NULL; }

            virtual std::vector<double> inflection_points(const double t0, const double t1) 
            const { return std::vector<double>(); }

            //决定在 [t0,t1] 内从哪里切（曲线分段时的切分参数）。
            virtual double get_cut_t(const double t0, const double t1, bool to_flip = false) const = 0;

            virtual std::string to_maple() const { return ""; }
            virtual std::string to_eps() const;

            //衡量 [t1,t2] 子段的弯曲度/“曲性”（用来判断简化、分段策略）。
            virtual double how_curve(double t1, double t2) const = 0;
            //这段特征是否几乎退化为一个点。
            virtual bool is_point() const = 0;
            //是否“太短”（相对某阈值），可在预处理简化里直接丢弃/合并。
            virtual bool is_too_short(const Eigen::MatrixXd& V, double eps) const = 0;

            virtual ~FeatureElement(){}
            FeatureElement(const FeatureElement& other){
                v_ids = other.v_ids;
                paras = other.paras;
                type = other.type;
                degree = other.degree;
                curve_id = other.curve_id;
            }
            FeatureElement(const std::vector<int> &v_ids_, const std::vector<double> &paras_,
                    std::string type_, int degree_, int curve_id_):
                    v_ids(v_ids_), paras(paras_), type(type_), degree(degree_), curve_id(curve_id_){}

            void merge_after(const FeatureElement &other);
        };

        class Line_Feature: public FeatureElement //直线段
        {
        public:
            Line_Feature(const std::vector<int> &v_ids_, const std::vector<double> &paras_,
                    const Point_2f &start_, const Point_2f &end_, int curve_id_) :
            FeatureElement(v_ids_, paras_, "Line", 1, curve_id_),
            start(start_), end(end_) {}

            Line_Feature(const Line_Feature& other): FeatureElement(other){
                start = other.start;
                end = other.end;
            }

            void print_info() const override { std::cout<<"start="<<start<<" end="<<end<<std::endl; }

            std::string to_eps() const override;

            //eval(t)：线性插值 start + t*(end-start)。
            Point_2f eval(double t) const override;
            //一二阶导：一阶导常量（方向向量），二阶导为零。
            Point_2f eval_first_derivative(double t) const override;
            Point_2f eval_second_derivative(double t) const override;

            //把点 p 在直线段上的投影参数并裁剪到 [t0,t1]。
            double inv_eval(Point_2f& p, double t, double t0, double t1) const override;
            // double distance(Point_2f& p, double t = 0) const override;

            //线性情况通常直接 (t0+t1)/2。
            double get_cut_t(const double t0, const double t1, bool to_flip = false) const override { return (t0 + t1)/2; }

            //返回 0（直线无曲性）。
            double how_curve(double t1, double t2) const override { return 0; }

            //判断 start≈end。
            bool is_point() const override;

            //与阈值比较长度决定是否可合并。
            bool is_too_short(const Eigen::MatrixXd& V, double eps) const override;

        protected:
            Point_2f start;
            Point_2f end;
        };

        ////有理贝塞尔曲线段，CAD 圆/圆弧/椭圆等常用以 NURBS/有理 Bezier 表示。
        class RationalBezierCurve_Feature: public FeatureElement 
        {
        public:
            RationalBezierCurve_Feature(const std::vector<int> &v_ids_, const std::vector<double> &paras_,
                        const ControlVector &poles_, const ControlVector &weights_, int curve_id_) :
                    FeatureElement(v_ids_, paras_, "RationalBezierCurve", 3, curve_id_),
                    poles(poles_), weights(weights_) {}

            RationalBezierCurve_Feature(const RationalBezierCurve_Feature& other): FeatureElement(other){
                poles = other.poles;
                weights = other.weights;
            }

            Point_2f eval(double t) const override;
            Point_2f eval_first_derivative(double t)  const override;
            Point_2f eval_second_derivative(double t)  const override;
            double inv_eval(Point_2f& p, double t, double t0, double t1) const override;
            // double distance(Point_2f& p, double t = 0) const override;

            double get_cut_t(const double t0, const double t1, bool to_flip = false) const override;
            double how_curve(double t1, double t2) const override;
            bool is_point() const override;
            bool is_too_short(const Eigen::MatrixXd& V, double eps) const override;

            void print_info() const override;

        protected:
            ControlVector poles;//控制点序列（二维点）。
            ControlVector weights;//对应的权重（有理曲线的关键）。
        };

        class BezierCurve_Feature: public FeatureElement //多项式贝塞尔曲线段
        {
        public:
            BezierCurve_Feature(const std::vector<int> &v_ids_, const std::vector<double> &paras_,
                            int degree_, const ControlVector &poles_, int curve_id_) :
                    FeatureElement(v_ids_, paras_, "BezierCurve", 3, curve_id_),
                    poles(poles_) {}

            BezierCurve_Feature(const BezierCurve_Feature& other): FeatureElement(other){
                poles = other.poles;
            }

            Point_2f eval(double t) const override;
            Point_2f eval_first_derivative(double t)  const override;
            Point_2f eval_second_derivative(double t)  const override;
            double inv_eval(Point_2f& p, double t, double t0, double t1) const override;
            // double distance(Point_2f& p, double t = 0) const override;

            std::vector<double> inflection_points(const double t0, const double t1) const  override;

            double get_cut_t(const double t0, const double t1, bool to_flip = false) const override;
            double how_curve(double t1, double t2) const override;
            bool is_point() const override;
            bool is_too_short(const Eigen::MatrixXd& V, double eps) const override;

            std::string to_maple() const override;
            std::string to_eps() const override;
            void print_info() const override;

            // void compute_roots_intervals();
            std::shared_ptr<FeatureElement> simplify(const Point_2f& start, const Point_2f& end) const override;

        protected:
            ControlVector poles;
        };

    }
}


#endif //TRIWILD_FEATUREELEMENTS_H
