// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef TRIWILD_RATIONAL_H
#define TRIWILD_RATIONAL_H

#include <gmp.h>
#include <iostream>

namespace triwild {

    class Rational {//https://cs.nyu.edu/acsys/cvc3/releases/1.5/doc/rational-gmp_8cpp-source.html
    public:
        //GMP 提供的有理数类型：分子分母都用任意精度整数存储。GMP 提供的有理数类型：分子分母都用任意精度整数存储。
        mpq_t value;
        //约分，把分子分母化成最简形式，比如 2/4 → 1/2。
        void canonicalize() {
            mpq_canonicalize(value);
        }
        //返回数的符号：正数、负数或零。
        int get_sign() {
            return mpq_sgn(value);
        }
        //默认构造，初始化为 0。
        Rational() {
            mpq_init(value);
            mpq_set_d(value, 0);
        }

        //可以从 double、已有的 GMP 有理数 mpq_t 或另一个 Rational 构造
        Rational(double d) {
            mpq_init(value);
            mpq_set_d(value, d);
//            canonicalize();
        }
        Rational(const mpq_t &v_) {
            mpq_init(value);
            mpq_set(value, v_);
//            canonicalize();
        }
        Rational(const Rational &other) {
            mpq_init(value);
            mpq_set(value, other.value);
        }
        //析构时释放内部 GMP 资源。
        ~Rational() {
            mpq_clear(value);
        }

//        //+, - another point
//        Rational operator+(const Rational &r) const {
//            Rational r_out;
//            mpq_add(r_out.value, value, r.value);
//            return r_out;
//        }
//
//        Rational operator-(const Rational &r) const {
//            Rational r_out;
//            mpq_sub(r_out.value, value, r.value);
//            return r_out;
//        }
//
//        //*, / double/rational
//        Rational operator*(const Rational &r) const {
//            Rational r_out;
//            mpq_mul(r_out.value, value, r.value);
//            return r_out;
//        }
//
//        Rational operator/(const Rational &r) const {
//            Rational r_out;
//            mpq_div(r_out.value, value, r.value);
//            return r_out;
//        }
//
//        //=
//        void operator=(const Rational &r) {
//            mpq_set(value, r.value);
//        }

        friend Rational operator+(const Rational& x, const Rational& y) {
            Rational r_out;
            mpq_add(r_out.value, x.value, y.value);
            return r_out;
        }

        friend Rational operator-(const Rational& x, const Rational& y) {
            Rational r_out;
            mpq_sub(r_out.value, x.value, y.value);
            return r_out;
        }

        friend Rational operator*(const Rational& x, const Rational& y) {
            Rational r_out;
            mpq_mul(r_out.value, x.value, y.value);
            return r_out;
        }

        friend Rational operator/(const Rational& x, const Rational& y) {
            Rational r_out;
            mpq_div(r_out.value, x.value, y.value);
            return r_out;
        }

        Rational& operator=(const Rational& x) {
            if (this == &x) return *this;
            mpq_set(value, x.value);
            return *this;
        }

        Rational& operator=(const double x) {
            mpq_set_d(value, x);
//            canonicalize();
            return *this;
        }

        //> < ==
        friend bool operator<(const Rational &r, const Rational &r1) {
            return mpq_cmp(r.value, r1.value) < 0;
        }

        friend bool operator>(const Rational &r, const Rational &r1) {
            return mpq_cmp(r.value, r1.value) > 0;
        }

        friend bool operator<=(const Rational &r, const Rational &r1) {
            return mpq_cmp(r.value, r1.value) <= 0;
        }

        friend bool operator>=(const Rational &r, const Rational &r1) {
            return mpq_cmp(r.value, r1.value) >= 0;
        }

        friend bool operator==(const Rational &r, const Rational &r1) {
            return mpq_equal(r.value, r1.value);
        }

        friend bool operator!=(const Rational &r, const Rational &r1) {
            return !mpq_equal(r.value, r1.value);
        }

        //to double 把高精度有理数转换成 double。
        double to_double() {
            return mpq_get_d(value);
        }

        //<<输出流重载，打印时直接转成 double 输出。
        friend std::ostream &operator<<(std::ostream &os, const Rational &r) {
            os << mpq_get_d(r.value);
            return os;
        }
    };
}

#endif //TRIWILD_RATIONAL_H
