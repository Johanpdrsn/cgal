#ifndef VISIBILITY_2_EXAMPLES_NUMERICALINTEGRALOpt_H
#define VISIBILITY_2_EXAMPLES_NUMERICALINTEGRALOpt_H

#include <cmath>
#include <iostream>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/gauss.hpp>

template<class T>
class Integral {
public:

    static T Compute(const T &T0X, const T &T0Y, const T &T1X, const T &T1Y, const T &B0X, const T &B0Y, const T &B1X,
                     const T &B1Y,
                     const T &PLX, const T &PLY, const T &PRX, const T &PRY, const T &V0X, const T V0Y, const T &V1X,
                     const T &V1Y, const T &rho0, const T &rho1) {

        auto f = [&](T rho) -> T {
            return CGAL::abs(
                    Integral<T>::Int(T0X, T0Y, T1X, T1Y, B0X, B0Y, B1X, B1Y, PRX,
                                     PLY,
                                     PRX, PRY, V0X, V0Y, V1X, V1Y, rho));
        };
        return trapezoidal(f, rho0, rho1);
    }

    static T
    ComputeOpt(const T &T0X, const T &T0Y, const T &T1X, const T &T1Y, const T &B0X, const T &B0Y, const T &B1X,
               const T &B1Y,
               const T &PLX, const T &PLY, const T &PRX, const T &PRY, const T &V0X, const T V0Y, const T &V1X,
               const T &V1Y, const T &rho0, const T &rho1) {

        auto f = [&](T rho) -> T {

            T t1 = V0Y * V0X;
            T t2 = V1Y * t1;
            T t3 = T0Y - T1Y;
            T t5 = 1 - rho;
            T t6 = t5 * V0X;
            T t7 = V1X * rho;
            T t8 = t6 + t7;
            T t11 = t5 * V0Y;
            T t12 = V1Y * rho;
            T t13 = -t11 - t12;
            T t15 = t8 * t3 - t13 * (-T0X + T1X);
            T t16 = 0.1e1 / t15;
            T t17 = t16 * t3 * rho;
            T t19 = V1Y * V1Y;
            T t20 = t19 * V0X;
            T t22 = V0Y * V0Y;
            T t23 = V1X * t22;
            T t25 = V0Y * V1X;
            T t33 = -t16 * t3 * V1Y * t1 + t17 * V1Y * t25 + t16 * t3 * t23 + t17 * t2 - t17 * t20 - t17 * t23;
            T t34 = B0Y - B1Y;
            T t35 = t34 * t34;
            T t40 = t8 * t34 - t13 * (-B0X + B1X);
            T t41 = t40 * t40;
            T t42 = 0.1e1 / t41;
            T t44 = t12 * t1;
            T t45 = rho * t20;
            T t46 = rho * t23;
            T t47 = t12 * t25;
            T t48 = t44 - t45 - t46 + t47 - t2 + t23;
            T t49 = t3 * t3;
            T t50 = t49 * t48;
            T t51 = t15 * t15;
            T t52 = 0.1e1 / t51;
            T t53 = t34 * t52;
            T t54 = 0.1e1 / t40;
            T t62 = (t11 + t12 + PRY) * PRX - (t6 + t7 + PRX) * PRY;
            T t63 = t62 * t62;
            T t64 = t63 * t63;
            T t65 = t13 * t13;
            T t66 = t65 * t65;
            T t67 = 0.1e1 / t66;
            T t69 = PLY * PLY;
            T t70 = t69 * t69;
            T t71 = t8 * t8;
            T t72 = t71 * t71;
            T t80 = T0X * T1Y - T0Y * T1X;
            T t82 = -0.1e1 / t13;
            T t84 = t82 * t16 * t13 * t80;
            T t94 = -t82 * t16 * t13 * t80 * t23 + t84 * t2 - t84 * t44 + t84 * t45 + t84 * t46 - t84 * t47;
            T t99 = B0X * B1Y - B0Y * B1X;
            T t102 = -t13 * t13 * t34;
            T t107 = 0.1e1 / t65;
            T t109 = t80 * t48;
            T t116 = t99 * t52;
            T t123 = t13 * t65;
            T t124 = -0.1e1 / t123;
            T t136 = t99 * t99;
            T t142 = t80 * t80;
            T t143 = t142 * t48;
            T t168 = t82 * t8 * PLY + t82 * t62;
            T t179 = (-t67 * t72 * t70 + t67 * t64) * (-t42 * t35 * t33 + t54 * t53 * t50) / 8 +
                     (t124 * t8 * t71 * PLY * t69 + t124 * t62 * t63) *
                     (-t107 * t42 * (-2 * t102 * t99 * t33 + t65 * t35 * t94) -
                      t54 * (2 * t53 * t82 * t3 * t13 * t109 + t13 * t116 * t82 * t50)) / 6 +
                     (-t107 * t71 * t69 + t107 * t63) *
                     (-t107 * t42 * (-2 * t102 * t99 * t94 + t65 * t136 * t33) - t54 *
                                                                                 (-2 *
                                                                                  t116 *
                                                                                  t107 *
                                                                                  t3 *
                                                                                  t65 *
                                                                                  t109 -
                                                                                  t34 *
                                                                                  t107 *
                                                                                  t52 *
                                                                                  t65 *
                                                                                  t143)) /
                     4 - t168 * t107 * t42 * t65 * t136 * t94 / 2 - t168 * t54 * t99 * t124 * t52 * t123 * t143 / 2;


            return CGAL::abs(t179);
        };
        return boost::math::quadrature::gauss<T, 100>::integrate(f, rho0, rho1);
//        return boost::math::quadrature::trapezoidal(f, rho0,rho1);;

    }

    static T IntOpt(const T &T0X, const T &T0Y, const T &T1X, const T &T1Y, const T &B0X, const T &B0Y, const T &B1X,
                    const T &B1Y,
                    const T &PLX, const T &PLY, const T &PRX, const T &PRY, const T &V0X, const T V0Y, const T &V1X,
                    const T &V1Y, const T &rho) {


        T t1 = V0Y * V0X;
        T t2 = V1Y * t1;
        T t3 = T0Y - T1Y;
        T t5 = 1 - rho;
        T t6 = t5 * V0X;
        T t7 = V1X * rho;
        T t8 = t6 + t7;
        T t11 = t5 * V0Y;
        T t12 = V1Y * rho;
        T t13 = -t11 - t12;
        T t15 = t8 * t3 - t13 * (-T0X + T1X);
        T t16 = 0.1e1 / t15;
        T t17 = t16 * t3 * rho;
        T t19 = V1Y * V1Y;
        T t20 = t19 * V0X;
        T t22 = V0Y * V0Y;
        T t23 = V1X * t22;
        T t25 = V0Y * V1X;
        T t33 = -t16 * t3 * V1Y * t1 + t17 * V1Y * t25 + t16 * t3 * t23 + t17 * t2 - t17 * t20 - t17 * t23;
        T t34 = B0Y - B1Y;
        T t35 = t34 * t34;
        T t40 = t8 * t34 - t13 * (-B0X + B1X);
        T t41 = t40 * t40;
        T t42 = 0.1e1 / t41;
        T t44 = t12 * t1;
        T t45 = rho * t20;
        T t46 = rho * t23;
        T t47 = t12 * t25;
        T t48 = t44 - t45 - t46 + t47 - t2 + t23;
        T t49 = t3 * t3;
        T t50 = t49 * t48;
        T t51 = t15 * t15;
        T t52 = 0.1e1 / t51;
        T t53 = t34 * t52;
        T t54 = 0.1e1 / t40;
        T t62 = (t11 + t12 + PRY) * PRX - (t6 + t7 + PRX) * PRY;
        T t63 = t62 * t62;
        T t64 = t63 * t63;
        T t65 = t13 * t13;
        T t66 = t65 * t65;
        T t67 = 0.1e1 / t66;
        T t69 = PLY * PLY;
        T t70 = t69 * t69;
        T t71 = t8 * t8;
        T t72 = t71 * t71;
        T t80 = T0X * T1Y - T0Y * T1X;
        T t82 = -0.1e1 / t13;
        T t84 = t82 * t16 * t13 * t80;
        T t94 = -t82 * t16 * t13 * t80 * t23 + t84 * t2 - t84 * t44 + t84 * t45 + t84 * t46 - t84 * t47;
        T t99 = B0X * B1Y - B0Y * B1X;
        T t102 = -t13 * t13 * t34;
        T t107 = 0.1e1 / t65;
        T t109 = t80 * t48;
        T t116 = t99 * t52;
        T t123 = t13 * t65;
        T t124 = -0.1e1 / t123;
        T t136 = t99 * t99;
        T t142 = t80 * t80;
        T t143 = t142 * t48;
        T t168 = t82 * t8 * PLY + t82 * t62;
        T t179 = (-t67 * t72 * t70 + t67 * t64) * (-t42 * t35 * t33 + t54 * t53 * t50) / 8 +
                 (t124 * t8 * t71 * PLY * t69 + t124 * t62 * t63) *
                 (-t107 * t42 * (-*t102 * t99 * t33 + t65 * t35 * t94) -
                  t54 * (2 * t53 * t82 * t3 * t13 * t109 + t13 * t116 * t82 * t50)) / 6 +
                 (-t107 * t71 * t69 + t107 * t63) * (-t107 * t42 * (-2 * t102 * t99 * t94 + t65 * t136 * t33) - t54 *
                                                                                                                (-2 *
                                                                                                                 t116 *
                                                                                                                 t107 *
                                                                                                                 t3 *
                                                                                                                 t65 *
                                                                                                                 t109 -
                                                                                                                 t34 *
                                                                                                                 t107 *
                                                                                                                 t52 *
                                                                                                                 t65 *
                                                                                                                 t143)) /
                 4 - t168 * t107 * t42 * t65 * t136 * t94 / 2 - t168 * t54 * t99 * t124 * t52 * t123 * t143 / 2;


        return CGAL::abs(t179);
    };

    static T Int(const T &T0X, const T &T0Y, const T &T1X, const T &T1Y, const T &B0X, const T &B0Y, const T &B1X,
                 const T &B1Y,
                 const T &PLX, const T &PLY, const T &PRX, const T &PRY, const T &V0X, const T V0Y, const T &V1X,
                 const T &V1Y, const T &rho) {
        return (-(V0X * V0Y * V1Y * rho * (T0Y - T1Y) /
                  ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) - (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho)) -
                  V0X * V1Y * V1Y * rho * (T0Y - T1Y) /
                  ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) - (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho)) -
                  V0Y * V0Y * V1X * rho * (T0Y - T1Y) /
                  ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) - (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho)) +
                  V0Y * V1X * V1Y * rho * (T0Y - T1Y) /
                  ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) - (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho)) -
                  V0X * V0Y * V1Y * (T0Y - T1Y) /
                  ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) - (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho)) +
                  V0Y * V0Y * V1X * (T0Y - T1Y) /
                  ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) - (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho))) *
                pow(B0Y - B1Y, 2) *
                pow((B0Y - B1Y) * (V0X * (1 - rho) + V1X * rho) - (-B0X + B1X) * (-V0Y * (1 - rho) - V1Y * rho), -2) /
                2 + (V0X * V0Y * V1Y * rho - V0X * V1Y * V1Y * rho - V0Y * V0Y * V1X * rho + V0Y * V1X * V1Y * rho -
                     V0X * V0Y * V1Y + V0Y * V0Y * V1X) * pow(T0Y - T1Y, 2) *
                    pow((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) - (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho),
                        -2) * (B0Y - B1Y) /
                    ((B0Y - B1Y) * (V0X * (1 - rho) + V1X * rho) - (-B0X + B1X) * (-V0Y * (1 - rho) - V1Y * rho)) / 2) *
               (pow(PRX * (V0Y * (1 - rho) + V1Y * rho + PRY) - PRY * (V0X * (1 - rho) + V1X * rho + PRX), 4) *
                pow(V0Y * (1 - rho) + V1Y * rho, -4) -
                pow(PLY, 4) * pow(V0X * (1 - rho) + V1X * rho, 4) * pow(V0Y * (1 - rho) + V1Y * rho, -4)) / 4 +
               (-((-V0X * V0Y * V1Y * rho * (T0X * T1Y - T0Y * T1X) * (-V0Y * (1 - rho) - V1Y * rho) /
                   ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) - (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho)) /
                   (V0Y * (1 - rho) + V1Y * rho) +
                   V0X * V1Y * V1Y * rho * (T0X * T1Y - T0Y * T1X) * (-V0Y * (1 - rho) - V1Y * rho) /
                   ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) - (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho)) /
                   (V0Y * (1 - rho) + V1Y * rho) +
                   V0Y * V0Y * V1X * rho * (T0X * T1Y - T0Y * T1X) * (-V0Y * (1 - rho) - V1Y * rho) /
                   ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) - (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho)) /
                   (V0Y * (1 - rho) + V1Y * rho) -
                   V0Y * V1X * V1Y * rho * (T0X * T1Y - T0Y * T1X) * (-V0Y * (1 - rho) - V1Y * rho) /
                   ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) - (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho)) /
                   (V0Y * (1 - rho) + V1Y * rho) +
                   V0X * V0Y * V1Y * (T0X * T1Y - T0Y * T1X) * (-V0Y * (1 - rho) - V1Y * rho) /
                   ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) - (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho)) /
                   (V0Y * (1 - rho) + V1Y * rho) -
                   V0Y * V0Y * V1X * (T0X * T1Y - T0Y * T1X) * (-V0Y * (1 - rho) - V1Y * rho) /
                   ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) - (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho)) /
                   (V0Y * (1 - rho) + V1Y * rho)) * pow(B0Y - B1Y, 2) * pow(V0Y * (1 - rho) + V1Y * rho, 2) - 2 * (V0X *
                                                                                                                   V0Y *
                                                                                                                   V1Y *
                                                                                                                   rho *
                                                                                                                   (T0Y -
                                                                                                                    T1Y) /
                                                                                                                   ((T0Y -
                                                                                                                     T1Y) *
                                                                                                                    (V0X *
                                                                                                                     (1 -
                                                                                                                      rho) +
                                                                                                                     V1X *
                                                                                                                     rho) -
                                                                                                                    (-T0X +
                                                                                                                     T1X) *
                                                                                                                    (-V0Y *
                                                                                                                     (1 -
                                                                                                                      rho) -
                                                                                                                     V1Y *
                                                                                                                     rho)) -
                                                                                                                   V0X *
                                                                                                                   V1Y *
                                                                                                                   V1Y *
                                                                                                                   rho *
                                                                                                                   (T0Y -
                                                                                                                    T1Y) /
                                                                                                                   ((T0Y -
                                                                                                                     T1Y) *
                                                                                                                    (V0X *
                                                                                                                     (1 -
                                                                                                                      rho) +
                                                                                                                     V1X *
                                                                                                                     rho) -
                                                                                                                    (-T0X +
                                                                                                                     T1X) *
                                                                                                                    (-V0Y *
                                                                                                                     (1 -
                                                                                                                      rho) -
                                                                                                                     V1Y *
                                                                                                                     rho)) -
                                                                                                                   V0Y *
                                                                                                                   V0Y *
                                                                                                                   V1X *
                                                                                                                   rho *
                                                                                                                   (T0Y -
                                                                                                                    T1Y) /
                                                                                                                   ((T0Y -
                                                                                                                     T1Y) *
                                                                                                                    (V0X *
                                                                                                                     (1 -
                                                                                                                      rho) +
                                                                                                                     V1X *
                                                                                                                     rho) -
                                                                                                                    (-T0X +
                                                                                                                     T1X) *
                                                                                                                    (-V0Y *
                                                                                                                     (1 -
                                                                                                                      rho) -
                                                                                                                     V1Y *
                                                                                                                     rho)) +
                                                                                                                   V0Y *
                                                                                                                   V1X *
                                                                                                                   V1Y *
                                                                                                                   rho *
                                                                                                                   (T0Y -
                                                                                                                    T1Y) /
                                                                                                                   ((T0Y -
                                                                                                                     T1Y) *
                                                                                                                    (V0X *
                                                                                                                     (1 -
                                                                                                                      rho) +
                                                                                                                     V1X *
                                                                                                                     rho) -
                                                                                                                    (-T0X +
                                                                                                                     T1X) *
                                                                                                                    (-V0Y *
                                                                                                                     (1 -
                                                                                                                      rho) -
                                                                                                                     V1Y *
                                                                                                                     rho)) -
                                                                                                                   V0X *
                                                                                                                   V0Y *
                                                                                                                   V1Y *
                                                                                                                   (T0Y -
                                                                                                                    T1Y) /
                                                                                                                   ((T0Y -
                                                                                                                     T1Y) *
                                                                                                                    (V0X *
                                                                                                                     (1 -
                                                                                                                      rho) +
                                                                                                                     V1X *
                                                                                                                     rho) -
                                                                                                                    (-T0X +
                                                                                                                     T1X) *
                                                                                                                    (-V0Y *
                                                                                                                     (1 -
                                                                                                                      rho) -
                                                                                                                     V1Y *
                                                                                                                     rho)) +
                                                                                                                   V0Y *
                                                                                                                   V0Y *
                                                                                                                   V1X *
                                                                                                                   (T0Y -
                                                                                                                    T1Y) /
                                                                                                                   ((T0Y -
                                                                                                                     T1Y) *
                                                                                                                    (V0X *
                                                                                                                     (1 -
                                                                                                                      rho) +
                                                                                                                     V1X *
                                                                                                                     rho) -
                                                                                                                    (-T0X +
                                                                                                                     T1X) *
                                                                                                                    (-V0Y *
                                                                                                                     (1 -
                                                                                                                      rho) -
                                                                                                                     V1Y *
                                                                                                                     rho))) *
                                                                                                              (B0X *
                                                                                                               B1Y -
                                                                                                               B0Y *
                                                                                                               B1X) *
                                                                                                              (-V0Y *
                                                                                                               (1 -
                                                                                                                rho) -
                                                                                                               V1Y *
                                                                                                               rho) *
                                                                                                              (B0Y -
                                                                                                               B1Y) *
                                                                                                              (V0Y *
                                                                                                               (1 -
                                                                                                                rho) +
                                                                                                               V1Y *
                                                                                                               rho)) *
                pow((B0Y - B1Y) * (V0X * (1 - rho) + V1X * rho) - (-B0X + B1X) * (-V0Y * (1 - rho) - V1Y * rho), -2) *
                pow(V0Y * (1 - rho) + V1Y * rho, -2) / 2 - (2 * (V0X * V0Y * V1Y * rho - V0X * V1Y * V1Y * rho -
                                                                 V0Y * V0Y * V1X * rho + V0Y * V1X * V1Y * rho -
                                                                 V0X * V0Y * V1Y + V0Y * V0Y * V1X) *
                                                            (T0X * T1Y - T0Y * T1X) * (-V0Y * (1 - rho) - V1Y * rho) *
                                                            (T0Y - T1Y) / (V0Y * (1 - rho) + V1Y * rho) *
                                                            pow((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) -
                                                                (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho), -2) *
                                                            (B0Y - B1Y) +
                                                            (V0X * V0Y * V1Y * rho - V0X * V1Y * V1Y * rho -
                                                             V0Y * V0Y * V1X * rho + V0Y * V1X * V1Y * rho -
                                                             V0X * V0Y * V1Y + V0Y * V0Y * V1X) * pow(T0Y - T1Y, 2) /
                                                            (V0Y * (1 - rho) + V1Y * rho) *
                                                            pow((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) -
                                                                (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho), -2) *
                                                            (B0X * B1Y - B0Y * B1X) * (-V0Y * (1 - rho) - V1Y * rho)) /
                                                           ((B0Y - B1Y) * (V0X * (1 - rho) + V1X * rho) -
                                                            (-B0X + B1X) * (-V0Y * (1 - rho) - V1Y * rho)) / 2) *
               (pow(PRX * (V0Y * (1 - rho) + V1Y * rho + PRY) - PRY * (V0X * (1 - rho) + V1X * rho + PRX), 3) *
                pow(V0Y * (1 - rho) + V1Y * rho, -3) +
                pow(PLY, 3) * pow(V0X * (1 - rho) + V1X * rho, 3) * pow(V0Y * (1 - rho) + V1Y * rho, -3)) / 3 + (-(-2 *
                                                                                                                   (-V0X *
                                                                                                                    V0Y *
                                                                                                                    V1Y *
                                                                                                                    rho *
                                                                                                                    (T0X *
                                                                                                                     T1Y -
                                                                                                                     T0Y *
                                                                                                                     T1X) *
                                                                                                                    (-V0Y *
                                                                                                                     (1 -
                                                                                                                      rho) -
                                                                                                                     V1Y *
                                                                                                                     rho) /
                                                                                                                    ((T0Y -
                                                                                                                      T1Y) *
                                                                                                                     (V0X *
                                                                                                                      (1 -
                                                                                                                       rho) +
                                                                                                                      V1X *
                                                                                                                      rho) -
                                                                                                                     (-T0X +
                                                                                                                      T1X) *
                                                                                                                     (-V0Y *
                                                                                                                      (1 -
                                                                                                                       rho) -
                                                                                                                      V1Y *
                                                                                                                      rho)) /
                                                                                                                    (V0Y *
                                                                                                                     (1 -
                                                                                                                      rho) +
                                                                                                                     V1Y *
                                                                                                                     rho) +
                                                                                                                    V0X *
                                                                                                                    V1Y *
                                                                                                                    V1Y *
                                                                                                                    rho *
                                                                                                                    (T0X *
                                                                                                                     T1Y -
                                                                                                                     T0Y *
                                                                                                                     T1X) *
                                                                                                                    (-V0Y *
                                                                                                                     (1 -
                                                                                                                      rho) -
                                                                                                                     V1Y *
                                                                                                                     rho) /
                                                                                                                    ((T0Y -
                                                                                                                      T1Y) *
                                                                                                                     (V0X *
                                                                                                                      (1 -
                                                                                                                       rho) +
                                                                                                                      V1X *
                                                                                                                      rho) -
                                                                                                                     (-T0X +
                                                                                                                      T1X) *
                                                                                                                     (-V0Y *
                                                                                                                      (1 -
                                                                                                                       rho) -
                                                                                                                      V1Y *
                                                                                                                      rho)) /
                                                                                                                    (V0Y *
                                                                                                                     (1 -
                                                                                                                      rho) +
                                                                                                                     V1Y *
                                                                                                                     rho) +
                                                                                                                    V0Y *
                                                                                                                    V0Y *
                                                                                                                    V1X *
                                                                                                                    rho *
                                                                                                                    (T0X *
                                                                                                                     T1Y -
                                                                                                                     T0Y *
                                                                                                                     T1X) *
                                                                                                                    (-V0Y *
                                                                                                                     (1 -
                                                                                                                      rho) -
                                                                                                                     V1Y *
                                                                                                                     rho) /
                                                                                                                    ((T0Y -
                                                                                                                      T1Y) *
                                                                                                                     (V0X *
                                                                                                                      (1 -
                                                                                                                       rho) +
                                                                                                                      V1X *
                                                                                                                      rho) -
                                                                                                                     (-T0X +
                                                                                                                      T1X) *
                                                                                                                     (-V0Y *
                                                                                                                      (1 -
                                                                                                                       rho) -
                                                                                                                      V1Y *
                                                                                                                      rho)) /
                                                                                                                    (V0Y *
                                                                                                                     (1 -
                                                                                                                      rho) +
                                                                                                                     V1Y *
                                                                                                                     rho) -
                                                                                                                    V0Y *
                                                                                                                    V1X *
                                                                                                                    V1Y *
                                                                                                                    rho *
                                                                                                                    (T0X *
                                                                                                                     T1Y -
                                                                                                                     T0Y *
                                                                                                                     T1X) *
                                                                                                                    (-V0Y *
                                                                                                                     (1 -
                                                                                                                      rho) -
                                                                                                                     V1Y *
                                                                                                                     rho) /
                                                                                                                    ((T0Y -
                                                                                                                      T1Y) *
                                                                                                                     (V0X *
                                                                                                                      (1 -
                                                                                                                       rho) +
                                                                                                                      V1X *
                                                                                                                      rho) -
                                                                                                                     (-T0X +
                                                                                                                      T1X) *
                                                                                                                     (-V0Y *
                                                                                                                      (1 -
                                                                                                                       rho) -
                                                                                                                      V1Y *
                                                                                                                      rho)) /
                                                                                                                    (V0Y *
                                                                                                                     (1 -
                                                                                                                      rho) +
                                                                                                                     V1Y *
                                                                                                                     rho) +
                                                                                                                    V0X *
                                                                                                                    V0Y *
                                                                                                                    V1Y *
                                                                                                                    (T0X *
                                                                                                                     T1Y -
                                                                                                                     T0Y *
                                                                                                                     T1X) *
                                                                                                                    (-V0Y *
                                                                                                                     (1 -
                                                                                                                      rho) -
                                                                                                                     V1Y *
                                                                                                                     rho) /
                                                                                                                    ((T0Y -
                                                                                                                      T1Y) *
                                                                                                                     (V0X *
                                                                                                                      (1 -
                                                                                                                       rho) +
                                                                                                                      V1X *
                                                                                                                      rho) -
                                                                                                                     (-T0X +
                                                                                                                      T1X) *
                                                                                                                     (-V0Y *
                                                                                                                      (1 -
                                                                                                                       rho) -
                                                                                                                      V1Y *
                                                                                                                      rho)) /
                                                                                                                    (V0Y *
                                                                                                                     (1 -
                                                                                                                      rho) +
                                                                                                                     V1Y *
                                                                                                                     rho) -
                                                                                                                    V0Y *
                                                                                                                    V0Y *
                                                                                                                    V1X *
                                                                                                                    (T0X *
                                                                                                                     T1Y -
                                                                                                                     T0Y *
                                                                                                                     T1X) *
                                                                                                                    (-V0Y *
                                                                                                                     (1 -
                                                                                                                      rho) -
                                                                                                                     V1Y *
                                                                                                                     rho) /
                                                                                                                    ((T0Y -
                                                                                                                      T1Y) *
                                                                                                                     (V0X *
                                                                                                                      (1 -
                                                                                                                       rho) +
                                                                                                                      V1X *
                                                                                                                      rho) -
                                                                                                                     (-T0X +
                                                                                                                      T1X) *
                                                                                                                     (-V0Y *
                                                                                                                      (1 -
                                                                                                                       rho) -
                                                                                                                      V1Y *
                                                                                                                      rho)) /
                                                                                                                    (V0Y *
                                                                                                                     (1 -
                                                                                                                      rho) +
                                                                                                                     V1Y *
                                                                                                                     rho)) *
                                                                                                                   (B0X *
                                                                                                                    B1Y -
                                                                                                                    B0Y *
                                                                                                                    B1X) *
                                                                                                                   (-V0Y *
                                                                                                                    (1 -
                                                                                                                     rho) -
                                                                                                                    V1Y *
                                                                                                                    rho) *
                                                                                                                   (B0Y -
                                                                                                                    B1Y) *
                                                                                                                   (V0Y *
                                                                                                                    (1 -
                                                                                                                     rho) +
                                                                                                                    V1Y *
                                                                                                                    rho) +
                                                                                                                   (V0X *
                                                                                                                    V0Y *
                                                                                                                    V1Y *
                                                                                                                    rho *
                                                                                                                    (T0Y -
                                                                                                                     T1Y) /
                                                                                                                    ((T0Y -
                                                                                                                      T1Y) *
                                                                                                                     (V0X *
                                                                                                                      (1 -
                                                                                                                       rho) +
                                                                                                                      V1X *
                                                                                                                      rho) -
                                                                                                                     (-T0X +
                                                                                                                      T1X) *
                                                                                                                     (-V0Y *
                                                                                                                      (1 -
                                                                                                                       rho) -
                                                                                                                      V1Y *
                                                                                                                      rho)) -
                                                                                                                    V0X *
                                                                                                                    V1Y *
                                                                                                                    V1Y *
                                                                                                                    rho *
                                                                                                                    (T0Y -
                                                                                                                     T1Y) /
                                                                                                                    ((T0Y -
                                                                                                                      T1Y) *
                                                                                                                     (V0X *
                                                                                                                      (1 -
                                                                                                                       rho) +
                                                                                                                      V1X *
                                                                                                                      rho) -
                                                                                                                     (-T0X +
                                                                                                                      T1X) *
                                                                                                                     (-V0Y *
                                                                                                                      (1 -
                                                                                                                       rho) -
                                                                                                                      V1Y *
                                                                                                                      rho)) -
                                                                                                                    V0Y *
                                                                                                                    V0Y *
                                                                                                                    V1X *
                                                                                                                    rho *
                                                                                                                    (T0Y -
                                                                                                                     T1Y) /
                                                                                                                    ((T0Y -
                                                                                                                      T1Y) *
                                                                                                                     (V0X *
                                                                                                                      (1 -
                                                                                                                       rho) +
                                                                                                                      V1X *
                                                                                                                      rho) -
                                                                                                                     (-T0X +
                                                                                                                      T1X) *
                                                                                                                     (-V0Y *
                                                                                                                      (1 -
                                                                                                                       rho) -
                                                                                                                      V1Y *
                                                                                                                      rho)) +
                                                                                                                    V0Y *
                                                                                                                    V1X *
                                                                                                                    V1Y *
                                                                                                                    rho *
                                                                                                                    (T0Y -
                                                                                                                     T1Y) /
                                                                                                                    ((T0Y -
                                                                                                                      T1Y) *
                                                                                                                     (V0X *
                                                                                                                      (1 -
                                                                                                                       rho) +
                                                                                                                      V1X *
                                                                                                                      rho) -
                                                                                                                     (-T0X +
                                                                                                                      T1X) *
                                                                                                                     (-V0Y *
                                                                                                                      (1 -
                                                                                                                       rho) -
                                                                                                                      V1Y *
                                                                                                                      rho)) -
                                                                                                                    V0X *
                                                                                                                    V0Y *
                                                                                                                    V1Y *
                                                                                                                    (T0Y -
                                                                                                                     T1Y) /
                                                                                                                    ((T0Y -
                                                                                                                      T1Y) *
                                                                                                                     (V0X *
                                                                                                                      (1 -
                                                                                                                       rho) +
                                                                                                                      V1X *
                                                                                                                      rho) -
                                                                                                                     (-T0X +
                                                                                                                      T1X) *
                                                                                                                     (-V0Y *
                                                                                                                      (1 -
                                                                                                                       rho) -
                                                                                                                      V1Y *
                                                                                                                      rho)) +
                                                                                                                    V0Y *
                                                                                                                    V0Y *
                                                                                                                    V1X *
                                                                                                                    (T0Y -
                                                                                                                     T1Y) /
                                                                                                                    ((T0Y -
                                                                                                                      T1Y) *
                                                                                                                     (V0X *
                                                                                                                      (1 -
                                                                                                                       rho) +
                                                                                                                      V1X *
                                                                                                                      rho) -
                                                                                                                     (-T0X +
                                                                                                                      T1X) *
                                                                                                                     (-V0Y *
                                                                                                                      (1 -
                                                                                                                       rho) -
                                                                                                                      V1Y *
                                                                                                                      rho))) *
                                                                                                                   pow(B0X *
                                                                                                                       B1Y -
                                                                                                                       B0Y *
                                                                                                                       B1X,
                                                                                                                       2) *
                                                                                                                   pow(-V0Y *
                                                                                                                       (1 -
                                                                                                                        rho) -
                                                                                                                       V1Y *
                                                                                                                       rho,
                                                                                                                       2)) *
                                                                                                                 pow((B0Y -
                                                                                                                      B1Y) *
                                                                                                                     (V0X *
                                                                                                                      (1 -
                                                                                                                       rho) +
                                                                                                                      V1X *
                                                                                                                      rho) -
                                                                                                                     (-B0X +
                                                                                                                      B1X) *
                                                                                                                     (-V0Y *
                                                                                                                      (1 -
                                                                                                                       rho) -
                                                                                                                      V1Y *
                                                                                                                      rho),
                                                                                                                     -2) *
                                                                                                                 pow(V0Y *
                                                                                                                     (1 -
                                                                                                                      rho) +
                                                                                                                     V1Y *
                                                                                                                     rho,
                                                                                                                     -2) /
                                                                                                                 2 -
                                                                                                                 (-(V0X *
                                                                                                                    V0Y *
                                                                                                                    V1Y *
                                                                                                                    rho -
                                                                                                                    V0X *
                                                                                                                    V1Y *
                                                                                                                    V1Y *
                                                                                                                    rho -
                                                                                                                    V0Y *
                                                                                                                    V0Y *
                                                                                                                    V1X *
                                                                                                                    rho +
                                                                                                                    V0Y *
                                                                                                                    V1X *
                                                                                                                    V1Y *
                                                                                                                    rho -
                                                                                                                    V0X *
                                                                                                                    V0Y *
                                                                                                                    V1Y +
                                                                                                                    V0Y *
                                                                                                                    V0Y *
                                                                                                                    V1X) *
                                                                                                                  pow(T0X *
                                                                                                                      T1Y -
                                                                                                                      T0Y *
                                                                                                                      T1X,
                                                                                                                      2) *
                                                                                                                  pow(-V0Y *
                                                                                                                      (1 -
                                                                                                                       rho) -
                                                                                                                      V1Y *
                                                                                                                      rho,
                                                                                                                      2) *
                                                                                                                  pow((T0Y -
                                                                                                                       T1Y) *
                                                                                                                      (V0X *
                                                                                                                       (1 -
                                                                                                                        rho) +
                                                                                                                       V1X *
                                                                                                                       rho) -
                                                                                                                      (-T0X +
                                                                                                                       T1X) *
                                                                                                                      (-V0Y *
                                                                                                                       (1 -
                                                                                                                        rho) -
                                                                                                                       V1Y *
                                                                                                                       rho),
                                                                                                                      -2) *
                                                                                                                  pow(V0Y *
                                                                                                                      (1 -
                                                                                                                       rho) +
                                                                                                                      V1Y *
                                                                                                                      rho,
                                                                                                                      -2) *
                                                                                                                  (B0Y -
                                                                                                                   B1Y) -
                                                                                                                  2 *
                                                                                                                  (V0X *
                                                                                                                   V0Y *
                                                                                                                   V1Y *
                                                                                                                   rho -
                                                                                                                   V0X *
                                                                                                                   V1Y *
                                                                                                                   V1Y *
                                                                                                                   rho -
                                                                                                                   V0Y *
                                                                                                                   V0Y *
                                                                                                                   V1X *
                                                                                                                   rho +
                                                                                                                   V0Y *
                                                                                                                   V1X *
                                                                                                                   V1Y *
                                                                                                                   rho -
                                                                                                                   V0X *
                                                                                                                   V0Y *
                                                                                                                   V1Y +
                                                                                                                   V0Y *
                                                                                                                   V0Y *
                                                                                                                   V1X) *
                                                                                                                  (T0X *
                                                                                                                   T1Y -
                                                                                                                   T0Y *
                                                                                                                   T1X) *
                                                                                                                  pow(-V0Y *
                                                                                                                      (1 -
                                                                                                                       rho) -
                                                                                                                      V1Y *
                                                                                                                      rho,
                                                                                                                      2) *
                                                                                                                  (T0Y -
                                                                                                                   T1Y) *
                                                                                                                  pow(V0Y *
                                                                                                                      (1 -
                                                                                                                       rho) +
                                                                                                                      V1Y *
                                                                                                                      rho,
                                                                                                                      -2) *
                                                                                                                  pow((T0Y -
                                                                                                                       T1Y) *
                                                                                                                      (V0X *
                                                                                                                       (1 -
                                                                                                                        rho) +
                                                                                                                       V1X *
                                                                                                                       rho) -
                                                                                                                      (-T0X +
                                                                                                                       T1X) *
                                                                                                                      (-V0Y *
                                                                                                                       (1 -
                                                                                                                        rho) -
                                                                                                                       V1Y *
                                                                                                                       rho),
                                                                                                                      -2) *
                                                                                                                  (B0X *
                                                                                                                   B1Y -
                                                                                                                   B0Y *
                                                                                                                   B1X)) /
                                                                                                                 ((B0Y -
                                                                                                                   B1Y) *
                                                                                                                  (V0X *
                                                                                                                   (1 -
                                                                                                                    rho) +
                                                                                                                   V1X *
                                                                                                                   rho) -
                                                                                                                  (-B0X +
                                                                                                                   B1X) *
                                                                                                                  (-V0Y *
                                                                                                                   (1 -
                                                                                                                    rho) -
                                                                                                                   V1Y *
                                                                                                                   rho)) /
                                                                                                                 2) *
                                                                                                                (pow(PRX *
                                                                                                                     (V0Y *
                                                                                                                      (1 -
                                                                                                                       rho) +
                                                                                                                      V1Y *
                                                                                                                      rho +
                                                                                                                      PRY) -
                                                                                                                     PRY *
                                                                                                                     (V0X *
                                                                                                                      (1 -
                                                                                                                       rho) +
                                                                                                                      V1X *
                                                                                                                      rho +
                                                                                                                      PRX),
                                                                                                                     2) *
                                                                                                                 pow(V0Y *
                                                                                                                     (1 -
                                                                                                                      rho) +
                                                                                                                     V1Y *
                                                                                                                     rho,
                                                                                                                     -2) -
                                                                                                                 PLY *
                                                                                                                 PLY *
                                                                                                                 pow(V0X *
                                                                                                                     (1 -
                                                                                                                      rho) +
                                                                                                                     V1X *
                                                                                                                     rho,
                                                                                                                     2) *
                                                                                                                 pow(V0Y *
                                                                                                                     (1 -
                                                                                                                      rho) +
                                                                                                                     V1Y *
                                                                                                                     rho,
                                                                                                                     -2)) /
                                                                                                                2 -
               (-V0X * V0Y * V1Y * rho * (T0X * T1Y - T0Y * T1X) * (-V0Y * (1 - rho) - V1Y * rho) /
                ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) - (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho)) /
                (V0Y * (1 - rho) + V1Y * rho) +
                V0X * V1Y * V1Y * rho * (T0X * T1Y - T0Y * T1X) * (-V0Y * (1 - rho) - V1Y * rho) /
                ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) - (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho)) /
                (V0Y * (1 - rho) + V1Y * rho) +
                V0Y * V0Y * V1X * rho * (T0X * T1Y - T0Y * T1X) * (-V0Y * (1 - rho) - V1Y * rho) /
                ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) - (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho)) /
                (V0Y * (1 - rho) + V1Y * rho) -
                V0Y * V1X * V1Y * rho * (T0X * T1Y - T0Y * T1X) * (-V0Y * (1 - rho) - V1Y * rho) /
                ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) - (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho)) /
                (V0Y * (1 - rho) + V1Y * rho) +
                V0X * V0Y * V1Y * (T0X * T1Y - T0Y * T1X) * (-V0Y * (1 - rho) - V1Y * rho) /
                ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) - (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho)) /
                (V0Y * (1 - rho) + V1Y * rho) -
                V0Y * V0Y * V1X * (T0X * T1Y - T0Y * T1X) * (-V0Y * (1 - rho) - V1Y * rho) /
                ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) - (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho)) /
                (V0Y * (1 - rho) + V1Y * rho)) * pow(B0X * B1Y - B0Y * B1X, 2) * pow(-V0Y * (1 - rho) - V1Y * rho, 2) *
               pow((B0Y - B1Y) * (V0X * (1 - rho) + V1X * rho) - (-B0X + B1X) * (-V0Y * (1 - rho) - V1Y * rho), -2) *
               pow(V0Y * (1 - rho) + V1Y * rho, -2) *
               ((PRX * (V0Y * (1 - rho) + V1Y * rho + PRY) - PRY * (V0X * (1 - rho) + V1X * rho + PRX)) /
                (V0Y * (1 - rho) + V1Y * rho) + PLY * (V0X * (1 - rho) + V1X * rho) / (V0Y * (1 - rho) + V1Y * rho)) /
               2 - (V0X * V0Y * V1Y * rho - V0X * V1Y * V1Y * rho - V0Y * V0Y * V1X * rho + V0Y * V1X * V1Y * rho -
                    V0X * V0Y * V1Y + V0Y * V0Y * V1X) * pow(T0X * T1Y - T0Y * T1X, 2) *
                   pow(-V0Y * (1 - rho) - V1Y * rho, 3) *
                   pow((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) - (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho),
                       -2) * pow(V0Y * (1 - rho) + V1Y * rho, -3) * (B0X * B1Y - B0Y * B1X) /
                   ((B0Y - B1Y) * (V0X * (1 - rho) + V1X * rho) - (-B0X + B1X) * (-V0Y * (1 - rho) - V1Y * rho)) *
                   ((PRX * (V0Y * (1 - rho) + V1Y * rho + PRY) - PRY * (V0X * (1 - rho) + V1X * rho + PRX)) /
                    (V0Y * (1 - rho) + V1Y * rho) +
                    PLY * (V0X * (1 - rho) + V1X * rho) / (V0Y * (1 - rho) + V1Y * rho)) / 2;
    }
};

#endif // VISIBILITY_2_EXAMPLES_NUMERICALINTEGRALOpt_H
