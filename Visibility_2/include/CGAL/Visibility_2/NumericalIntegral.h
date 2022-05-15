#ifndef VISIBILITY_2_EXAMPLES_NUMERICALINTEGRALOpt_H
#define VISIBILITY_2_EXAMPLES_NUMERICALINTEGRALOpt_H

#include <cmath>
#include <iostream>

template<class T>
class Integral {
public:
    static T Compute(T T0X, T T0Y, T T1X, T T1Y, T B0X, T B0Y, T B1X, T B1Y,
                     T PLX, T PLY, T PRX, T PRY, T V0X, T V0Y, T V1X, T V1Y, T rho) {

//        std::cout << "T0: "
//                  << "(" << T0X << "," << T0Y << ")" << std::endl;
//        std::cout << "T1: "
//                  << "(" << T1X << "," << T1Y << ")" << std::endl;
//        std::cout << "B0: "
//                  << "(" << B0X << "," << B0Y << ")" << std::endl;
//        std::cout << "B1: "
//                  << "(" << B1X << "," << B1Y << ")" << std::endl;
//        std::cout << "PL: "
//                  << "(" << PLX << "," << PLX << ")" << std::endl;
//        std::cout << "PR: "
//                  << "(" << PRX << "," << PRY << ")" << std::endl;
//        std::cout << "V0: "
//                  << "(" << V0X << "," << V0Y << ")" << std::endl;
//        std::cout << "V1: "
//                  << "(" << V1X << "," << V1Y << ")" << std::endl;

        T t1 = V0Y * V0X;
        T t2 = V1Y * t1;
        T t3 = (T0Y - T1Y);
        T t5 = 1 - rho;
        T t6 = t5 * V0X;
        T t7 = V1X * rho;
        T t8 = t6 + t7;
        T t11 = t5 * V0Y;
        T t12 = V1Y * rho;
        T t13 = -t11 - t12;
        T t15 = t8 * t3 - t13 * (-T0X + T1X);
        T t16 = 1 / t15;
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
        T t42 = 1 / t41;
        T t44 = t12 * t1;
        T t45 = rho * t20;
        T t46 = rho * t23;
        T t47 = t12 * t25;
        T t48 = -t44 + t45 + t46 - t47 + t2 - t23;
        T t49 = t3 * t3;
        T t50 = t49 * t48;
        T t51 = t15 * t15;
        T t52 = 1 / t51;
        T t53 = t34 * t52;
        T t54 = 1 / t40;
        T t62 = (t11 + t12 + PRY) * PRX - (t6 + t7 + PRX) * PRY;
        T t63 = t62 * t62;
        T t64 = t63 * t63;
        T t65 = t13 * t13;
        T t66 = t65 * t65;
        T t67 = 1 / t66;
        T t69 = PLY * PLY;
        T t70 = t69 * t69;
        T t71 = t8 * t8;
        T t72 = t71 * t71;
        T t80 = T0X * T1Y - T0Y * T1X;
        T t82 = -1 / t13;
        T t84 = t82 * t16 * t13 * t80;
        T t94 = -t82 * t16 * t13 * t80 * t23 + t84 * t2 - t84 * t44 + t84 * t45 + t84 * t46 - t84 * t47;
        T t99 = B0X * B1Y - B0Y * B1X;
        T t102 = -t13 * t13 * t34;
        T t107 = 1 / t65;
        T t109 = t80 * t48;
        T t116 = t99 * t52;
        T t123 = t65 * t13;
        T t124 = -1 / t123;
        T t136 = t99 * t99;
        T t142 = t80 * t80;
        T t143 = t142 * t48;
        T t168 = t82 * t8 * PLY + t82 * t62;
        T t179 = ((-t67 * t72 * t70 + t67 * t64) * (-t42 * t35 * t33 + t54 * t53 * t50)) / 0.8e1 +
                 ((t124 * t71 * t8 * t69 * PLY + t124 * t63 * t62) *
                  (-t107 * t42 * (-2 * t102 * t99 * t33 + t65 * t35 * t94) -
                   t54 * (2 * t53 * t82 * t3 * t13 * t109 + t13 * t116 * t82 * t50))) / 0.6e1 +
                 ((-t107 * t71 * t69 + t107 * t63) * (-t107 * t42 * (-2 * t102 * t99 * t94 + t65 * t136 * t33) - t54 *
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
                                                                                                                  t143))) /
                 0.4e1 - (t168 * t107 * t42 * t65 * t136 * t94) / 0.2e1 -
                 (t168 * t54 * t99 * t124 * t52 * t123 * t143) / 0.2e1;
        return t179;
    };

    static T ComputeUnOpt(T T0X, T T0Y, T T1X, T T1Y, T B0X, T B0Y, T B1X, T B1Y,
                          T PLX, T PLY, T PRX, T PRY, T V0X, T V0Y, T V1X, T V1Y, T rho) {

        return (-(double) (V0X * V0Y * V1Y * rho * (T0Y - T1Y) / ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) -
                                                                  (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho)) -
                           V0X * V1Y * V1Y * rho * (T0Y - T1Y) / ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) -
                                                                  (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho)) -
                           V0Y * V0Y * V1X * rho * (T0Y - T1Y) / ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) -
                                                                  (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho)) +
                           V0Y * V1X * V1Y * rho * (T0Y - T1Y) / ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) -
                                                                  (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho)) -
                           V0X * V0Y * V1Y * (T0Y - T1Y) / ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) -
                                                            (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho)) +
                           V0Y * V0Y * V1X * (T0Y - T1Y) / ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) -
                                                            (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho))) *
                pow(B0Y - B1Y, 0.2e1) * pow((B0Y - B1Y) * (double) (V0X * (1 - rho) + V1X * rho) -
                                            (double) ((-B0X + B1X) * (-V0Y * (1 - rho) - V1Y * rho)), -0.2e1) / 0.2e1 +
                (double) (-V0X * V0Y * V1Y * rho + V0X * V1Y * V1Y * rho + V0Y * V0Y * V1X * rho -
                          V0Y * V1X * V1Y * rho + V0X * V0Y * V1Y - V0Y * V0Y * V1X) *
                (double) (int) pow((double) (T0Y - T1Y), (double) 2) * (double) (int) pow(
                        (double) ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) -
                                  (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho)), (double) (-2)) * (B0Y - B1Y) /
                ((B0Y - B1Y) * (double) (V0X * (1 - rho) + V1X * rho) -
                 (double) ((-B0X + B1X) * (-V0Y * (1 - rho) - V1Y * rho))) / 0.2e1) * ((double) ((int) pow(
                (double) (PRX * (V0Y * (1 - rho) + V1Y * rho + PRY) - PRY * (V0X * (1 - rho) + V1X * rho + PRX)),
                (double) 4) * (int) pow((double) (V0Y * (1 - rho) + V1Y * rho), (double) (-4))) - pow(PLY, 0.4e1) *
                                                                                                  (double) (int) pow(
                                                                                                          (double) (
                                                                                                                  V0X *
                                                                                                                  (1 -
                                                                                                                   rho) +
                                                                                                                  V1X *
                                                                                                                  rho),
                                                                                                          (double) 4) *
                                                                                                  (double) (int) pow(
                                                                                                          (double) (
                                                                                                                  V0Y *
                                                                                                                  (1 -
                                                                                                                   rho) +
                                                                                                                  V1Y *
                                                                                                                  rho),
                                                                                                          (double) (-4))) /
               0.4e1 + (-((double) (-V0X * V0Y * V1Y * rho * (T0X * T1Y - T0Y * T1X) * (-V0Y * (1 - rho) - V1Y * rho) /
                                    ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) -
                                     (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho)) / (V0Y * (1 - rho) + V1Y * rho) +
                                    V0X * V1Y * V1Y * rho * (T0X * T1Y - T0Y * T1X) * (-V0Y * (1 - rho) - V1Y * rho) /
                                    ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) -
                                     (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho)) / (V0Y * (1 - rho) + V1Y * rho) +
                                    V0Y * V0Y * V1X * rho * (T0X * T1Y - T0Y * T1X) * (-V0Y * (1 - rho) - V1Y * rho) /
                                    ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) -
                                     (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho)) / (V0Y * (1 - rho) + V1Y * rho) -
                                    V0Y * V1X * V1Y * rho * (T0X * T1Y - T0Y * T1X) * (-V0Y * (1 - rho) - V1Y * rho) /
                                    ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) -
                                     (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho)) / (V0Y * (1 - rho) + V1Y * rho) +
                                    V0X * V0Y * V1Y * (T0X * T1Y - T0Y * T1X) * (-V0Y * (1 - rho) - V1Y * rho) /
                                    ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) -
                                     (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho)) / (V0Y * (1 - rho) + V1Y * rho) -
                                    V0Y * V0Y * V1X * (T0X * T1Y - T0Y * T1X) * (-V0Y * (1 - rho) - V1Y * rho) /
                                    ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) -
                                     (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho)) / (V0Y * (1 - rho) + V1Y * rho)) *
                          pow(B0Y - B1Y, 0.2e1) *
                          (double) (int) pow((double) (V0Y * (1 - rho) + V1Y * rho), (double) 2) - 0.2e1 * (double) (
                V0X * V0Y * V1Y * rho * (T0Y - T1Y) /
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
                                                                                                   ((double) B0X * B1Y -
                                                                                                    B0Y *
                                                                                                    (double) B1X) *
                                                                                                   (double) (-V0Y *
                                                                                                             (1 - rho) -
                                                                                                             V1Y *
                                                                                                             rho) *
                                                                                                   (B0Y - B1Y) *
                                                                                                   (double) (V0Y *
                                                                                                             (1 - rho) +
                                                                                                             V1Y *
                                                                                                             rho)) *
                        pow((B0Y - B1Y) * (double) (V0X * (1 - rho) + V1X * rho) -
                            (double) ((-B0X + B1X) * (-V0Y * (1 - rho) - V1Y * rho)), -0.2e1) *
                        (double) (int) pow((double) (V0Y * (1 - rho) + V1Y * rho), (double) (-2)) / 0.2e1 - (0.2e1 *
                                                                                                             (double) (
                                                                                                                     -V0X *
                                                                                                                     V0Y *
                                                                                                                     V1Y *
                                                                                                                     rho +
                                                                                                                     V0X *
                                                                                                                     V1Y *
                                                                                                                     V1Y *
                                                                                                                     rho +
                                                                                                                     V0Y *
                                                                                                                     V0Y *
                                                                                                                     V1X *
                                                                                                                     rho -
                                                                                                                     V0Y *
                                                                                                                     V1X *
                                                                                                                     V1Y *
                                                                                                                     rho +
                                                                                                                     V0X *
                                                                                                                     V0Y *
                                                                                                                     V1Y -
                                                                                                                     V0Y *
                                                                                                                     V0Y *
                                                                                                                     V1X) *
                                                                                                             (double) (
                                                                                                                     T0X *
                                                                                                                     T1Y -
                                                                                                                     T0Y *
                                                                                                                     T1X) *
                                                                                                             (double) (
                                                                                                                     -V0Y *
                                                                                                                     (1 -
                                                                                                                      rho) -
                                                                                                                     V1Y *
                                                                                                                     rho) *
                                                                                                             (double) (
                                                                                                                     T0Y -
                                                                                                                     T1Y) /
                                                                                                             (double) (
                                                                                                                     V0Y *
                                                                                                                     (1 -
                                                                                                                      rho) +
                                                                                                                     V1Y *
                                                                                                                     rho) *
                                                                                                             (double) (int) pow(
                                                                                                                     (double) (
                                                                                                                             (T0Y -
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
                                                                                                                              rho)),
                                                                                                                     (double) (-2)) *
                                                                                                             (B0Y -
                                                                                                              B1Y) +
                                                                                                             (double) (
                                                                                                                     -V0X *
                                                                                                                     V0Y *
                                                                                                                     V1Y *
                                                                                                                     rho +
                                                                                                                     V0X *
                                                                                                                     V1Y *
                                                                                                                     V1Y *
                                                                                                                     rho +
                                                                                                                     V0Y *
                                                                                                                     V0Y *
                                                                                                                     V1X *
                                                                                                                     rho -
                                                                                                                     V0Y *
                                                                                                                     V1X *
                                                                                                                     V1Y *
                                                                                                                     rho +
                                                                                                                     V0X *
                                                                                                                     V0Y *
                                                                                                                     V1Y -
                                                                                                                     V0Y *
                                                                                                                     V0Y *
                                                                                                                     V1X) *
                                                                                                             (double) (int) pow(
                                                                                                                     (double) (
                                                                                                                             T0Y -
                                                                                                                             T1Y),
                                                                                                                     (double) 2) /
                                                                                                             (double) (
                                                                                                                     V0Y *
                                                                                                                     (1 -
                                                                                                                      rho) +
                                                                                                                     V1Y *
                                                                                                                     rho) *
                                                                                                             (double) (int) pow(
                                                                                                                     (double) (
                                                                                                                             (T0Y -
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
                                                                                                                              rho)),
                                                                                                                     (double) (-2)) *
                                                                                                             ((double) B0X *
                                                                                                              B1Y -
                                                                                                              B0Y *
                                                                                                              (double) B1X) *
                                                                                                             (double) (
                                                                                                                     -V0Y *
                                                                                                                     (1 -
                                                                                                                      rho) -
                                                                                                                     V1Y *
                                                                                                                     rho)) /
                                                                                                            ((B0Y -
                                                                                                              B1Y) *
                                                                                                             (double) (
                                                                                                                     V0X *
                                                                                                                     (1 -
                                                                                                                      rho) +
                                                                                                                     V1X *
                                                                                                                     rho) -
                                                                                                             (double) (
                                                                                                                     (-B0X +
                                                                                                                      B1X) *
                                                                                                                     (-V0Y *
                                                                                                                      (1 -
                                                                                                                       rho) -
                                                                                                                      V1Y *
                                                                                                                      rho))) /
                                                                                                            0.2e1) *
                       ((double) ((int) pow((double) (PRX * (V0Y * (1 - rho) + V1Y * rho + PRY) -
                                                      PRY * (V0X * (1 - rho) + V1X * rho + PRX)), (double) 3) *
                                  (int) pow((double) (V0Y * (1 - rho) + V1Y * rho), (double) (-3))) +
                        pow(PLY, 0.3e1) * (double) (int) pow((double) (V0X * (1 - rho) + V1X * rho), (double) 3) *
                        (double) (int) pow((double) (V0Y * (1 - rho) + V1Y * rho), (double) (-3))) / 0.3e1 + (-(-0.2e1 *
                                                                                                                (double) (
                                                                                                                        -V0X *
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
                                                                                                                ((double) B0X *
                                                                                                                 B1Y -
                                                                                                                 B0Y *
                                                                                                                 (double) B1X) *
                                                                                                                (double) (
                                                                                                                        -V0Y *
                                                                                                                        (1 -
                                                                                                                         rho) -
                                                                                                                        V1Y *
                                                                                                                        rho) *
                                                                                                                (B0Y -
                                                                                                                 B1Y) *
                                                                                                                (double) (
                                                                                                                        V0Y *
                                                                                                                        (1 -
                                                                                                                         rho) +
                                                                                                                        V1Y *
                                                                                                                        rho) +
                                                                                                                (double) (
                                                                                                                        V0X *
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
                                                                                                                pow((double) B0X *
                                                                                                                    B1Y -
                                                                                                                    B0Y *
                                                                                                                    (double) B1X,
                                                                                                                    0.2e1) *
                                                                                                                (double) (int) pow(
                                                                                                                        (double) (
                                                                                                                                -V0Y *
                                                                                                                                (1 -
                                                                                                                                 rho) -
                                                                                                                                V1Y *
                                                                                                                                rho),
                                                                                                                        (double) 2)) *
                                                                                                              pow((B0Y -
                                                                                                                   B1Y) *
                                                                                                                  (double) (
                                                                                                                          V0X *
                                                                                                                          (1 -
                                                                                                                           rho) +
                                                                                                                          V1X *
                                                                                                                          rho) -
                                                                                                                  (double) (
                                                                                                                          (-B0X +
                                                                                                                           B1X) *
                                                                                                                          (-V0Y *
                                                                                                                           (1 -
                                                                                                                            rho) -
                                                                                                                           V1Y *
                                                                                                                           rho)),
                                                                                                                  -0.2e1) *
                                                                                                              (double) (int) pow(
                                                                                                                      (double) (
                                                                                                                              V0Y *
                                                                                                                              (1 -
                                                                                                                               rho) +
                                                                                                                              V1Y *
                                                                                                                              rho),
                                                                                                                      (double) (-2)) /
                                                                                                              0.2e1 -
                                                                                                              (-(double) (
                                                                                                                      -V0X *
                                                                                                                      V0Y *
                                                                                                                      V1Y *
                                                                                                                      rho +
                                                                                                                      V0X *
                                                                                                                      V1Y *
                                                                                                                      V1Y *
                                                                                                                      rho +
                                                                                                                      V0Y *
                                                                                                                      V0Y *
                                                                                                                      V1X *
                                                                                                                      rho -
                                                                                                                      V0Y *
                                                                                                                      V1X *
                                                                                                                      V1Y *
                                                                                                                      rho +
                                                                                                                      V0X *
                                                                                                                      V0Y *
                                                                                                                      V1Y -
                                                                                                                      V0Y *
                                                                                                                      V0Y *
                                                                                                                      V1X) *
                                                                                                               (double) (int) pow(
                                                                                                                       (double) (
                                                                                                                               T0X *
                                                                                                                               T1Y -
                                                                                                                               T0Y *
                                                                                                                               T1X),
                                                                                                                       (double) 2) *
                                                                                                               (double) (int) pow(
                                                                                                                       (double) (
                                                                                                                               -V0Y *
                                                                                                                               (1 -
                                                                                                                                rho) -
                                                                                                                               V1Y *
                                                                                                                               rho),
                                                                                                                       (double) 2) *
                                                                                                               (double) (int) pow(
                                                                                                                       (double) (
                                                                                                                               (T0Y -
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
                                                                                                                                rho)),
                                                                                                                       (double) (-2)) *
                                                                                                               (double) (int) pow(
                                                                                                                       (double) (
                                                                                                                               V0Y *
                                                                                                                               (1 -
                                                                                                                                rho) +
                                                                                                                               V1Y *
                                                                                                                               rho),
                                                                                                                       (double) (-2)) *
                                                                                                               (B0Y -
                                                                                                                B1Y) -
                                                                                                               0.2e1 *
                                                                                                               (double) (
                                                                                                                       -V0X *
                                                                                                                       V0Y *
                                                                                                                       V1Y *
                                                                                                                       rho +
                                                                                                                       V0X *
                                                                                                                       V1Y *
                                                                                                                       V1Y *
                                                                                                                       rho +
                                                                                                                       V0Y *
                                                                                                                       V0Y *
                                                                                                                       V1X *
                                                                                                                       rho -
                                                                                                                       V0Y *
                                                                                                                       V1X *
                                                                                                                       V1Y *
                                                                                                                       rho +
                                                                                                                       V0X *
                                                                                                                       V0Y *
                                                                                                                       V1Y -
                                                                                                                       V0Y *
                                                                                                                       V0Y *
                                                                                                                       V1X) *
                                                                                                               (double) (
                                                                                                                       T0X *
                                                                                                                       T1Y -
                                                                                                                       T0Y *
                                                                                                                       T1X) *
                                                                                                               (double) (int) pow(
                                                                                                                       (double) (
                                                                                                                               -V0Y *
                                                                                                                               (1 -
                                                                                                                                rho) -
                                                                                                                               V1Y *
                                                                                                                               rho),
                                                                                                                       (double) 2) *
                                                                                                               (double) (
                                                                                                                       T0Y -
                                                                                                                       T1Y) *
                                                                                                               (double) (int) pow(
                                                                                                                       (double) (
                                                                                                                               V0Y *
                                                                                                                               (1 -
                                                                                                                                rho) +
                                                                                                                               V1Y *
                                                                                                                               rho),
                                                                                                                       (double) (-2)) *
                                                                                                               (double) (int) pow(
                                                                                                                       (double) (
                                                                                                                               (T0Y -
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
                                                                                                                                rho)),
                                                                                                                       (double) (-2)) *
                                                                                                               ((double) B0X *
                                                                                                                B1Y -
                                                                                                                B0Y *
                                                                                                                (double) B1X)) /
                                                                                                              ((B0Y -
                                                                                                                B1Y) *
                                                                                                               (double) (
                                                                                                                       V0X *
                                                                                                                       (1 -
                                                                                                                        rho) +
                                                                                                                       V1X *
                                                                                                                       rho) -
                                                                                                               (double) (
                                                                                                                       (-B0X +
                                                                                                                        B1X) *
                                                                                                                       (-V0Y *
                                                                                                                        (1 -
                                                                                                                         rho) -
                                                                                                                        V1Y *
                                                                                                                        rho))) /
                                                                                                              0.2e1) *
                                                                                                             ((double) (
                                                                                                                     (int) pow(
                                                                                                                             (double) (
                                                                                                                                     PRX *
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
                                                                                                                                      PRX)),
                                                                                                                             (double) 2) *
                                                                                                                     (int) pow(
                                                                                                                             (double) (
                                                                                                                                     V0Y *
                                                                                                                                     (1 -
                                                                                                                                      rho) +
                                                                                                                                     V1Y *
                                                                                                                                     rho),
                                                                                                                             (double) (-2))) -
                                                                                                              PLY *
                                                                                                              PLY *
                                                                                                              (double) (int) pow(
                                                                                                                      (double) (
                                                                                                                              V0X *
                                                                                                                              (1 -
                                                                                                                               rho) +
                                                                                                                              V1X *
                                                                                                                              rho),
                                                                                                                      (double) 2) *
                                                                                                              (double) (int) pow(
                                                                                                                      (double) (
                                                                                                                              V0Y *
                                                                                                                              (1 -
                                                                                                                               rho) +
                                                                                                                              V1Y *
                                                                                                                              rho),
                                                                                                                      (double) (-2))) /
                                                                                                             0.2e1 -
               (double) (-V0X * V0Y * V1Y * rho * (T0X * T1Y - T0Y * T1X) * (-V0Y * (1 - rho) - V1Y * rho) /
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
                         (V0Y * (1 - rho) + V1Y * rho)) * pow((double) B0X * B1Y - B0Y * (double) B1X, 0.2e1) *
               (double) (int) pow((double) (-V0Y * (1 - rho) - V1Y * rho), (double) 2) *
               pow((B0Y - B1Y) * (double) (V0X * (1 - rho) + V1X * rho) -
                   (double) ((-B0X + B1X) * (-V0Y * (1 - rho) - V1Y * rho)), -0.2e1) *
               (double) (int) pow((double) (V0Y * (1 - rho) + V1Y * rho), (double) (-2)) *
               ((double) ((PRX * (V0Y * (1 - rho) + V1Y * rho + PRY) - PRY * (V0X * (1 - rho) + V1X * rho + PRX)) /
                          (V0Y * (1 - rho) + V1Y * rho)) +
                PLY * (double) (V0X * (1 - rho) + V1X * rho) / (double) (V0Y * (1 - rho) + V1Y * rho)) / 0.2e1 -
               (double) (-V0X * V0Y * V1Y * rho + V0X * V1Y * V1Y * rho + V0Y * V0Y * V1X * rho -
                         V0Y * V1X * V1Y * rho + V0X * V0Y * V1Y - V0Y * V0Y * V1X) *
               (double) (int) pow((double) (T0X * T1Y - T0Y * T1X), (double) 2) *
               (double) (int) pow((double) (-V0Y * (1 - rho) - V1Y * rho), (double) 3) * (double) (int) pow(
                       (double) ((T0Y - T1Y) * (V0X * (1 - rho) + V1X * rho) -
                                 (-T0X + T1X) * (-V0Y * (1 - rho) - V1Y * rho)), (double) (-2)) *
               (double) (int) pow((double) (V0Y * (1 - rho) + V1Y * rho), (double) (-3)) *
               ((double) B0X * B1Y - B0Y * (double) B1X) / ((B0Y - B1Y) * (double) (V0X * (1 - rho) + V1X * rho) -
                                                            (double) ((-B0X + B1X) * (-V0Y * (1 - rho) - V1Y * rho))) *
               ((double) ((PRX * (V0Y * (1 - rho) + V1Y * rho + PRY) - PRY * (V0X * (1 - rho) + V1X * rho + PRX)) /
                          (V0Y * (1 - rho) + V1Y * rho)) +
                PLY * (double) (V0X * (1 - rho) + V1X * rho) / (double) (V0Y * (1 - rho) + V1Y * rho)) / 0.2e1;

    }


    static T ComputeOrig(T T0X, T T0Y, T T1X, T T1Y, T B0X, T B0Y, T B1X, T B1Y,
                         T PLX, T PLY, T PRX, T PRY, T V0X, T V0Y, T V1X, T V1Y, T phi) {


        T cg = -sin(phi) / (0.2e1 * ((T0X - T0Y - T1X + T1Y) * B0X + (-T0X + T0Y + T1X - T1Y) * B1X -
                                     (T0X + T0Y - T1X - T1Y) * (B0Y - B1Y)) *
                            ((T0X + T0Y - T1X - T1Y) * B0X + (-T0X - T0Y + T1X + T1Y) * B1X +
                             (T0X - T0Y - T1X + T1Y) * (B0Y - B1Y)) * pow(cos(phi), 0.4e1) +
                            0.4e1 * ((T0X - T1X) * B0X + (-T0X + T1X) * B1X - (T0Y - T1Y) * (B0Y - B1Y)) * sin(phi) *
                            ((T0Y - T1Y) * B0X + (-T0Y + T1Y) * B1X + (T0X - T1X) * (B0Y - B1Y)) *
                            pow(cos(phi), 0.3e1) +
                            ((-0.4e1 * T0X * T0X + 0.8e1 * T0X * T1X - 0.4e1 * T1X * T1X +
                              0.2e1 * pow(T0Y - T1Y, 0.2e1)) * B0X * B0X +
                             ((0.8e1 * T0X * T0X - 0.16e2 * T0X * T1X + 0.8e1 * T1X * T1X -
                               0.4e1 * pow(T0Y - T1Y, 0.2e1)) * B1X + 0.8e1 * (T0Y - T1Y) * (T0X - T1X) * (B0Y - B1Y)) *
                             B0X + (-0.4e1 * T0X * T0X + 0.8e1 * T0X * T1X - 0.4e1 * T1X * T1X +
                                    0.2e1 * pow(T0Y - T1Y, 0.2e1)) * B1X * B1X -
                             0.8e1 * (T0Y - T1Y) * (T0X - T1X) * (B0Y - B1Y) * B1X +
                             0.2e1 * pow(T0X - T1X, 0.2e1) * pow(B0Y - B1Y, 0.2e1)) * pow(cos(phi), 0.2e1) -
                            0.4e1 * (T0X - T1X) * (B0X - B1X) * sin(phi) *
                            ((T0Y - T1Y) * B0X + (-T0Y + T1Y) * B1X + (T0X - T1X) * (B0Y - B1Y)) * cos(phi) +
                            0.2e1 * pow(T0X - T1X, 0.2e1) * pow(B0X - B1X, 0.2e1)) *
               (((T0X + T0Y - T1X - T1Y) * (T0X - T0Y - T1X + T1Y) * pow(cos(phi), 0.2e1) +
                 0.2e1 * sin(phi) * (T0Y - T1Y) * (T0X - T1X) * cos(phi) - pow(T0X - T1X, 0.2e1)) *
                sqrt(-pow(T0X * T1Y - T0Y * T1X + (PRX - PRY / tan(phi)) * T0Y - (PRX - PRY / tan(phi)) * T1Y, 0.2e1) /
                     (0.2e1 * sin(phi) * cos(phi) * T0X * T0Y - 0.2e1 * sin(phi) * cos(phi) * T0X * T1Y -
                      0.2e1 * sin(phi) * cos(phi) * T0Y * T1X + 0.2e1 * sin(phi) * cos(phi) * T1X * T1Y +
                      T0X * T0X * pow(cos(phi), 0.2e1) - 0.2e1 * T0X * T1X * pow(cos(phi), 0.2e1) -
                      pow(cos(phi), 0.2e1) * T0Y * T0Y + 0.2e1 * pow(cos(phi), 0.2e1) * T0Y * T1Y +
                      T1X * T1X * pow(cos(phi), 0.2e1) - pow(cos(phi), 0.2e1) * T1Y * T1Y - T0X * T0X +
                      0.2e1 * T0X * T1X - T1X * T1X)) * (PRX - PRY / tan(phi)) *
                (-0.24e2 * B0X * B0Y * B1X * B1Y * T0X * T1Y + 0.24e2 * B0X * B0Y * B1X * B1Y * T0Y * T1X -
                 0.3e1 * B0Y * B0Y * T1Y * pow(PRX - PRY / tan(phi), 0.3e1) -
                 0.3e1 * B1Y * B1Y * T1Y * pow(PRX - PRY / tan(phi), 0.3e1) +
                 0.3e1 * B1Y * B1Y * T0Y * pow(PRX - PRY / tan(phi), 0.3e1) +
                 0.3e1 * B0Y * B0Y * T0Y * pow(PRX - PRY / tan(phi), 0.3e1) +
                 0.6e1 * B0X * B0X * B1Y * B1Y * T0Y * (PRX - PRY / tan(phi)) -
                 0.6e1 * B0X * B0X * B1Y * B1Y * T1Y * (PRX - PRY / tan(phi)) -
                 0.8e1 * B0X * B1Y * B1Y * T0Y * pow(PRX - PRY / tan(phi), 0.2e1) +
                 0.8e1 * B0X * B1Y * B1Y * T1Y * pow(PRX - PRY / tan(phi), 0.2e1) +
                 0.6e1 * B0Y * B0Y * B1X * B1X * T0Y * (PRX - PRY / tan(phi)) -
                 0.6e1 * B0Y * B0Y * B1X * B1X * T1Y * (PRX - PRY / tan(phi)) -
                 0.8e1 * B0Y * B0Y * B1X * T0Y * pow(PRX - PRY / tan(phi), 0.2e1) +
                 0.8e1 * B0Y * B0Y * B1X * T1Y * pow(PRX - PRY / tan(phi), 0.2e1) +
                 0.4e1 * B0Y * B0Y * T0X * T1Y * pow(PRX - PRY / tan(phi), 0.2e1) -
                 0.4e1 * B0Y * B0Y * T0Y * T1X * pow(PRX - PRY / tan(phi), 0.2e1) -
                 0.6e1 * B0Y * B1Y * T0Y * pow(PRX - PRY / tan(phi), 0.3e1) +
                 0.6e1 * B0Y * B1Y * T1Y * pow(PRX - PRY / tan(phi), 0.3e1) +
                 0.4e1 * B1Y * B1Y * T0X * T1Y * pow(PRX - PRY / tan(phi), 0.2e1) -
                 0.4e1 * B1Y * B1Y * T0Y * T1X * pow(PRX - PRY / tan(phi), 0.2e1) -
                 0.12e2 * B0X * B0Y * B1X * B1Y * T0Y * (PRX - PRY / tan(phi)) +
                 0.12e2 * B0X * B0Y * B1X * B1Y * T1Y * (PRX - PRY / tan(phi)) +
                 0.12e2 * B0X * B0Y * B1Y * T0X * T1Y * (PRX - PRY / tan(phi)) -
                 0.12e2 * B0X * B0Y * B1Y * T0Y * T1X * (PRX - PRY / tan(phi)) +
                 0.12e2 * B0Y * B1X * B1Y * T0X * T1Y * (PRX - PRY / tan(phi)) -
                 0.12e2 * B0Y * B1X * B1Y * T0Y * T1X * (PRX - PRY / tan(phi)) +
                 0.8e1 * B0X * B0Y * B1Y * T0Y * pow(PRX - PRY / tan(phi), 0.2e1) -
                 0.8e1 * B0X * B0Y * B1Y * T1Y * pow(PRX - PRY / tan(phi), 0.2e1) -
                 0.12e2 * B0X * B1Y * B1Y * T0X * T1Y * (PRX - PRY / tan(phi)) +
                 0.12e2 * B0X * B1Y * B1Y * T0Y * T1X * (PRX - PRY / tan(phi)) -
                 0.12e2 * B0Y * B0Y * B1X * T0X * T1Y * (PRX - PRY / tan(phi)) +
                 0.12e2 * B0Y * B0Y * B1X * T0Y * T1X * (PRX - PRY / tan(phi)) +
                 0.8e1 * B0Y * B1X * B1Y * T0Y * pow(PRX - PRY / tan(phi), 0.2e1) -
                 0.8e1 * B0Y * B1X * B1Y * T1Y * pow(PRX - PRY / tan(phi), 0.2e1) -
                 0.8e1 * B0Y * B1Y * T0X * T1Y * pow(PRX - PRY / tan(phi), 0.2e1) +
                 0.8e1 * B0Y * B1Y * T0Y * T1X * pow(PRX - PRY / tan(phi), 0.2e1) +
                 0.12e2 * B0X * B0X * B1Y * B1Y * T0X * T1Y - 0.12e2 * B0X * B0X * B1Y * B1Y * T0Y * T1X +
                 0.12e2 * B0Y * B0Y * B1X * B1X * T0X * T1Y - 0.12e2 * B0Y * B0Y * B1X * B1X * T0Y * T1X) /
                (T0X * T1Y - T0Y * T1X + (PRX - PRY / tan(phi)) * T0Y - (PRX - PRY / tan(phi)) * T1Y) / 0.12e2 +
                ((B0X + B0Y - B1X - B1Y) * (B0X - B0Y - B1X + B1Y) * pow(cos(phi), 0.2e1) +
                 0.2e1 * sin(phi) * (B0Y - B1Y) * (B0X - B1X) * cos(phi) - pow(B0X - B1X, 0.2e1)) *
                sqrt(-pow(B0X * B1Y - B0Y * B1X + (PRX - PRY / tan(phi)) * B0Y - (PRX - PRY / tan(phi)) * B1Y, 0.2e1) /
                     (B0X * B0X * pow(cos(phi), 0.2e1) + 0.2e1 * B0X * B0Y * sin(phi) * cos(phi) -
                      0.2e1 * B0X * B1X * pow(cos(phi), 0.2e1) - 0.2e1 * B0X * B1Y * sin(phi) * cos(phi) -
                      B0Y * B0Y * pow(cos(phi), 0.2e1) - 0.2e1 * B0Y * B1X * sin(phi) * cos(phi) +
                      0.2e1 * B0Y * B1Y * pow(cos(phi), 0.2e1) + B1X * B1X * pow(cos(phi), 0.2e1) +
                      0.2e1 * B1X * B1Y * sin(phi) * cos(phi) - B1Y * B1Y * pow(cos(phi), 0.2e1) - B0X * B0X +
                      0.2e1 * B0X * B1X - B1X * B1X)) * (PRX - PRY / tan(phi)) *
                (-0.24e2 * B0X * B1Y * T0X * T0Y * T1X * T1Y + 0.24e2 * B0Y * B1X * T0X * T0Y * T1X * T1Y +
                 0.3e1 * B0Y * T1Y * T1Y * pow(PRX - PRY / tan(phi), 0.3e1) -
                 0.3e1 * B1Y * T0Y * T0Y * pow(PRX - PRY / tan(phi), 0.3e1) -
                 0.3e1 * B1Y * T1Y * T1Y * pow(PRX - PRY / tan(phi), 0.3e1) -
                 0.8e1 * B0Y * T0X * T1Y * T1Y * pow(PRX - PRY / tan(phi), 0.2e1) +
                 0.6e1 * B0Y * T0Y * T0Y * T1X * T1X * (PRX - PRY / tan(phi)) -
                 0.8e1 * B0Y * T0Y * T0Y * T1X * pow(PRX - PRY / tan(phi), 0.2e1) -
                 0.6e1 * B0Y * T0Y * T1Y * pow(PRX - PRY / tan(phi), 0.3e1) -
                 0.6e1 * B1Y * T0X * T0X * T1Y * T1Y * (PRX - PRY / tan(phi)) +
                 0.8e1 * B1Y * T0X * T1Y * T1Y * pow(PRX - PRY / tan(phi), 0.2e1) -
                 0.6e1 * B1Y * T0Y * T0Y * T1X * T1X * (PRX - PRY / tan(phi)) +
                 0.8e1 * B1Y * T0Y * T0Y * T1X * pow(PRX - PRY / tan(phi), 0.2e1) +
                 0.6e1 * B1Y * T0Y * T1Y * pow(PRX - PRY / tan(phi), 0.3e1) +
                 0.3e1 * B0Y * T0Y * T0Y * pow(PRX - PRY / tan(phi), 0.3e1) +
                 0.12e2 * B0X * B1Y * T0X * T0X * T1Y * T1Y +
                 0.12e2 * B0X * B1Y * T0Y * T0Y * T1X * T1X - 0.12e2 * B0Y * B1X * T0X * T0X * T1Y * T1Y -
                 0.12e2 * B0Y * B1X * T0Y * T0Y * T1X * T1X +
                 0.4e1 * B0X * B1Y * T0Y * T0Y * pow(PRX - PRY / tan(phi), 0.2e1) +
                 0.4e1 * B0X * B1Y * T1Y * T1Y * pow(PRX - PRY / tan(phi), 0.2e1) -
                 0.4e1 * B0Y * B1X * T0Y * T0Y * pow(PRX - PRY / tan(phi), 0.2e1) -
                 0.4e1 * B0Y * B1X * T1Y * T1Y * pow(PRX - PRY / tan(phi), 0.2e1) +
                 0.6e1 * B0Y * T0X * T0X * T1Y * T1Y * (PRX - PRY / tan(phi)) +
                 0.12e2 * B0X * B1Y * T0X * T0Y * T1Y * (PRX - PRY / tan(phi)) +
                 0.12e2 * B0X * B1Y * T0Y * T1X * T1Y * (PRX - PRY / tan(phi)) -
                 0.12e2 * B0Y * B1X * T0X * T0Y * T1Y * (PRX - PRY / tan(phi)) -
                 0.12e2 * B0Y * B1X * T0Y * T1X * T1Y * (PRX - PRY / tan(phi)) -
                 0.12e2 * B0Y * T0X * T0Y * T1X * T1Y * (PRX - PRY / tan(phi)) +
                 0.12e2 * B1Y * T0X * T0Y * T1X * T1Y * (PRX - PRY / tan(phi)) -
                 0.12e2 * B0X * B1Y * T0X * T1Y * T1Y * (PRX - PRY / tan(phi)) -
                 0.12e2 * B0X * B1Y * T0Y * T0Y * T1X * (PRX - PRY / tan(phi)) -
                 0.8e1 * B0X * B1Y * T0Y * T1Y * pow(PRX - PRY / tan(phi), 0.2e1) +
                 0.12e2 * B0Y * B1X * T0X * T1Y * T1Y * (PRX - PRY / tan(phi)) +
                 0.12e2 * B0Y * B1X * T0Y * T0Y * T1X * (PRX - PRY / tan(phi)) +
                 0.8e1 * B0Y * B1X * T0Y * T1Y * pow(PRX - PRY / tan(phi), 0.2e1) +
                 0.8e1 * B0Y * T0X * T0Y * T1Y * pow(PRX - PRY / tan(phi), 0.2e1) +
                 0.8e1 * B0Y * T0Y * T1X * T1Y * pow(PRX - PRY / tan(phi), 0.2e1) -
                 0.8e1 * B1Y * T0X * T0Y * T1Y * pow(PRX - PRY / tan(phi), 0.2e1) -
                 0.8e1 * B1Y * T0Y * T1X * T1Y * pow(PRX - PRY / tan(phi), 0.2e1)) /
                (B0X * B1Y - B0Y * B1X + (PRX - PRY / tan(phi)) * B0Y - (PRX - PRY / tan(phi)) * B1Y) / 0.12e2) +
               sin(phi) / (0.2e1 * ((T0X - T0Y - T1X + T1Y) * B0X + (-T0X + T0Y + T1X - T1Y) * B1X -
                                    (T0X + T0Y - T1X - T1Y) * (B0Y - B1Y)) *
                           ((T0X + T0Y - T1X - T1Y) * B0X + (-T0X - T0Y + T1X + T1Y) * B1X +
                            (T0X - T0Y - T1X + T1Y) * (B0Y - B1Y)) * pow(cos(phi), 0.4e1) +
                           0.4e1 * ((T0X - T1X) * B0X + (-T0X + T1X) * B1X - (T0Y - T1Y) * (B0Y - B1Y)) * sin(phi) *
                           ((T0Y - T1Y) * B0X + (-T0Y + T1Y) * B1X + (T0X - T1X) * (B0Y - B1Y)) * pow(cos(phi), 0.3e1) +
                           ((-0.4e1 * T0X * T0X + 0.8e1 * T0X * T1X - 0.4e1 * T1X * T1X +
                             0.2e1 * pow(T0Y - T1Y, 0.2e1)) *
                            B0X * B0X + ((0.8e1 * T0X * T0X - 0.16e2 * T0X * T1X + 0.8e1 * T1X * T1X -
                                          0.4e1 * pow(T0Y - T1Y, 0.2e1)) * B1X +
                                         0.8e1 * (T0Y - T1Y) * (T0X - T1X) * (B0Y - B1Y)) * B0X +
                            (-0.4e1 * T0X * T0X + 0.8e1 * T0X * T1X - 0.4e1 * T1X * T1X +
                             0.2e1 * pow(T0Y - T1Y, 0.2e1)) *
                            B1X * B1X - 0.8e1 * (T0Y - T1Y) * (T0X - T1X) * (B0Y - B1Y) * B1X +
                            0.2e1 * pow(T0X - T1X, 0.2e1) * pow(B0Y - B1Y, 0.2e1)) * pow(cos(phi), 0.2e1) -
                           0.4e1 * (T0X - T1X) * (B0X - B1X) * sin(phi) *
                           ((T0Y - T1Y) * B0X + (-T0Y + T1Y) * B1X + (T0X - T1X) * (B0Y - B1Y)) * cos(phi) +
                           0.2e1 * pow(T0X - T1X, 0.2e1) * pow(B0X - B1X, 0.2e1)) *
               (((T0X + T0Y - T1X - T1Y) * (T0X - T0Y - T1X + T1Y) * pow(cos(phi), 0.2e1) +
                 0.2e1 * sin(phi) * (T0Y - T1Y) * (T0X - T1X) * cos(phi) - pow(T0X - T1X, 0.2e1)) *
                sqrt(-pow(T0X * T1Y - T0Y * T1X + (PLX - PLY / tan(phi)) * T0Y - (PLX - PLY / tan(phi)) * T1Y, 0.2e1) /
                     (0.2e1 * sin(phi) * cos(phi) * T0X * T0Y - 0.2e1 * sin(phi) * cos(phi) * T0X * T1Y -
                      0.2e1 * sin(phi) * cos(phi) * T0Y * T1X + 0.2e1 * sin(phi) * cos(phi) * T1X * T1Y +
                      T0X * T0X * pow(cos(phi), 0.2e1) - 0.2e1 * T0X * T1X * pow(cos(phi), 0.2e1) -
                      pow(cos(phi), 0.2e1) * T0Y * T0Y + 0.2e1 * pow(cos(phi), 0.2e1) * T0Y * T1Y +
                      T1X * T1X * pow(cos(phi), 0.2e1) - pow(cos(phi), 0.2e1) * T1Y * T1Y - T0X * T0X +
                      0.2e1 * T0X * T1X - T1X * T1X)) * (PLX - PLY / tan(phi)) *
                (-0.24e2 * B0X * B0Y * B1X * B1Y * T0X * T1Y + 0.24e2 * B0X * B0Y * B1X * B1Y * T0Y * T1X -
                 0.3e1 * B0Y * B0Y * T1Y * pow(PLX - PLY / tan(phi), 0.3e1) +
                 0.3e1 * B1Y * B1Y * T0Y * pow(PLX - PLY / tan(phi), 0.3e1) +
                 0.3e1 * B0Y * B0Y * T0Y * pow(PLX - PLY / tan(phi), 0.3e1) +
                 0.6e1 * B0X * B0X * B1Y * B1Y * T0Y * (PLX - PLY / tan(phi)) -
                 0.6e1 * B0X * B0X * B1Y * B1Y * T1Y * (PLX - PLY / tan(phi)) -
                 0.8e1 * B0X * B1Y * B1Y * T0Y * pow(PLX - PLY / tan(phi), 0.2e1) +
                 0.8e1 * B0X * B1Y * B1Y * T1Y * pow(PLX - PLY / tan(phi), 0.2e1) +
                 0.6e1 * B0Y * B0Y * B1X * B1X * T0Y * (PLX - PLY / tan(phi)) -
                 0.6e1 * B0Y * B0Y * B1X * B1X * T1Y * (PLX - PLY / tan(phi)) -
                 0.8e1 * B0Y * B0Y * B1X * T0Y * pow(PLX - PLY / tan(phi), 0.2e1) +
                 0.8e1 * B0Y * B0Y * B1X * T1Y * pow(PLX - PLY / tan(phi), 0.2e1) +
                 0.4e1 * B0Y * B0Y * T0X * T1Y * pow(PLX - PLY / tan(phi), 0.2e1) -
                 0.4e1 * B0Y * B0Y * T0Y * T1X * pow(PLX - PLY / tan(phi), 0.2e1) -
                 0.6e1 * B0Y * B1Y * T0Y * pow(PLX - PLY / tan(phi), 0.3e1) +
                 0.6e1 * B0Y * B1Y * T1Y * pow(PLX - PLY / tan(phi), 0.3e1) +
                 0.4e1 * B1Y * B1Y * T0X * T1Y * pow(PLX - PLY / tan(phi), 0.2e1) -
                 0.4e1 * B1Y * B1Y * T0Y * T1X * pow(PLX - PLY / tan(phi), 0.2e1) -
                 0.3e1 * B1Y * B1Y * T1Y * pow(PLX - PLY / tan(phi), 0.3e1) -
                 0.12e2 * B0X * B0Y * B1X * B1Y * T0Y * (PLX - PLY / tan(phi)) +
                 0.12e2 * B0X * B0Y * B1X * B1Y * T1Y * (PLX - PLY / tan(phi)) +
                 0.12e2 * B0X * B0Y * B1Y * T0X * T1Y * (PLX - PLY / tan(phi)) -
                 0.12e2 * B0X * B0Y * B1Y * T0Y * T1X * (PLX - PLY / tan(phi)) +
                 0.12e2 * B0Y * B1X * B1Y * T0X * T1Y * (PLX - PLY / tan(phi)) -
                 0.12e2 * B0Y * B1X * B1Y * T0Y * T1X * (PLX - PLY / tan(phi)) +
                 0.8e1 * B0X * B0Y * B1Y * T0Y * pow(PLX - PLY / tan(phi), 0.2e1) -
                 0.8e1 * B0X * B0Y * B1Y * T1Y * pow(PLX - PLY / tan(phi), 0.2e1) -
                 0.12e2 * B0X * B1Y * B1Y * T0X * T1Y * (PLX - PLY / tan(phi)) +
                 0.12e2 * B0X * B1Y * B1Y * T0Y * T1X * (PLX - PLY / tan(phi)) -
                 0.12e2 * B0Y * B0Y * B1X * T0X * T1Y * (PLX - PLY / tan(phi)) +
                 0.12e2 * B0Y * B0Y * B1X * T0Y * T1X * (PLX - PLY / tan(phi)) +
                 0.8e1 * B0Y * B1X * B1Y * T0Y * pow(PLX - PLY / tan(phi), 0.2e1) -
                 0.8e1 * B0Y * B1X * B1Y * T1Y * pow(PLX - PLY / tan(phi), 0.2e1) -
                 0.8e1 * B0Y * B1Y * T0X * T1Y * pow(PLX - PLY / tan(phi), 0.2e1) +
                 0.8e1 * B0Y * B1Y * T0Y * T1X * pow(PLX - PLY / tan(phi), 0.2e1) +
                 0.12e2 * B0X * B0X * B1Y * B1Y * T0X * T1Y - 0.12e2 * B0X * B0X * B1Y * B1Y * T0Y * T1X +
                 0.12e2 * B0Y * B0Y * B1X * B1X * T0X * T1Y - 0.12e2 * B0Y * B0Y * B1X * B1X * T0Y * T1X) /
                (T0X * T1Y - T0Y * T1X + (PLX - PLY / tan(phi)) * T0Y - (PLX - PLY / tan(phi)) * T1Y) / 0.12e2 +
                ((B0X + B0Y - B1X - B1Y) * (B0X - B0Y - B1X + B1Y) * pow(cos(phi), 0.2e1) +
                 0.2e1 * sin(phi) * (B0Y - B1Y) * (B0X - B1X) * cos(phi) - pow(B0X - B1X, 0.2e1)) *
                sqrt(-pow(B0X * B1Y - B0Y * B1X + (PLX - PLY / tan(phi)) * B0Y - (PLX - PLY / tan(phi)) * B1Y, 0.2e1) /
                     (B0X * B0X * pow(cos(phi), 0.2e1) + 0.2e1 * B0X * B0Y * sin(phi) * cos(phi) -
                      0.2e1 * B0X * B1X * pow(cos(phi), 0.2e1) - 0.2e1 * B0X * B1Y * sin(phi) * cos(phi) -
                      B0Y * B0Y * pow(cos(phi), 0.2e1) - 0.2e1 * B0Y * B1X * sin(phi) * cos(phi) +
                      0.2e1 * B0Y * B1Y * pow(cos(phi), 0.2e1) + B1X * B1X * pow(cos(phi), 0.2e1) +
                      0.2e1 * B1X * B1Y * sin(phi) * cos(phi) - B1Y * B1Y * pow(cos(phi), 0.2e1) - B0X * B0X +
                      0.2e1 * B0X * B1X - B1X * B1X)) * (PLX - PLY / tan(phi)) *
                (-0.24e2 * B0X * B1Y * T0X * T0Y * T1X * T1Y + 0.24e2 * B0Y * B1X * T0X * T0Y * T1X * T1Y +
                 0.3e1 * B0Y * T0Y * T0Y * pow(PLX - PLY / tan(phi), 0.3e1) -
                 0.3e1 * B1Y * T0Y * T0Y * pow(PLX - PLY / tan(phi), 0.3e1) +
                 0.3e1 * B0Y * T1Y * T1Y * pow(PLX - PLY / tan(phi), 0.3e1) +
                 0.4e1 * B0X * B1Y * T0Y * T0Y * pow(PLX - PLY / tan(phi), 0.2e1) +
                 0.4e1 * B0X * B1Y * T1Y * T1Y * pow(PLX - PLY / tan(phi), 0.2e1) -
                 0.4e1 * B0Y * B1X * T0Y * T0Y * pow(PLX - PLY / tan(phi), 0.2e1) -
                 0.4e1 * B0Y * B1X * T1Y * T1Y * pow(PLX - PLY / tan(phi), 0.2e1) +
                 0.6e1 * B0Y * T0X * T0X * T1Y * T1Y * (PLX - PLY / tan(phi)) -
                 0.8e1 * B0Y * T0X * T1Y * T1Y * pow(PLX - PLY / tan(phi), 0.2e1) +
                 0.6e1 * B0Y * T0Y * T0Y * T1X * T1X * (PLX - PLY / tan(phi)) -
                 0.8e1 * B0Y * T0Y * T0Y * T1X * pow(PLX - PLY / tan(phi), 0.2e1) -
                 0.6e1 * B0Y * T0Y * T1Y * pow(PLX - PLY / tan(phi), 0.3e1) -
                 0.6e1 * B1Y * T0X * T0X * T1Y * T1Y * (PLX - PLY / tan(phi)) +
                 0.8e1 * B1Y * T0X * T1Y * T1Y * pow(PLX - PLY / tan(phi), 0.2e1) -
                 0.6e1 * B1Y * T0Y * T0Y * T1X * T1X * (PLX - PLY / tan(phi)) +
                 0.8e1 * B1Y * T0Y * T0Y * T1X * pow(PLX - PLY / tan(phi), 0.2e1) +
                 0.6e1 * B1Y * T0Y * T1Y * pow(PLX - PLY / tan(phi), 0.3e1) -
                 0.3e1 * B1Y * T1Y * T1Y * pow(PLX - PLY / tan(phi), 0.3e1) +
                 0.12e2 * B0X * B1Y * T0X * T0X * T1Y * T1Y +
                 0.12e2 * B0X * B1Y * T0Y * T0Y * T1X * T1X - 0.12e2 * B0Y * B1X * T0X * T0X * T1Y * T1Y -
                 0.12e2 * B0Y * B1X * T0Y * T0Y * T1X * T1X +
                 0.12e2 * B0X * B1Y * T0X * T0Y * T1Y * (PLX - PLY / tan(phi)) +
                 0.12e2 * B0X * B1Y * T0Y * T1X * T1Y * (PLX - PLY / tan(phi)) -
                 0.12e2 * B0Y * B1X * T0X * T0Y * T1Y * (PLX - PLY / tan(phi)) -
                 0.12e2 * B0Y * B1X * T0Y * T1X * T1Y * (PLX - PLY / tan(phi)) -
                 0.12e2 * B0Y * T0X * T0Y * T1X * T1Y * (PLX - PLY / tan(phi)) +
                 0.12e2 * B1Y * T0X * T0Y * T1X * T1Y * (PLX - PLY / tan(phi)) +
                 0.12e2 * B0Y * B1X * T0X * T1Y * T1Y * (PLX - PLY / tan(phi)) +
                 0.12e2 * B0Y * B1X * T0Y * T0Y * T1X * (PLX - PLY / tan(phi)) +
                 0.8e1 * B0Y * B1X * T0Y * T1Y * pow(PLX - PLY / tan(phi), 0.2e1) +
                 0.8e1 * B0Y * T0X * T0Y * T1Y * pow(PLX - PLY / tan(phi), 0.2e1) +
                 0.8e1 * B0Y * T0Y * T1X * T1Y * pow(PLX - PLY / tan(phi), 0.2e1) -
                 0.8e1 * B1Y * T0X * T0Y * T1Y * pow(PLX - PLY / tan(phi), 0.2e1) -
                 0.8e1 * B1Y * T0Y * T1X * T1Y * pow(PLX - PLY / tan(phi), 0.2e1) -
                 0.12e2 * B0X * B1Y * T0X * T1Y * T1Y * (PLX - PLY / tan(phi)) -
                 0.12e2 * B0X * B1Y * T0Y * T0Y * T1X * (PLX - PLY / tan(phi)) -
                 0.8e1 * B0X * B1Y * T0Y * T1Y * pow(PLX - PLY / tan(phi), 0.2e1)) /
                (B0X * B1Y - B0Y * B1X + (PLX - PLY / tan(phi)) * B0Y - (PLX - PLY / tan(phi)) * B1Y) / 0.12e2);
    }
};

#endif // VISIBILITY_2_EXAMPLES_NUMERICALINTEGRALOpt_H
