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

        T t1 = (int) ((double) V0Y * (double) V0X);
        T t2 = (int) ((double) V1Y * (double) t1);
        T t3 = (int) ((double) T0Y - (double) T1Y);
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
        T t25 = V1X * V0Y;
        T t33 = -t16 * t3 * V1Y * t1 + t17 * V1Y * t25 + t16 * t3 * t23 + t17 * t2 - t17 * t20 - t17 * t23;
        T t34 = (int) ((double) B0Y - (double) B1Y);
        T t35 = (int) ((double) t34 * (double) t34);
        T t40 = t8 * t34 - t13 * (-B0X + B1X);
        T t41 = t40 * t40;
        T t42 = 1 / t41;
        T t44 = t12 * t1;
        T t45 = rho * t20;
        T t46 = rho * t23;
        T t47 = t12 * t25;
        T t48 = t44 - t45 - t46 + t47 - t2 + t23;
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
        T t69 = (int) ((double) PLY * (double) PLY);
        T t70 = (int) ((double) t69 * (double) t69);
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
        T t179 = (double) ((-t67 * t72 * t70 + t67 * t64) * (-t42 * t35 * t33 + t54 * t53 * t50)) / 0.8e1 +
                 (double) ((t124 * t71 * t8 * t69 * PLY + t124 * t63 * t62) *
                           (-t107 * t42 * (-2 * t102 * t99 * t33 + t65 * t35 * t94) -
                            t54 * (2 * t53 * t82 * t3 * t13 * t109 + t13 * t116 * t82 * t50))) / 0.6e1 +
                 (double) ((-t107 * t71 * t69 + t107 * t63) *
                           (-t107 * t42 * (-2 * t102 * t99 * t94 + t65 * t136 * t33) -
                            t54 * (-2 * t116 * t107 * t3 * t65 * t109 - t34 * t107 * t52 * t65 * t143))) / 0.4e1 -
                 (double) (t168 * t107 * t42 * t65 * t136 * t94) / 0.2e1 -
                 (double) (t168 * t54 * t99 * t124 * t52 * t123 * t143) / 0.2e1;

        return t179;
    };

    static T ComputeZero(T T0X, T T0Y, T T1X, T T1Y, T B0X, T B0Y, T B1X, T B1Y,
                         T PLX, T PLY, T PRX, T PRY, T V0X, T V0Y, T V1X, T V1Y, T rho) {

        T t1 = V0Y * V0X;
        T t2 = T0Y - T1Y;
        T t7 = V0X * t2 + (-T0X + T1X) * V0Y;
        T t8 = 0.1e1 / t7;
        T t11 = V0Y * V0Y;
        T t12 = V1X * t11;
        T t15 = -t8 * t2 * V1Y * t1 + t8 * t2 * t12;
        T t16 = (int) (B0Y - B1Y);
        T t17 = (double) t16 * (double) t16;
        T t22 = V0X * (double) t16 + (-B0X + B1X) * V0Y;
        T t23 = t22 * t22;
        T t24 = 0.1e1 / t23;
        T t27 = -V1Y * t1 + t12;
        T t28 = t2 * t2;
        T t29 = t28 * t27;
        T t30 = t7 * t7;
        T t31 = 0.1e1 / t30;
        T t32 = (double) t16 * t31;
        T t33 = 0.1e1 / t22;
        T t41 = PRX * (V0Y + PRY) - PRY * (V0X + PRX);
        T t42 = t41 * t41;
        T t43 = t42 * t42;
        T t44 = t11 * t11;
        T t45 = 0.1e1 / t44;
        T t47 = PLY * PLY;
        T t48 = t47 * t47;
        T t49 = V0X * V0X;
        T t50 = t49 * t49;
        T t58 = T0X * T1Y - T0Y * T1X;
        T t64 = -t8 * t58 * V1Y * t1 + t8 * t58 * t12;
        T t69 = B0X * B1Y - B0Y * B1X;
        T t71 = (double) t16 * t11;
        T t76 = 0.1e1 / t11;
        T t78 = (int) (t58 * t27);
        T t79 = (int) (t31 * t2);
        T t90 = 0.1e1 / t11 / V0Y;
        T t102 = t69 * t69;
        T t108 = t58 * t58;
        T t109 = t108 * t27;
        T t124 = 0.1e1 / V0Y;
        T t128 = t124 * V0X * PLY + t124 * t41;
        return (-t45 * t50 * t48 + t45 * t43) * (-t17 * t15 * t24 + t33 * t32 * t29) / 0.8e1 +
               (t90 * t49 * V0X * t47 * PLY + t90 * t42 * t41) *
               (-t76 * t24 * (t11 * t17 * t64 + 0.2e1 * t71 * t69 * t15) -
                t33 * (-(double) (2 * t16 * t79 * t78) - t69 * t31 * t29)) / 0.6e1 + (-t76 * t49 * t47 + t76 * t42) *
                                                                                     (-t76 * t24 * (t11 * t102 * t15 +
                                                                                                    0.2e1 * t71 * t69 *
                                                                                                    t64) - t33 *
                                                                                                           (-0.2e1 *
                                                                                                            t69 *
                                                                                                            (double) t79 *
                                                                                                            (double) t78 -
                                                                                                            t32 *
                                                                                                            t109)) /
                                                                                     0.4e1 -
               t128 * t24 * t102 * t64 / 0.2e1 + t128 * t33 * t69 * t31 * t109 / 0.2e1;

    };

    static T ComputeOne(T T0X, T T0Y, T T1X, T T1Y, T B0X, T B0Y, T B1X, T B1Y,
                        T PLX, T PLY, T PRX, T PRY, T V0X, T V0Y, T V1X, T V1Y, T rho) {

        T t1 = V1Y * V1Y;
        T t2 = t1 * V0X;
        T t3 = T0Y - T1Y;
        T t7 = V1X * t3 + (-T0X + T1X) * V1Y;
        T t8 = 0.1e1 / t7;
        T t11 = V0Y * V1X;
        T t15 = t8 * t3 * V1Y * t11 - t8 * t3 * t2;
        T t16 = (int) (B0Y - B1Y);
        T t17 = (double) t16 * (double) t16;
        T t22 = V1X * (double) t16 + (-B0X + B1X) * V1Y;
        T t23 = t22 * t22;
        T t24 = 0.1e1 / t23;
        T t27 = V1Y * t11 - t2;
        T t28 = t3 * t3;
        T t29 = t28 * t27;
        T t30 = t7 * t7;
        T t31 = 0.1e1 / t30;
        T t32 = (double) t16 * t31;
        T t33 = 0.1e1 / t22;
        T t41 = PRX * (V1Y + PRY) - PRY * (V1X + PRX);
        T t42 = t41 * t41;
        T t43 = t42 * t42;
        T t44 = t1 * t1;
        T t45 = 0.1e1 / t44;
        T t47 = PLY * PLY;
        T t48 = t47 * t47;
        T t49 = V1X * V1X;
        T t50 = t49 * t49;
        T t58 = T0X * T1Y - T0Y * T1X;
        T t64 = t8 * t58 * V1Y * t11 - t8 * t58 * t2;
        T t69 = B0X * B1Y - B0Y * B1X;
        T t71 = (double) t16 * t1;
        T t76 = 0.1e1 / t1;
        T t78 = (int) (t58 * t27);
        T t79 = (int) (t31 * t3);
        T t90 = 0.1e1 / t1 / V1Y;
        T t102 = t69 * t69;
        T t108 = t58 * t58;
        T t109 = t108 * t27;
        T t124 = 0.1e1 / V1Y;
        T t128 = t124 * V1X * PLY + t124 * t41;
        return (-t45 * t50 * t48 + t45 * t43) * (-t17 * t15 * t24 + t33 * t32 * t29) / 0.8e1 +
               (t90 * t49 * V1X * t47 * PLY + t90 * t42 * t41) *
               (-t76 * t24 * (t1 * t17 * t64 + 0.2e1 * t71 * t69 * t15) -
                t33 * (-(double) (2 * t16 * t79 * t78) - t69 * t31 * t29)) / 0.6e1 + (-t76 * t49 * t47 + t76 * t42) *
                                                                                     (-t76 * t24 * (t1 * t102 * t15 +
                                                                                                    0.2e1 * t71 * t69 *
                                                                                                    t64) - t33 *
                                                                                                           (-0.2e1 *
                                                                                                            t69 *
                                                                                                            (double) t79 *
                                                                                                            (double) t78 -
                                                                                                            t32 *
                                                                                                            t109)) /
                                                                                     0.4e1 -
               t128 * t24 * t102 * t64 / 0.2e1 + t128 * t33 * t69 * t31 * t109 / 0.2e1;

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
                (double) (V0X * V0Y * V1Y * rho - V0X * V1Y * V1Y * rho - V0Y * V0Y * V1X * rho +
                          V0Y * V1X * V1Y * rho - V0X * V0Y * V1Y + V0Y * V0Y * V1X) *
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
                                                                                                                     V0X *
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
                                                                                                                     V0X *
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
                                                                                                                      V0X *
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
                                                                                                                       V0X *
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
               (double) (V0X * V0Y * V1Y * rho - V0X * V1Y * V1Y * rho - V0Y * V0Y * V1X * rho + V0Y * V1X * V1Y * rho -
                         V0X * V0Y * V1Y + V0Y * V0Y * V1X) *
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


};

#endif // VISIBILITY_2_EXAMPLES_NUMERICALINTEGRALOpt_H
