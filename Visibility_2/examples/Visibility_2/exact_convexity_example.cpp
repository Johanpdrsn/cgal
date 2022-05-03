//
// Created by Johan Pedersen on 25/04/2022.
//

#include <iostream>
#include <cmath>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Visibility_2/RedBlackTree.h>
#include <CGAL/Visibility_2/BinaryTree.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Rotational_sweep_visibility_2.h>
#include <CGAL/Cartesian.h>


//typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Cartesian<double> Kernel;
typedef CGAL::Polygon_2<Kernel> Polygon_2;
typedef Kernel::Point_2 Point_2;

typedef CGAL::Polygon_2<Kernel> Polygon_2;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Triangle_2 Triangle_2;
typedef Kernel::Segment_2 Segment_2;

using namespace std;

struct FaceInfo2 {
    FaceInfo2() = default;

    int nesting_level;

    bool in_domain() const {
        return nesting_level % 2 == 1;
    }
};

typedef CGAL::Triangulation_vertex_base_2<Kernel> Vb;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, Kernel> Fbb;
typedef CGAL::Constrained_triangulation_face_base_2<Kernel, Fbb> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> TDS;
typedef CGAL::Exact_intersections_tag Itag;
typedef CGAL::Constrained_triangulation_2<Kernel, TDS, Itag> CT;
typedef CGAL::Constrained_triangulation_plus_2<CT> CTP;

class Convexity_measure_exact_2 final {
public:


    Polygon_2 polygon;
    vector<CTP::Face> triangles;
    vector<Segment_2> diagonals;

    CTP triangulation;
    RedBlackTree<Polygon_2> bst;
    BinaryTree<CTP> tree;


    explicit Convexity_measure_exact_2(Polygon_2 &input_polygon) : polygon(input_polygon) {
        triangulate();
        sweep_diagonals();

        tree = BinaryTree<CTP>(triangles);
        generate_tree_rec(tree.Root());
//        generate_diagonals();
    }


private:

    static void mark_domains(CTP &ctp,
                             CTP::Face_handle start,
                             int index,
                             std::list<CTP::Edge> &border) {
        if (start->info().nesting_level != -1) {
            return;
        }
        std::list<CTP::Face_handle> queue;
        queue.push_back(start);
        while (!queue.empty()) {
            CTP::Face_handle fh = queue.front();
            queue.pop_front();
            if (fh->info().nesting_level == -1) {
                fh->info().nesting_level = index;
                for (int i = 0; i < 3; ++i) {
                    CTP::Edge e(fh, i);
                    CTP::Face_handle n = fh->neighbor(i);
                    if (n->info().nesting_level == -1) {
                        if (ctp.is_constrained(e)) border.push_back(e);
                        else queue.push_back(n);
                    }
                }
            }
        }
    }

    static void mark_domains(CTP &ctp) {

        for (CTP::Face_handle f: ctp.all_face_handles()) {
            f->info().nesting_level = -1;
        }
        std::list<CTP::Edge> border;
        mark_domains(ctp, ctp.infinite_face(), 0, border);
        while (!border.empty()) {
            const CTP::Edge e = border.front();
            border.pop_front();
            CTP::Face_handle n = e.first->neighbor(e.second);
            if (n->info().nesting_level == -1) {
                mark_domains(ctp, n, e.first->info().nesting_level + 1, border);
            }
        }
    }

    void triangulate() {

        triangulation.insert_constraint(polygon.vertices_begin(), polygon.vertices_end(), true);
        mark_domains(triangulation);

        for (const auto &face: triangulation.finite_face_handles()) {
            if (!face->info().in_domain()) {
                triangulation.delete_face(face);
            } else {
                triangles.emplace_back(*face);
            }
        }
    }


    void generate_tree_rec(BinaryTree<CTP>::Node *node) {
        auto leftNode = tree.EmptyNode();
        auto rightNode = tree.EmptyNode();

        if (node->data.size() <= 1) {
            return;
        } else {
            auto size = static_cast<double>(node->data.size());
            long bound = std::floor(size / 2.0);

            for (const auto &it: node->data) {
                if (leftNode->data.size() < bound) {
                    leftNode->data.emplace_back(it);
                } else {
                    rightNode->data.emplace_back(it);
                }
            }
        }

        node->left = leftNode;
        generate_tree_rec(leftNode);
        node->right = rightNode;
        generate_tree_rec(rightNode);
    }

    void sweep_diagonals() {
        typedef CGAL::Arr_segment_traits_2<Kernel> Traits_2;
        typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;
        typedef CGAL::Rotational_sweep_visibility_2<Arrangement_2> RSV;


        // insert geometry into the arrangement
        Arrangement_2 env;
        CGAL::insert_non_intersecting_curves(env, polygon.edges_begin(), polygon.edges_end());


        for (auto he: env.halfedge_handles()) {
            if (!he->face()->is_unbounded()) {
                Arrangement_2 output_arr;
                RSV rsv(env);

                auto p = he->target()->point();
                rsv.compute_visibility(p, he, output_arr);

                for (auto v: output_arr.vertex_handles()) {
                    bool vis_point_in_poly = find(polygon.vertices_begin(), polygon.vertices_end(), v->point()) !=
                                             polygon.vertices_end();
                    if (vis_point_in_poly && p != v->point()) {
                        auto seg = Segment_2(p, v->point());
                        bool c = find(polygon.edges_begin(), polygon.edges_end(), seg) != polygon.edges_end();
                        bool d =
                                find(polygon.edges_begin(), polygon.edges_end(), seg.opposite()) != polygon.edges_end();
                        bool e = find(diagonals.begin(), diagonals.end(), seg) != diagonals.end();
                        bool f = find(diagonals.begin(), diagonals.end(), seg.opposite()) != diagonals.end();

                        if (!(c || d || e || f)) {
                            diagonals.emplace_back(seg);
                        }
                    }
                }
            }
        }

        std::sort(diagonals.begin(), diagonals.end(),
                  [&](Segment_2 &A, Segment_2 &B) {
                      return segment_angle(A) < segment_angle(B);
                  });
    }

    static double segment_angle(const Segment_2 &seg) {
        if (seg.is_horizontal()) {
            return 0;
        } else if (seg.is_vertical()) {
            return M_PI_2;
        } else {
            auto slope = CGAL::to_double(seg.direction().dy() / seg.direction().dx());
            return atan(slope) > 0 ? atan(slope) : atan(slope) + M_PI;
        }
    }
};


int main() {
    CGAL::IO::set_pretty_mode(std::cout);

    Polygon_2 polygon;
    polygon.push_back(Point_2{0.0, 0.0});
    polygon.push_back(Point_2{-2.0, -1.0});
    polygon.push_back(Point_2{0.0, -2.0});
    polygon.push_back(Point_2{1.0, -4.0});
    polygon.push_back(Point_2{2.0, -2.0});
    polygon.push_back(Point_2{4.0, -1.0});
    polygon.push_back(Point_2{2.0, 0.0});
    polygon.push_back(Point_2{1.0, 2.0});


    auto test = Convexity_measure_exact_2(polygon);

    test.tree.prettyPrint();

    return 0;
}