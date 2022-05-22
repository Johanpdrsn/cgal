//
// Created by Johan Pedersen on 18/05/2022.
//

#ifndef VISIBILITY_2_EXAMPLES_EXACT_CONVEXITY_MEASURE_2_H
#define VISIBILITY_2_EXAMPLES_EXACT_CONVEXITY_MEASURE_2_H

#include <CGAL/Polygon_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Rotational_sweep_visibility_2.h>
#include <CGAL/aff_transformation_tags.h>

#include <CGAL/Visibility_2/BinaryTree.h>
#include <CGAL/Visibility_2/NumericalIntegral.h>
#include <CGAL/Visibility_2/CentroidDecomposition.h>

using namespace std;

class Convexity_measure_exact_2 final {
public:
    typedef CGAL::Cartesian<double> K;

    typedef typename K::Point_2 Point_2;
    typedef typename K::Triangle_2 Triangle_2;
    typedef typename K::Segment_2 Segment_2;
    typedef typename K::Direction_2 Direction_2;
    typedef typename K::FT FT;

    typedef CGAL::Polygon_2<K> Polygon_2;
    typedef CGAL::Arr_segment_traits_2<K> Traits_2;
    typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;
    typedef CGAL::Rotational_sweep_visibility_2<Arrangement_2> RSV;


    explicit Convexity_measure_exact_2(Polygon_2 &input_polygon) : polygon(input_polygon), tree() {
        triangulate();
        sweep_diagonals();
    }

    FT generate_tree() {
        tree.SetRoot(triangles);
        return decompose_tree_rec(tree.Root(), startingDiagonals);
    }

private:
    struct FaceInfo2 {
        FaceInfo2() = default;

        int nesting_level;
        int id;

        bool in_domain() const {
            return nesting_level % 2 == 1;
        }
    };

    typedef typename CGAL::Triangulation_vertex_base_2<K> Vb;
    typedef typename CGAL::Triangulation_face_base_with_info_2<FaceInfo2, K> Fbb;
    typedef typename CGAL::Constrained_triangulation_face_base_2<K, Fbb> Fb;
    typedef CGAL::Triangulation_data_structure_2<Vb, Fb> TDS;
    typedef CGAL::Exact_intersections_tag Itag;
    typedef typename CGAL::Constrained_triangulation_2<K, TDS, Itag> CT;
    typedef CGAL::Constrained_triangulation_plus_2<CT> CTP;
    typedef typename CTP::Face_handle Face_handle;


    struct SegmentHash {
        size_t operator()(const Segment_2 &t) const {
            auto h = std::hash<double>{};
            return h(CGAL::to_double((t.source().x()))) ^ h(CGAL::to_double((t.source().y())))
                   ^ h(CGAL::to_double((t.target().x()))) ^ h(CGAL::to_double((t.target().y())));
        }
    };

    struct PointHashFunction {
        size_t operator()(const typename CTP::Point &t) const {
            return std::hash<double>()(CGAL::to_double((t.x()))) ^
                   std::hash<double>()(CGAL::to_double((t.y())));
        }
    };

    struct Corner {
        Point_2 target;
        Segment_2 top;
        Segment_2 bot;
        Direction_2 direction{};

        bool operator==(const Corner &A) const {
            return this->target == A.target;
        }
    };

    using DiagonalMap = std::unordered_map<Segment_2, Segment_2, SegmentHash>;


    Polygon_2 polygon;
    std::vector<Face_handle> triangles;

    CTP triangulation;
    BinaryTree<CTP> tree;
    DiagonalMap startingDiagonals;

    // check if 2 faces are equal
    static bool face_equality(const Face_handle &A, const Face_handle &B) {
        std::unordered_set<typename CTP::Point, PointHashFunction> s{A->vertex(0)->point(), A->vertex(1)->point(),
                                                                     A->vertex(2)->point()};
        for (int i = 0; i < 3; i++) {
            if (s.count(B->vertex(i)->point()) == 0)
                return false;
        }
        return true;
    }

    // check if face is contained in list of faces
    static bool
    face_in_list(const std::vector<Face_handle> &list, const Face_handle &face) {
        return std::find_if(list.begin(), list.end(),
                            [&](const Face_handle &A) { return face_equality(A, face); }) != list.end();
    }

    // count vertices in list of faces
    static double count_vertices(const std::vector<Face_handle> &list) {
        std::unordered_set<typename CTP::Point, PointHashFunction> s;
        for (const auto &face: list) {
            for (int i = 0; i < 3; i++) {
                s.insert(face->vertex(i)->point());
            }
        }
        return static_cast<double>(s.size());
    }

    // check if edge is in list of faces
    static bool is_edge_in_faces(const std::vector<Face_handle> &faces, const Segment_2 &seg) {
        for (const auto &face: faces) {
            for (int i = 0; i < 3; i++) {
                if (seg == Segment_2(face->vertex(i)->point(), face->vertex((i + 1) % 3)->point()))
                    return true;
            }
        }
        return false;
    }

    static void mark_domains(CTP &ctp,
                             Face_handle start,
                             int index,
                             std::list<typename CTP::Edge> &border) {
        if (start->info().nesting_level != -1) {
            return;
        }
        std::list<Face_handle> queue;
        queue.push_back(start);
        while (!queue.empty()) {
            Face_handle fh = queue.front();
            queue.pop_front();
            if (fh->info().nesting_level == -1) {
                fh->info().nesting_level = index;
                for (int i = 0; i < 3; ++i) {
                    typename CTP::Edge e(fh, i);
                    Face_handle n = fh->neighbor(i);
                    if (n->info().nesting_level == -1) {
                        if (ctp.is_constrained(e)) border.push_back(e);
                        else queue.push_back(n);
                    }
                }
            }
        }
    }

    static void mark_domains(CTP &ctp) {
        for (Face_handle f: ctp.all_face_handles()) {
            f->info().nesting_level = -1;
        }
        std::list<typename CTP::Edge> border;
        mark_domains(ctp, ctp.infinite_face(), 0, border);
        while (!border.empty()) {
            const typename CTP::Edge e = border.front();
            border.pop_front();
            Face_handle n = e.first->neighbor(e.second);
            if (n->info().nesting_level == -1) {
                mark_domains(ctp, n, e.first->info().nesting_level + 1, border);
            }
        }
    }

    void triangulate() {
        triangulation.insert_constraint(polygon.vertices_begin(), polygon.vertices_end(), true);
        mark_domains(triangulation);

        int id = 0;
        for (const auto &face: triangulation.finite_face_handles()) {
            if (!face->info().in_domain()) {
                triangulation.delete_face(face);
            } else {
                face->info().id = id++;
                triangles.emplace_back(face);
            }
        }

        std::cout << "TRIANGULATION" << std::endl;
        for (auto a: triangles) {
            std::cout << triangulation.triangle(a) << a->info().id << std::endl;
        }
    }

    // find the closest intersecting segment in p1/p2
    Segment_2 find_nearest_intersecting_segment(const Segment_2 &s) const {
        typename K::Ray_2 r{s.target(), s.direction()};
        FT min_y = INT_MAX;
        Segment_2 left_plus;

        for (auto it = polygon.edges_begin(); it != polygon.edges_end(); it++) {
            const auto result = intersection(r, *it);
            if (result) {
                if (const Point_2 *p = boost::get<Point_2>(&*result)) {
                    FT diff = abs(p->y() - s.target().y());
                    if (diff > 0 && diff < min_y) {
                        min_y = diff;
                        left_plus = *it;
                    }
                }
            }
        }
        return left_plus;
    }

    // find eta edge for startingDiagonal
    Segment_2 find_eta(const Segment_2 &s) const {
        bool is_edge = std::any_of(polygon.edges_begin(), polygon.edges_end(),
                                   [&s](const Segment_2 &a) { return a == s; });
        bool is_opposite_edge = std::any_of(polygon.edges_begin(), polygon.edges_end(),
                                            [&s](const Segment_2 &a) { return a == s.opposite(); });

        auto cur = find(polygon.vertices_circulator(), polygon.vertices_circulator() - 1, s.target());
        auto prev = cur - 1;
        auto next = cur + 1;

        // clockwise edge
        if (is_opposite_edge) {
            if (right_turn(*next, *cur, *prev)) {
                return s;
            } else {
                return {s.target(), *next};
            }
            // counterclockwise edge
        } else if (is_edge) {
            if (!right_turn(*prev, *cur, *next)) {
                return {s.target(), *next};
            } else {
                return find_nearest_intersecting_segment(s);
            }
            // non edge
        } else {
            if (!right_turn(s.source(), s.target(), *next)) {
                return {s.target(), *next};
            } else {
                return find_nearest_intersecting_segment(s);
            }
        }
    }

    // find all diagonals in polygon
    void sweep_diagonals() {
        // insert geometry into the arrangement
        Arrangement_2 env;
        CGAL::insert_non_intersecting_curves(env, polygon.edges_begin(), polygon.edges_end());

        for (const auto &he: env.halfedge_handles()) {

            if (!he->face()->is_unbounded()) {
                Arrangement_2 output_arr;
                RSV rsv(env);

                auto p = he->target()->point();
                rsv.compute_visibility(p, he, output_arr);

                for (const auto &v: output_arr.vertex_handles()) {
                    bool vis_point_in_poly = find(polygon.vertices_begin(), polygon.vertices_end(), v->point()) !=
                                             polygon.vertices_end();

                    // CHeck if point in the visible polygon exists in the input polygon
                    if (vis_point_in_poly && p != v->point()) {
                        auto seg = Segment_2(p, v->point());
                        startingDiagonals.insert(std::make_tuple(seg, find_eta(seg)));
                    }
                }
            }
        }
    }

    // find the startingDiagonal splitting two nodes
    Segment_2 find_diagonal(const Face_handle &leftFace, typename BinaryTree<CTP>::Node *rightNode) const {
        int res;
        for (const auto &rightFace: rightNode->data) {
            if (leftFace->has_neighbor(rightFace, res)) {
                return Segment_2(triangulation.segment(typename CTP::Edge{leftFace, res}));
            }
        }

        throw std::logic_error{"No shared startingDiagonal between polygons"};
    }

    // find startingDiagonal subset that cross the splitting startingDiagonal
    std::tuple<DiagonalMap, DiagonalMap, DiagonalMap>
    find_diagonal_subsets(typename BinaryTree<CTP>::Node *leftNode, typename BinaryTree<CTP>::Node *rightNode,
                          DiagonalMap diagonals, const Segment_2 &diagonal) {
        DiagonalMap leftDiagonals;
        DiagonalMap rightDiagonals;
        DiagonalMap crossingDiagonals;

        std::unordered_set<Point_2, PointHashFunction> leftPoints;
        std::unordered_set<Point_2, PointHashFunction> rightPoints;

        for (const auto &f: leftNode->data) {
            for (int i = 0; i < 3; i++) {
                leftPoints.insert(f->vertex(i)->point());
            }
        }
        for (const auto &f: rightNode->data) {
            for (int i = 0; i < 3; i++) {
                rightPoints.insert(f->vertex(i)->point());
            }
        }

        for (const auto &diag: diagonals) {
            if ((leftPoints.count(diag.first.source()) == 1 && rightPoints.count(diag.first.target()) == 1) ||
                (rightPoints.count(diag.first.source()) == 1 && leftPoints.count(diag.first.target())) == 1) {
                crossingDiagonals[diag.first] = diag.second;
            } else if (leftPoints.count(diag.first.source()) == 1 && leftPoints.count(diag.first.target()) == 1) {
                leftDiagonals[diag.first] = diag.second;
            } else if (rightPoints.count(diag.first.source()) == 1 && rightPoints.count(diag.first.target()) == 1) {
                rightDiagonals[diag.first] = diag.second;
            }
        }

        // update eta edges for left subset
        for (auto &s: leftDiagonals) {
            if (!is_edge_in_faces(leftNode->data, diagonals[s.first]))
                leftDiagonals[s.first] = diagonal;
        }

        // update eta edges for right subset
        for (const auto &s: rightDiagonals) {
            if (!is_edge_in_faces(rightNode->data, diagonals[s.first]))
                rightDiagonals[s.first] = diagonal;
        }

        return make_tuple(leftDiagonals, crossingDiagonals, rightDiagonals);
    }

    std::vector<Segment_2> sort_diagonals(const DiagonalMap &diagMap, const Segment_2 &splittingDiagonal) {
        std::vector<Segment_2> crossingDiagonals;

        for (const auto &v: diagMap) {
            crossingDiagonals.emplace_back(v.first);
        }

        std::sort(crossingDiagonals.begin(), crossingDiagonals.end(),
                  [&splittingDiagonal](const Segment_2 &A, const Segment_2 &B) {
                      if (A == splittingDiagonal || A == splittingDiagonal.opposite())
                          return true;
                      else if (B == splittingDiagonal || B == splittingDiagonal.opposite())
                          return false;
                      else {
                          auto d = std::max(splittingDiagonal.direction(), splittingDiagonal.opposite().direction());

                          auto AMin = std::min(A.direction(), A.opposite().direction());
                          auto AMax = std::max(A.direction(), A.opposite().direction());

                          auto BMin = std::min(B.direction(), B.opposite().direction());
                          auto BMax = std::max(B.direction(), B.opposite().direction());

                          return AMin.counterclockwise_in_between(d, BMin) || AMax.counterclockwise_in_between(d, BMax);
                      }
                  });

        return crossingDiagonals;
    }

    // handles the event list
    FT create_event_list(typename BinaryTree<CTP>::Node *leftNode, typename BinaryTree<CTP>::Node *rightNode,
                         DiagonalMap &diagMap, const Segment_2 &splittingDiagonal) {

        auto crossingDiagonals = sort_diagonals(diagMap, splittingDiagonal);
        cout << "SPLIT" << splittingDiagonal << endl;

        FT final_res = 0;
        auto eventList = std::list<Corner>();

        auto f = crossingDiagonals.front();
        auto b = crossingDiagonals[1];
        crossingDiagonals.erase(crossingDiagonals.begin(), next(crossingDiagonals.begin(), 2));

        if (f.target().x() > b.target().x()) {
            std::swap(f, b);
        }

        eventList.push_front({f.target(), diagMap[b], diagMap[f], b.direction()});
        eventList.push_back({b.target(), diagMap[b], diagMap[f], b.direction()});


        auto lastPlaced = eventList.begin();

        for (auto crossingDiagonal = crossingDiagonals.begin();
             crossingDiagonal != crossingDiagonals.end(); crossingDiagonal++) {

            if (crossingDiagonal->direction() > crossingDiagonal->opposite().direction()) {
                continue;
            }

            Corner newCorner = {crossingDiagonal->target(), diagMap[*crossingDiagonal],
                                diagMap[crossingDiagonal->opposite()].opposite(),
                                crossingDiagonal->direction()};
            auto target_in_list = find(eventList.begin(), eventList.end(), newCorner);
            auto source_in_list = find(eventList.begin(), eventList.end(), Corner{crossingDiagonal->source()});


            for (auto a: eventList) {
                cout << a.target << endl;
            }
            cout << "Last placed: " << lastPlaced->target << endl;
            cout << crossingDiagonal->direction() << endl;


            if (target_in_list == eventList.end()) {
                final_res += calc(*(lastPlaced), *next(lastPlaced), newCorner.direction, splittingDiagonal);

                lastPlaced = (eventList.insert(next(lastPlaced), newCorner));

                prev(lastPlaced)->direction = newCorner.direction;
                lastPlaced->direction = newCorner.direction;

                if (source_in_list != eventList.end()) {
                    source_in_list->top = newCorner.top;
                    source_in_list->bot = newCorner.bot;
                }
                prev(lastPlaced)->top = lastPlaced->top;

            } else if (source_in_list == eventList.end()) {
                newCorner.target = crossingDiagonal->source();

                final_res += calc(*(lastPlaced), *next(lastPlaced), newCorner.direction, splittingDiagonal);

                lastPlaced = (eventList.insert(next(lastPlaced), newCorner));

                prev(lastPlaced)->direction = newCorner.direction;
                lastPlaced->direction = newCorner.direction;

                if (target_in_list != eventList.end()) {
                    target_in_list->top = newCorner.top;
                    target_in_list->bot = newCorner.bot;
                }

                lastPlaced->bot = prev(lastPlaced)->bot;

            } else if (find_if(next(crossingDiagonal), crossingDiagonals.end(),
                               [&](const Segment_2 &A) {
                                   return (A.target() == crossingDiagonal->target() ||
                                           A.source() == crossingDiagonal->target()) &&
                                          A.opposite() != *crossingDiagonal;
                               }) == crossingDiagonals.end() &&
                       crossingDiagonal->target() != splittingDiagonal.target() &&
                       crossingDiagonal->target() != splittingDiagonal.source()) {

                final_res += calc(*prev(target_in_list), *(target_in_list), newCorner.direction, splittingDiagonal);
                final_res += calc(*(target_in_list), *next(target_in_list, 1), newCorner.direction, splittingDiagonal);

                lastPlaced = prev(target_in_list);

                prev(target_in_list)->direction = newCorner.direction;
                next(target_in_list)->direction = newCorner.direction;

                next(target_in_list)->top = prev(target_in_list)->top;

                eventList.erase(target_in_list);

            } else if (find_if(next(crossingDiagonal), crossingDiagonals.end(),
                               [&crossingDiagonal](const Segment_2 &A) {
                                   return (A.target() == crossingDiagonal->source() ||
                                           A.source() == crossingDiagonal->source()) &&
                                          A.opposite() != *crossingDiagonal;
                               }) == crossingDiagonals.end() &&
                       crossingDiagonal->source() != splittingDiagonal.target() &&
                       crossingDiagonal->source() != splittingDiagonal.source()) {

                final_res += calc(*prev(source_in_list), *(source_in_list), newCorner.direction, splittingDiagonal);
                final_res += calc(*(source_in_list), *next(source_in_list, 1), newCorner.direction, splittingDiagonal);

                lastPlaced = prev(source_in_list);

                lastPlaced->direction = newCorner.direction;
                next(lastPlaced)->direction = newCorner.direction;

                source_in_list->bot = next(source_in_list)->bot;

                eventList.erase(source_in_list);
            } else {
                final_res += calc(*prev(lastPlaced, 2), *prev(lastPlaced), newCorner.direction, splittingDiagonal);
                final_res += calc(*prev(lastPlaced), *(lastPlaced), newCorner.direction, splittingDiagonal);
                final_res += calc(*(lastPlaced), *next(lastPlaced), newCorner.direction, splittingDiagonal);

//                prev(target_in_list, 2)->direction = newCorner.direction;
                prev(target_in_list)->direction = newCorner.direction;
                target_in_list->direction = newCorner.direction;
                next(target_in_list)->direction = newCorner.direction;


                target_in_list->bot = next(source_in_list)->bot;
                source_in_list->top = target_in_list->top;

                iter_swap(target_in_list, source_in_list);
                lastPlaced = prev(lastPlaced);
            }
        }
        final_res += calc(*(lastPlaced), *next(lastPlaced), f.direction(), splittingDiagonal);
        cout << "Sweep: " << final_res << endl;
        return final_res;
    }

    // compute contribution from trapezoid
    FT calc(const Corner &c1, const Corner &c2, const Direction_2 &newDir, const Segment_2 &diagonal) const {
        auto dia = std::min(diagonal.direction(), diagonal.opposite().direction());

        FT t = dia.dy() == 0 ? 0 : atan(dia.dy() / dia.dx());
        typename K::Aff_transformation_2 rotate(CGAL::ROTATION, sin(-t), cos(-t));

        auto PL = rotate(c1.target);
        auto PR = rotate(c2.target);

        auto TOP0 = rotate(c2.top.source());
        auto TOP1 = rotate(c2.top.target());
        auto BOT0 = rotate(c2.bot.source());
        auto BOT1 = rotate(c2.bot.target());

        auto V0 = rotate(c1.direction);
        auto V1 = rotate(newDir);

        typename K::Aff_transformation_2 trans(CGAL::TRANSLATION, typename K::Vector_2{-PL.x(), 0});

        PL = trans(PL);
        PR = trans(PR);
        TOP0 = trans(TOP0);
        TOP1 = trans(TOP1);
        BOT0 = trans(BOT0);
        BOT1 = trans(BOT1);

        return 2 * Integral<FT>::ComputeOpt(TOP0.x(), TOP0.y(), TOP1.x(), TOP1.y(),
                                            BOT0.x(), BOT0.y(), BOT1.x(), BOT1.y(),
                                            0, PL.y(), PR.x(), PR.y(),
                                            V0.dx(), V0.dy(), V1.dx(), V1.dy(),
                                            0.0, 1.0);
    }


    FT decompose_tree_rec(typename BinaryTree<CTP>::Node *node, const DiagonalMap &diagonals) {
        auto leftNode = tree.EmptyNode();
        auto rightNode = tree.EmptyNode();

        CentroidDecomposition<Face_handle> centroidDecomposition{node->data};
        auto centroid = std::find_if(node->data.begin(), node->data.end(),
                                     [&centroidDecomposition](const Face_handle &A) {
                                         return A->info().id == centroidDecomposition.id;
                                     });
        FT res = 0;
        if (node->data.size() == 1) {
            return 2 * triangulation.triangle(node->data.front()).area();
        } else {
            for (auto face: node->data) {
                if (centroidDecomposition.left.count(face->info().id) > 0) {
                    leftNode->data.emplace_back(face);
                } else if (centroidDecomposition.right.count(face->info().id) > 0) {
                    rightNode->data.emplace_back(face);
                } else {
                    throw std::logic_error("Face should belong to left or right set");
                }
            }

            node->left = leftNode;
            node->right = rightNode;
            auto diagonal = find_diagonal(*centroid, rightNode);

            auto diagonalSets = find_diagonal_subsets(leftNode, rightNode, diagonals, diagonal);
            
            res += decompose_tree_rec(leftNode, get<0>(diagonalSets));
            res += decompose_tree_rec(rightNode, get<2>(diagonalSets));

            res += create_event_list(leftNode, rightNode, get<1>(diagonalSets), diagonal);
            return res;
        }
    }
};

#endif //VISIBILITY_2_EXAMPLES_EXACT_CONVEXITY_MEASURE_2_H
