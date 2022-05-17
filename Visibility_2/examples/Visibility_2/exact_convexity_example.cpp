//
// Created by Johan Pedersen on 25/04/2022.
//

#include <iostream>
#include <cmath>
#include <CGAL/Cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Visibility_2/BinaryTree.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Rotational_sweep_visibility_2.h>
#include <CGAL/aff_transformation_tags.h>
#include <CGAL/Visibility_2/NumericalIntegral.h>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <CGAL/Aff_transformation_2.h>

//typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Cartesian<double> Kernel;
typedef CGAL::Polygon_2<Kernel> Polygon_2;
typedef Kernel::Point_2 Point_2;

typedef CGAL::Polygon_2<Kernel> Polygon_2;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Triangle_2 Triangle_2;
typedef Kernel::Segment_2 Segment_2;
typedef CGAL::Arr_segment_traits_2<Kernel> Traits_2;
typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;
typedef CGAL::Rotational_sweep_visibility_2<Arrangement_2> RSV;
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


private:
    struct PointHashFunction {
        size_t operator()(const CTP::Point &t) const {
            return hash<double>()(CGAL::to_double(CGAL::exact(t.x()))) ^
                   hash<double>()(CGAL::to_double(CGAL::exact(t.y())));
        }
    };

    struct SegmentHash {
        size_t operator()(const Segment_2 &t) const {
            auto h = hash<double>{};
            return h(CGAL::to_double((t.source().x()))) ^ h(CGAL::to_double((t.source().y())))
                   ^ h(CGAL::to_double((t.target().x()))) ^ h(CGAL::to_double((t.target().y())));
        }
    };

    struct Corner {
        Point_2 target;
        Point_2 source;
        Segment_2 top;
        Segment_2 bot;
        Kernel::Direction_2 direction{};

        bool operator==(const Corner &A) const {
            return this->target == A.target;
        }
    };


    static Kernel::FT calc(const Corner &c1, const Corner &c2, const Kernel::Direction_2 &newDir) {
        Kernel::FT PLX = c1.target.x();
        Kernel::FT PLY = c1.target.y();

        auto diffX = PLX;
        PLX = 0;


        Kernel::FT PRX = c2.target.x() - diffX;
        Kernel::FT PRY = c2.target.y();

        Kernel::FT TOP0X = c2.top.source().x() - diffX;
        Kernel::FT TOP0Y = c2.top.source().y();
        Kernel::FT TOP1X = c2.top.target().x() - diffX;
        Kernel::FT TOP1Y = c2.top.target().y();

        Kernel::FT BOT0X = c2.bot.source().x() - diffX;
        Kernel::FT BOT0Y = c2.bot.source().y();
        Kernel::FT BOT1X = c2.bot.target().x() - diffX;
        Kernel::FT BOT1Y = c2.bot.target().y();


        Kernel::FT V0X = c1.direction.dx();
        Kernel::FT V0Y = c1.direction.dy();
        Kernel::FT V1X = newDir.dx();
        Kernel::FT V1Y = newDir.dy();


        using boost::math::quadrature::trapezoidal;
        auto f = [&](Kernel::FT rho) {
            if (rho == 0.0) {
                return rho;
            } else {
                return CGAL::abs(
                        Integral<Kernel::FT>::Compute(TOP0X, TOP0Y, TOP1X, TOP1Y, BOT0X, BOT0Y, BOT1X, BOT1Y, PLX,
                                                      PLY,
                                                      PRX, PRY, V0X, V0Y, V1X, V1Y, rho));
            }

        };

        cout << "TOP0: " << "(" << TOP0X << "," << TOP0Y << ")" << endl;
        cout << "TOP1: " << "(" << TOP1X << "," << TOP1Y << ")" << endl;
        cout << "BOT0: " << "(" << BOT0X << "," << BOT0Y << ")" << endl;
        cout << "BOT1: " << "(" << BOT1X << "," << BOT1Y << ")" << endl;
        cout << "PL: " << "(" << PLX << "," << PLY << ")" << endl;
        cout << "PR: " << "(" << PRX << "," << PRY << ")" << endl;
        cout << "V0: " << "(" << V0X << "," << V0Y << ")" << endl;
        cout << "V1: " << "(" << V1X << "," << V1Y << ")" << endl;

        Kernel::FT intRes = trapezoidal(f, 0.0, 1.0-1e-10);
//        Kernel::FT intRes = Integral<Kernel::FT>::Compute(TOP0X, TOP0Y, TOP1X, TOP1Y, BOT0X, BOT0Y, BOT1X, BOT1Y, PLX,
//                                                          PLY, PRX, PRY, V0X, V0Y, V1X, V1Y);
        cout << intRes << endl;
        return intRes;
    }

    static bool face_equality(const CTP::Face_handle &A, const CTP::Face_handle &B) {
        std::unordered_set<CTP::Point, PointHashFunction> s{A->vertex(0)->point(), A->vertex(1)->point(),
                                                            A->vertex(2)->point()};
        for (int i = 0; i < 3; i++) {
            if (s.count(B->vertex(i)->point()) == 0)
                return false;
        }
        return true;
    }

    static bool face_in_list(const vector<CTP::Face_handle> &list, const CTP::Face_handle &face) {
        return std::find_if(list.begin(), list.end(),
                            [&](const CTP::Face_handle &A) { return face_equality(A, face); }) != list.end();
    }

    static double count_vertices(const vector<CTP::Face_handle> &list) {
        std::unordered_set<CTP::Point, PointHashFunction> s;
        for (auto face: list) {
            for (int i = 0; i < 3; i++) {
                s.insert(face->vertex(i)->point());
            }
        }
        return static_cast<double>(s.size());
    }

    static bool is_edge_in_faces(vector<CTP::Face_handle> &faces, Segment_2 &seg) {
        for (auto face: faces) {
            for (int i = 0; i < 3; i++) {
                if (seg == Segment_2(face->vertex(i)->point(), face->vertex((i + 1) % 3)->point()))
                    return true;
            }
        }
        return false;
    }


public:

    Polygon_2 polygon;
    vector<CTP::Face_handle> triangles;

    CTP triangulation;
    BinaryTree<CTP> tree;
    Segment_2 diagonal;
    unordered_map<Segment_2, Segment_2, SegmentHash> diagMap;

    explicit Convexity_measure_exact_2(Polygon_2 &input_polygon) : polygon(input_polygon), tree() {
        triangulate();
        sweep_diagonals();
        generate_tree(triangles);

        create_event_list();
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
                triangles.emplace_back(face);
            }
        }
    }

    Segment_2 *find_diagonal(BinaryTree<CTP>::Node *leftNode, BinaryTree<CTP>::Node *rightNode) {
        int res;
        for (auto leftFace: leftNode->data) {
            for (auto rightFace: rightNode->data) {
                if (leftFace->has_neighbor(rightFace, res)) {
                    return new Segment_2(triangulation.segment(CTP::Edge(leftFace, res)));
                }
            }
        }
        return nullptr;
    }


    tuple<vector<Segment_2>, deque<Segment_2>, vector<Segment_2>>
    find_diagonal_subsets(BinaryTree<CTP>::Node *leftNode, BinaryTree<CTP>::Node *rightNode) {
        vector<Segment_2> leftDiagonals;
        vector<Segment_2> rightDiagonals;
        deque<Segment_2> crossingDiagonals;


        unordered_set<Point_2, PointHashFunction> leftPoints;
        unordered_set<Point_2, PointHashFunction> rightPoints;

        diagonal = *find_diagonal(leftNode, rightNode);

        for (auto f: leftNode->data) {
            for (int i = 0; i < 3; i++) {
                leftPoints.insert(f->vertex(i)->point());
            }
        }
        for (auto f: rightNode->data) {
            for (int i = 0; i < 3; i++) {
                rightPoints.insert(f->vertex(i)->point());
            }
        }

        for (const auto &diag: diagMap) {
            if ((leftPoints.count(diag.first.source()) == 1 && rightPoints.count(diag.first.target()) == 1) ||
                (rightPoints.count(diag.first.source()) == 1 && leftPoints.count(diag.first.target())) == 1) {
                crossingDiagonals.emplace_back(diag.first);
            } else if (leftPoints.count(diag.first.source()) == 1 && leftPoints.count(diag.first.target()) == 1) {
                leftDiagonals.emplace_back(diag.first);
            } else if (rightPoints.count(diag.first.source()) == 1 && rightPoints.count(diag.first.target()) == 1) {
                rightDiagonals.emplace_back(diag.first);
            }
        }

        for (const auto &s: leftDiagonals) {
            if (!is_edge_in_faces(leftNode->data, diagMap[s]))
                diagMap[s] = diagonal;
        }

        for (const auto &s: rightDiagonals) {
            if (!is_edge_in_faces(rightNode->data, diagMap[s]))
                diagMap[s] = diagonal;
        }

        cout << "DIAG: " << diagonal << endl;

        std::sort(leftDiagonals.begin(), leftDiagonals.end(),
                  [](const Segment_2 &A, const Segment_2 &B) { return A.direction() < B.direction(); }
        );
        std::sort(rightDiagonals.begin(), rightDiagonals.end(),
                  [](const Segment_2 &A, const Segment_2 &B) { return A.direction() < B.direction(); }
        );
        std::sort(crossingDiagonals.begin(), crossingDiagonals.end(),
                  [this](const Segment_2 &A, const Segment_2 &B) {
                      if (A == diagonal || A == diagonal.opposite())
                          return true;
                      else if (B == diagonal || B == diagonal.opposite())
                          return false;
                      else {
                          auto d = std::max(diagonal.direction(), diagonal.opposite().direction());


                          auto AMin = std::min(A.direction(), A.opposite().direction());
                          auto AMax = std::max(A.direction(), A.opposite().direction());

                          auto BMin = std::min(B.direction(), B.opposite().direction());
                          auto BMax = std::max(B.direction(), B.opposite().direction());


                          return AMin.counterclockwise_in_between(d, BMin) || AMax.counterclockwise_in_between(d, BMax);
                      }
                  }
        );

        return make_tuple(leftDiagonals, crossingDiagonals, rightDiagonals);
    }

    void decompose_tree_rec(BinaryTree<CTP>::Node *node) {
        auto leftNode = tree.EmptyNode();
        auto rightNode = tree.EmptyNode();

        auto n = count_vertices(node->data);
        auto lowerBound = std::floor((n - 1) / 3.0);
        auto upperBound = std::floor((2 * n - 5) / 3.0);

        if (node->data.size() <= 1) {
            return;
        } else {
            for (auto face: node->data) {
                leftNode->data.emplace_back(face);
                face->tds_data().mark_processed();
                while (true) {

                    // Find face that is in the input list and has not been added to the left node
                    for (int i = 0; i < 3; i++) {
                        face = face->neighbor(i);
                        if (face_in_list(node->data, face) && !face->tds_data().processed())
                            break;
                    }

                    // If the size of the leftNode is smaller than bound add it and mark the face
                    if (static_cast<double>(leftNode->data.size()) + 1 <= upperBound) {
                        face->tds_data().mark_processed();
                        leftNode->data.emplace_back(face);
                    } else {
                        break;
                    }
                }

                // If the left node list is larger than the lower bound add unmark faced to right
                if (static_cast<double>(leftNode->data.size()) >= lowerBound) {
                    for (auto f: node->data) {
                        f->tds_data().clear();
                        if (!face_in_list(leftNode->data, f))
                            rightNode->data.emplace_back(f);
                    }
                    break;
                } else {
                    leftNode->data.clear();
                    rightNode->data.clear();
                }
            }
            node->left = leftNode;
            decompose_tree_rec(leftNode);
            node->right = rightNode;
            decompose_tree_rec(rightNode);
        }
    }

    void generate_tree(vector<CTP::Face_handle> &data) {
        tree.SetRoot(data);
        decompose_tree_rec(tree.Root());
    }


    void sweep_diagonals() {
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

                    // CHeck if point in the visible polygon exists in the input polygon
                    if (vis_point_in_poly && p != v->point()) {
                        auto seg = Segment_2(p, v->point());
                        diagMap.insert(make_tuple(seg, find_eta(seg)));
                    }
                }
            }
        }
    }

    Segment_2 find_nearest_intersecting_segment(const Segment_2 &s) const {
        auto r = Kernel::Ray_2(s.target(), s.direction());
        Kernel::FT min_y = INT_MAX;
        Segment_2 left_plus;

        for (auto it = polygon.edges_begin(); it != polygon.edges_end(); it++) {
            const auto result = intersection(r, *it);
            if (result) {
                if (const Point_2 *p = boost::get<Point_2>(&*result)) {
                    Kernel::FT diff = abs(p->y() - s.target().y());
                    if (diff > 0 && diff < min_y) {
                        min_y = diff;
                        left_plus = *it;
                    }
                }
            }
        }
        return left_plus;
    }

    Segment_2 find_eta(Segment_2 &s) const {
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

    void create_event_list() {

        auto current_diagonals = find_diagonal_subsets(tree.Root()->left, tree.Root()->right);
//        auto leftDiagonals = get<0>(current_diagonals);
        auto crossingDiagonals = get<1>(current_diagonals);
//        auto rightDiagonals = get<2>(current_diagonals);

//
//        for (auto a: diagMap) {
//            cout << a.first << " : " << a.second << endl;
//        }

//        for (auto a: crossingDiagonals) {
//            cout << a << endl;
//        }

        Kernel::FT final_res = 0;
        auto eventList = std::list<Corner>();


        auto f = crossingDiagonals.front();
        crossingDiagonals.pop_front();
        auto b = crossingDiagonals.front();
        crossingDiagonals.pop_front();

        if (f.direction() < b.direction()) {
            std::swap(f, b);
        }


        eventList.push_front({f.target(), f.source(), diagMap[b], diagMap[f], b.direction()});
        eventList.push_back({b.target(), b.source(), diagMap[b], diagMap[f], b.direction()});


        auto lastPlaced = eventList.begin();

        for (auto crossingDiagonal = crossingDiagonals.begin();
             crossingDiagonal != crossingDiagonals.end(); crossingDiagonal++) {


            auto s = crossingDiagonal->direction() < crossingDiagonal->opposite().direction() ? *crossingDiagonal
                                                                                              : crossingDiagonal->opposite();


            Corner newCorner = {crossingDiagonal->target(), crossingDiagonal->source(), diagMap[s],
                                diagMap[s.opposite()],
                                s.direction()};
            auto target_in_list = find(eventList.begin(), eventList.end(), newCorner);
            auto source_in_list = find(eventList.begin(), eventList.end(), Corner{newCorner.source});


            bool onTop = crossingDiagonal->direction().dy() > 0;


            if (eventList.size() < 2) {
                if (newCorner.direction.dx() > 0) {
                    lastPlaced = eventList.insert(next(lastPlaced), newCorner);
                } else {
                    lastPlaced = eventList.insert(lastPlaced, newCorner);
                }
            } else {
                if (target_in_list == eventList.end()) {
                    cout << "APPEARANCE: " << newCorner.target << endl;

                    for (auto &a: eventList) {
                        cout << a.target << " : " << a.top << " : " << a.bot << endl;
                    }

                    cout << "Last placed: " << lastPlaced->target << endl;
                    final_res += calc(*(lastPlaced), *next(lastPlaced), newCorner.direction);


                    lastPlaced = (eventList.insert(next(lastPlaced), newCorner));

                    prev(lastPlaced)->direction = newCorner.direction;
                    lastPlaced -> direction = newCorner.direction;

                    if (source_in_list != eventList.end()) {
                        source_in_list->top = newCorner.top;
                        source_in_list->bot = newCorner.bot;
                    }

                    if (onTop) {
                        prev(lastPlaced)->top = lastPlaced->top;
                    } else {
                        lastPlaced->bot = prev(lastPlaced)->bot;
                    }


                    for (auto &a: eventList) {
                        cout << a.target << " : " << a.top << " : " << a.bot << endl;
                    }
                    cout << "Last placed: " << lastPlaced->target << endl;

                } else if (source_in_list != eventList.end()) {

                    if (crossingDiagonal->source().y() < 0 && crossingDiagonal->target().y() > 0) {
                        cout << "SWAP:" << crossingDiagonal->source() << "<->" << crossingDiagonal->target() << endl;
                        for (auto &a: eventList) {
                            cout << a.target << " : " << a.top << " : " << a.bot << endl;
                        }
                        cout << "Last placed: " << lastPlaced->target << endl;
                        cout << "TARGET IN LIST: " << target_in_list->target << endl;

                        final_res += calc(*prev(lastPlaced, 2), *prev(lastPlaced), newCorner.direction);
                        final_res += calc(*prev(lastPlaced), *(lastPlaced), newCorner.direction);
                        final_res += calc(*(lastPlaced), *next(lastPlaced), newCorner.direction);

                        if (target_in_list->target.y() >= next(target_in_list)->target.y()){
                            (target_in_list)->top = prev(target_in_list)->top;
                            next(target_in_list, 0)->bot = next(target_in_list, 2)->bot;
                        } else{
                            target_in_list->bot = next(target_in_list,2)->bot;
                            next(target_in_list)->top = target_in_list->top;
                        }

//                        std::swap(target_in_list->direction, next(target_in_list)->direction);
                        prev(target_in_list,2)-> direction = newCorner.direction;
                        prev(target_in_list)-> direction = newCorner.direction;
                        target_in_list->direction = newCorner.direction;
                        next(target_in_list)->direction = newCorner.direction;

                        iter_swap(target_in_list, source_in_list);
                        lastPlaced = prev(lastPlaced);

                        for (auto &a: eventList) {
                            cout << a.target << " : " << a.top << " : " << a.bot << " : " << a.direction << endl;
                        }
                        cout << "Last placed: " << lastPlaced->target << endl;
                    } else if (count_if(next(crossingDiagonal), crossingDiagonals.end(),
                                        [&crossingDiagonal](const Segment_2 &A) {
                                            return A.target() == crossingDiagonal->target();
                                        }) == 0 && (crossingDiagonal->target() != diagonal.target() &&
                                                    crossingDiagonal->target() != diagonal.source())) {
                        cout << "DISAPPEARANCE: " << crossingDiagonal->target() << endl;


                        for (auto &a: eventList) {
                            cout << a.target << " : " << a.top << " : " << a.bot << " : " << a.direction << endl;
                        }

                        final_res += calc(*prev(target_in_list), *(target_in_list), newCorner.direction);
                        final_res += calc(*(target_in_list), *next(target_in_list, 1), newCorner.direction);

                        lastPlaced = prev(target_in_list);


                        lastPlaced->direction = newCorner.direction;
                        next(lastPlaced)->direction = newCorner.direction;


                        if (onTop) {
                            prev(target_in_list)->top = target_in_list->top;
                            next(target_in_list)->top = newCorner.top;
                        } else
                            target_in_list->bot = next(target_in_list)->bot;


                        eventList.erase(target_in_list);

                        for (auto &a: eventList) {
                            cout << a.target << " : " << a.top << " : " << a.bot << " : " << a.direction << endl;
                        }

                        cout << "Last placed: " << lastPlaced->target << endl;
                    }
                }
            }
        }
        cout << "END" << endl;
        final_res += calc(*(lastPlaced), *next(lastPlaced), f.direction());

        cout << "Final: " << final_res << endl;
    }
};


int main() {
    CGAL::IO::set_pretty_mode(std::cout);

    Polygon_2 polygon;
//    polygon.push_back(Point_2{0, 0});
//    polygon.push_back(Point_2{-2, -1});
//    polygon.push_back(Point_2{0, -2});
//    polygon.push_back(Point_2{1, -4});
//    polygon.push_back(Point_2{2, -2});
//    polygon.push_back(Point_2{4, -1});
//    polygon.push_back(Point_2{2, 0});
//    polygon.push_back(Point_2{1, 2});

    polygon.push_back(Point_2{0, 0});
    polygon.push_back(Point_2{1, -2.0});
    polygon.push_back(Point_2{2, 0});
    polygon.push_back(Point_2{2, 2});
//
//    polygon.push_back(Point_2{0.0, 0.0});
//    polygon.push_back(Point_2{-2.0, -1.0});
//    polygon.push_back(Point_2{1.0, -1.0});
//    polygon.push_back(Point_2{1.0, 2.0});


    auto test = Convexity_measure_exact_2(polygon);


//    auto a = Point_2(1, 0);
//
//    Kernel::Aff_transformation_2 rational_rotate(CGAL::ROTATION, Kernel::Direction_2(2, 1), 1, 100);
//    Kernel::Aff_transformation_2 rotate(CGAL::ROTATION, sin(CGAL_PI), cos(CGAL_PI));
//
//    auto b = rational_rotate(a);
//
//    auto c = rotate(a);
//
//    cout << a << endl;
//    cout << b << endl;
//    cout << c << endl;



//    test.tree.prettyPrint();
//

    return 0;
}