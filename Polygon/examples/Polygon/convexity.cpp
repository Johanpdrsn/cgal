#include <CGAL/Polygon_2.h>
#include <CGAL/Triangular_expansion_visibility_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <fstream>
#include <random>

//#include <CGAL/Cartesian.h>
//typedef CGAL::Cartesian<double> Kernel;
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polygon_2<Kernel> Polygon_2;
typedef Kernel::Triangle_2 Triangle_2;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Segment_2 Segment_2;


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

using std::cout;
using std::endl;

void mark_domains(CT &ct,
                  CT::Face_handle start,
                  int index,
                  std::list<CT::Edge> &border) {
    if (start->info().nesting_level != -1) {
        return;
    }
    std::list<CT::Face_handle> queue;
    queue.push_back(start);
    while (!queue.empty()) {
        CT::Face_handle fh = queue.front();
        queue.pop_front();
        if (fh->info().nesting_level == -1) {
            fh->info().nesting_level = index;
            for (int i = 0; i < 3; i++) {
                CT::Edge e(fh, i);
                CT::Face_handle n = fh->neighbor(i);
                if (n->info().nesting_level == -1) {
                    if (ct.is_constrained(e)) border.push_back(e);
                    else queue.push_back(n);
                }
            }
        }
    }
}

//explore set of facets connected with non-constrained edges,
//and attribute to each such set a nesting level.
//We start from facets incident to the infinite vertex, with a nesting
//level of 0. Then we recursively consider the non-explored facets incident
//to constrained edges bounding the former set and increase the nesting level by 1.
//Facets in the domain are those with an odd nesting level.
void mark_domains(CT &cdt) {
    for (CT::Face_handle f: cdt.all_face_handles()) {
        f->info().nesting_level = -1;
    }
    std::list<CT::Edge> border;
    mark_domains(cdt, cdt.infinite_face(), 0, border);
    while (!border.empty()) {
        CT::Edge e = border.front();
        border.pop_front();
        CT::Face_handle n = e.first->neighbor(e.second);
        if (n->info().nesting_level == -1) {
            mark_domains(cdt, n, e.first->info().nesting_level + 1, border);
        }
    }
}

std::vector<Triangle_2> triangulate(Polygon_2 &polygon) {

    CT ct;
    ct.insert_constraint(polygon.vertices_begin(), polygon.vertices_end(), true);
    mark_domains(ct);

    std::vector<Triangle_2> triangles;
    for (auto var: ct.finite_face_handles()) {
        if (var->info().in_domain()) {
            triangles.emplace_back(ct.triangle(var));
        }
    }

    return triangles;
}

void read_polygon(const std::string &fileName, Polygon_2 &polygon) {
    std::ifstream in{fileName};

    if (!in.is_open()) {
        throw std::runtime_error("Could not open polygon file: " + fileName);
    }

    Kernel::FT x{};
    Kernel::FT y{};

    int n;
    in >> n;

    for (int i = 1; i < n; i++) {
        in >> x;
        in >> y;
        polygon.push_back(Point_2{x, y});
    }
    in.close();
}

Point_2 kraemer_method(CGAL::Triangle_2<Kernel> &T) {

    std::random_device rand;
    std::mt19937 generator(rand());
    std::uniform_real_distribution<double> distribution(0, 1);

    double a = distribution(generator);
    double b = distribution(generator);

    double q = CGAL::abs(a - b);

    double t = 0.5 * (a + b - q);
    double u = 1 - 0.5 * (q + a + b);

    return Point_2{q * T[0].x() + t * T[1].x() + u * T[2].x(), q * T[0].y() + t * T[1].y() + u * T[2].y()};
}


Kernel::FT two_point_visibility_sample(const Polygon_2 &polygon,
                                       CGAL::Random_points_in_triangles_2<Point_2> &triangle_point_generator,
                                       const Kernel::FT &n) {

    Kernel::FT sum = 0;
    for (int i = 0; i < n; i++) {
        Point_2 p1 = *triangle_point_generator++;
        Point_2 p2 = *triangle_point_generator++;

//        Triangle_2 A3 = Triangle_2{Point_2{1, 2}, Point_2{0, 0}, Point_2{1, 0.5}};
//        Triangle_2 A1a = Triangle_2{Point_2{1, 0.5}, Point_2{0, 0}, Point_2{1, -1}};
//        Triangle_2 A1b = Triangle_2{Point_2{1, -1}, Point_2{0, 0}, Point_2{-0.5, -1}};
//        Triangle_2 A2 = Triangle_2{Point_2{0, 0}, Point_2{-0.5, -1}, Point_2{-2, -1}};
//        Point_2 p2 = kraemer_method(A3);

        const Segment_2 seg{p2, p1};

        for (auto it = polygon.edges_begin(); it != polygon.edges_end(); it++) {
            if (CGAL::do_intersect(seg, it.make_value_type(CGAL::Tag_true()))) {
                sum += 1;
                break;
            }
        }
    }

    return (n - sum) / n;
}

Kernel::FT visibility_polygon_sample(const Polygon_2 &polygon,
                                     CGAL::Random_points_in_triangles_2<Point_2> &point_generator,
                                     const Kernel::FT &n) {
    typedef CGAL::Arr_segment_traits_2<Kernel> Traits_2;
    typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;
    typedef CGAL::Triangular_expansion_visibility_2<Arrangement_2, CGAL::Tag_false> RSPV;

//    Triangle_2 A3 = Triangle_2{Point_2{1, 2}, Point_2{0, 0}, Point_2{1, 0.5}};
//    Triangle_2 A1a = Triangle_2{Point_2{1, 0.5}, Point_2{0, 0}, Point_2{1, -1}};
//    Triangle_2 A1b = Triangle_2{Point_2{1, -1}, Point_2{0, 0}, Point_2{-0.5, -1}};
//    Triangle_2 A2 = Triangle_2{Point_2{0, 0}, Point_2{-0.5, -1}, Point_2{-2, -1}};

    Arrangement_2 env;
    CGAL::insert(env, polygon.edges_begin(), polygon.edges_end());
    RSPV regular_visibility{env};

    Arrangement_2::Face_const_handle *face;
    CGAL::Arr_naive_point_location<Arrangement_2> pl{env};
    CGAL::Arr_point_location_result<Arrangement_2>::Type obj;
    Arrangement_2 regular_output;

    Kernel::FT sum = 0;

    for (int i = 0; i < n; i++) {
        // find the face of the query point
        // (usually you may know that by other means)
        Point_2 point_sample = *point_generator++;
//        Point_2 point_sample = kraemer_method(A3);

        obj = pl.locate(point_sample);
        // The query point locates in the interior of a face
        face = boost::get<Arrangement_2::Face_const_handle>(&obj);

        // compute regularized visibility_sample area
        // Define visibility object type that computes regularized visibility_sample area
        regular_visibility.compute_visibility(point_sample, *face, regular_output);

        std::vector<Point_2> visible_points;
        for (const auto &vertex: regular_output.vertex_handles()) {
            visible_points.emplace_back(vertex->point());
        }

        Kernel::FT area{};
        CGAL::area_2(visible_points.begin(), visible_points.end(), area);

        sum += CGAL::abs(area);
    }

    return sum / (polygon.area() * n);
}


Kernel::FT simulate(const std::function<Kernel::FT(
        const Polygon_2 &, CGAL::Random_points_in_triangles_2<Point_2> &, const Kernel::FT &)> &vis_func,
                    const std::string &fileName,
                    const Kernel::FT &n = 10000) {

    Polygon_2 polygon;

    polygon.push_back(Point_2{0, 0});
    polygon.push_back(Point_2{-2, -1});
    polygon.push_back(Point_2{1, -1});
    polygon.push_back(Point_2{1, 2});

//    read_polygon(fileName, polygon);

    const std::vector<Triangle_2> triangles = triangulate(polygon);
    CGAL::Random_points_in_triangles_2<Point_2> triangle_point_generator{triangles};

    return vis_func(polygon, triangle_point_generator, n);
}


int main() {
    const Kernel::FT N = 100000;
    const std::string fileName{"rpg.line"};

    CGAL::set_pretty_mode(cout);

    std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double, std::milli>> t1;
    std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double, std::milli>> t2;
    std::chrono::duration<double, std::milli> ms_double{};

    t1 = std::chrono::high_resolution_clock::now();
    Kernel::FT polygon = simulate(visibility_polygon_sample, fileName, N);
    t2 = std::chrono::high_resolution_clock::now();
    ms_double = t2 - t1;
    cout << "Polygon prob: " << polygon << " in " << ms_double.count() << "ms" << endl;

    t1 = std::chrono::high_resolution_clock::now();
    Kernel::FT two_points = simulate(two_point_visibility_sample, fileName, N);
    t2 = std::chrono::high_resolution_clock::now();
    ms_double = t2 - t1;
    cout << "Two points prob: " << two_points << " in " << ms_double.count() << "ms" << endl;

    return 0;
}
