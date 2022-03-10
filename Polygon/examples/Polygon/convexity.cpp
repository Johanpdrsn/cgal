#include <CGAL/partition_2.h>
#include <CGAL/Simple_polygon_visibility_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/point_generators_2.h>
#include <fstream>

//#include <CGAL/Cartesian.h>
//typedef CGAL::Cartesian<double> Kernel;
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polygon_2<Kernel> Polygon_2;
typedef Kernel::Triangle_2 Triangle_2;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Segment_2 Segment_2;

using std::cout;
using std::endl;
using std::vector;

vector<Triangle_2> triangulate(Polygon_2 &poly) {
    typedef CGAL::Partition_traits_2<Kernel> Traits_T;
    typedef Traits_T::Polygon_2 Polygon_T_2;
    typedef std::list<Polygon_T_2> Polygon_list;

    Polygon_list partition_polys;

    if (!poly.is_counterclockwise_oriented()) {
        poly.reverse_orientation();
        CGAL::approx_convex_partition_2(poly.vertices_begin(),
                                        poly.vertices_end(), back_inserter(partition_polys));
        poly.reverse_orientation();

    } else {
        CGAL::approx_convex_partition_2(poly.vertices_begin(),
                                        poly.vertices_end(), back_inserter(partition_polys));
    }

    vector<Triangle_2> res{};
    for (const auto &fac: partition_polys) {
        res.emplace_back(Triangle_2{fac.vertex(0), fac.vertex(1), fac.vertex(2)});
    }

    return res;
}

Kernel::FT two_point_visibility_sample(const Polygon_2 &polygon,
                                       CGAL::Random_points_in_triangles_2<Point_2> &triangle_point_generator) {

    const Point_2 point_sample1 = *triangle_point_generator++;
    const Point_2 point_sample2 = *triangle_point_generator++;

    const Segment_2 seg{point_sample1, point_sample2};

    for (auto it = polygon.edges_begin(); it != polygon.edges_end(); it++) {

        if (CGAL::do_intersect(seg, it.make_value_type(CGAL::Tag_true()))) {
            return 0;
        }
    }

    return 1;
}

Kernel::FT visibility_polygon_sample(const Polygon_2 &polygon,
                                     CGAL::Random_points_in_triangles_2<Point_2> &point_generator) {

    typedef CGAL::Arr_segment_traits_2<Kernel> Traits_2;
    typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;
    typedef CGAL::Simple_polygon_visibility_2<Arrangement_2, CGAL::Tag_true> RSPV;

    vector<Segment_2> p{polygon.edges_begin(), polygon.edges_end()};

    Point_2 point_sample = *point_generator++;

    Arrangement_2 env;
    CGAL::insert_non_intersecting_curves(env, p.begin(), p.end());
    // find the face of the query point
    // (usually you may know that by other means)
    Arrangement_2::Face_const_handle *face;
    CGAL::Arr_naive_point_location<Arrangement_2> pl{env};
    CGAL::Arr_point_location_result<Arrangement_2>::Type obj = pl.locate(point_sample);
    // The query point locates in the interior of a face
    face = boost::get<Arrangement_2::Face_const_handle>(&obj);

    // compute regularized visibility_sample area
    // Define visibility object type that computes non-regularized visibility_sample area
    Arrangement_2 regular_output;
    RSPV regular_visibility{env};
    regular_visibility.compute_visibility(point_sample, *face, regular_output);

    vector<Point_2> visible_points;
    for (const auto &vertex: regular_output.vertex_handles()) {
        visible_points.emplace_back(vertex->point());
    }

    Kernel::FT area{};
    CGAL::area_2(visible_points.begin(), visible_points.end(), area);

    return CGAL::abs(area / polygon.area());
}

void read_polygon(const std::string &fileName, Polygon_2 &polygon) {
    std::ifstream in{fileName};

    if (!in.is_open()) {
        throw std::runtime_error("Could not open polygon file: " + fileName);
    }

    std::string message;

    int n;
    in >> n;

    Kernel::FT x{};
    Kernel::FT y{};

    for (int i = 1; i < n; i++) {
        in >> x;
        in >> y;
        polygon.push_back(Point_2{x, y});
    }

    in.close();
}


Kernel::FT
simulate(const std::function<Kernel::FT(const Polygon_2 &, CGAL::Random_points_in_triangles_2<Point_2> &)> &vis_func,
         const std::string &fileName, int n = 10000) {

//    Polygon_2 test_polygon;
//    test_polygon.push_back(Point_2{0, 0});
//    test_polygon.push_back(Point_2{-1, -1});
//    test_polygon.push_back(Point_2{2, -1});
//    test_polygon.push_back(Point_2{-1, 2});


    Polygon_2 file_polygon;
    read_polygon(fileName, file_polygon);

    Polygon_2 *polygon = &file_polygon;

    const vector<Triangle_2> triangles = triangulate(*polygon);
    CGAL::Random_points_in_triangles_2<Point_2> triangle_point_generator{triangles};

    Kernel::FT sum = 0;

    for (int i = 0; i < n; i++) {
        sum += vis_func(*polygon, triangle_point_generator);
    }

    return sum / n;
}


int main() {
    const int N = 10000;
    const std::string fileName{"rpg.line"};

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
