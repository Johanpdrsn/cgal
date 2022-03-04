#include <CGAL/partition_2.h>
#include <CGAL/Simple_polygon_visibility_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Cartesian.h>

//typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Cartesian<double> Kernel;
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

    vector<Triangle_2> res;
    for (auto fac: partition_polys) {
        res.emplace_back(Triangle_2(fac.vertex(0), fac.vertex(1), fac.vertex(2)));
    }

    return res;
}

int two_point_visibility_sample(Polygon_2 &polygon,
                                CGAL::Random_points_in_triangles_2<Point_2> &triangle_point_generator) {

    Point_2 point_sample1 = *triangle_point_generator++;
    Point_2 point_sample2 = *triangle_point_generator++;

    Segment_2 seg(point_sample1, point_sample2);

    for (auto it = polygon.edges_begin(); it != polygon.edges_end(); it++) {
        if (CGAL::do_intersect(seg, it.make_value_type(CGAL::Tag_true()))) {
            return 0;
        }
    }

    return 1;
}

double visibility_polygon_sample(Polygon_2 &polygon,
                                 CGAL::Random_points_in_triangles_2<Point_2> &point_generator) {

    typedef CGAL::Arr_segment_traits_2<Kernel> Traits_2;
    typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;
    typedef CGAL::Simple_polygon_visibility_2<Arrangement_2, CGAL::Tag_false> NSPV;

    vector<Segment_2> p(polygon.edges_begin(), polygon.edges_end());

    Point_2 point_sample = *point_generator++;

    Arrangement_2 env;
    CGAL::insert_non_intersecting_curves(env, p.begin(), p.end());
    // find the face of the query point
    // (usually you may know that by other means)
    Arrangement_2::Face_const_handle *face;
    CGAL::Arr_naive_point_location<Arrangement_2> pl(env);
    CGAL::Arr_point_location_result<Arrangement_2>::Type obj = pl.locate(point_sample);
    // The query point locates in the interior of a face
    face = boost::get<Arrangement_2::Face_const_handle>(&obj);

    // compute non regularized visibility_sample area
    // Define visibility object type that computes non-regularized visibility_sample area
    Arrangement_2 non_regular_output;
    NSPV non_regular_visibility(env);
    non_regular_visibility.compute_visibility(point_sample, *face, non_regular_output);

    vector<Point_2> visible_points;
    for (auto vertex: non_regular_output.vertex_handles()) {
        visible_points.emplace_back(vertex->point());
    }

    double area;
    CGAL::area_2(visible_points.begin(), visible_points.end(), area);

    return CGAL::abs(area / polygon.area());
}


double simulate(int n, const std::function<double(Polygon_2 &,
                                                  CGAL::Random_points_in_triangles_2<Point_2> &)> &vis_func) {

    Polygon_2 polygon;
    polygon.push_back(Point_2(0, 0));
    polygon.push_back(Point_2(-2, -1));
    polygon.push_back(Point_2(1, -1));
    polygon.push_back(Point_2(1, 2));

    vector<Triangle_2> triangles = triangulate(polygon);
    CGAL::Random_points_in_triangles_2<Point_2> triangle_point_generator(triangles);

    double sum = 0;

    for (int i = 0; i < n; i++) {
        sum += vis_func(polygon, triangle_point_generator);
    }

    return sum / n;
}


int main() {
    int n = 1000000;

    std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double, std::milli>> t1;
    std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double, std::milli>> t2;
    std::chrono::duration<double, std::milli> ms_double{};

    t1 = std::chrono::high_resolution_clock::now();
    double polygon = simulate(n, visibility_polygon_sample);
    t2 = std::chrono::high_resolution_clock::now();
    ms_double = t2 - t1;
    cout << "Polygon prob: " << polygon << " in " << ms_double.count() << "ms" << endl;


    t1 = std::chrono::high_resolution_clock::now();
    double two_points = simulate(n, two_point_visibility_sample);
    t2 = std::chrono::high_resolution_clock::now();
    ms_double = t2 - t1;
    cout << "Two points prob: " << two_points << " in " << ms_double.count() << "ms" << endl;


    return 0;
}
