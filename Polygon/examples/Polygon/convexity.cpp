#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <iostream>
#include <c++/v1/random>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_polygon_visibility_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_naive_point_location.h>
using namespace std;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Partition_traits_2<K> Traits;
typedef CGAL::Point_2<K> Point;
typedef Traits::Point_2 Point_T_2;
typedef Traits::Polygon_2 Polygon_T_2;
typedef std::list<Polygon_T_2> Polygon_list;


typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Segment_2 Segment_2;
typedef CGAL::Arr_segment_traits_2<Kernel> Traits_2;
typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;
typedef Arrangement_2::Edge_const_iterator Edge_const_iterator;

typedef CGAL::Simple_polygon_visibility_2<Arrangement_2, CGAL::Tag_false> NSPV;


# define nice(os) ((os == CGAL::ON_ORIENTED_BOUNDARY) ? "on boundary" :  \
                   (os == CGAL::NEGATIVE) ? "inside" : "outside")


unsigned long triangle_index_search(double sample, std::vector<double> cdf) {
    unsigned long l = 0;
    unsigned long h = cdf.size();
    while (l <= h) {
        unsigned long m = l + (h - l) / 2;
        if (sample <= cdf[m]) {
            if (m == 0 || m > 0 && sample > cdf[m - 1]) {
                return m;
            }
            h = m - 1;
        } else {
            l = m + 1;
        }
    }
    return -1;
}

Point sample_triangle_point(CGAL::Triangle_2<K> &t) {

    random_device rand;
    mt19937 generator(rand());
    uniform_real_distribution<double> distribution(0, 1);

    double a = distribution(generator);
    double b = distribution(generator);

    Point v1 = Point(0, 3);
    Point v2 = Point(1, 1);

    auto res = Point{(v1.x() + v2.x()) * a, (v1.y() + v2.y()) * b};
    while (!t.bounded_side(res)) {
        a = distribution(generator);
        b = distribution(generator);
        res = Point{(v1.x() + v2.x()) * a, (v1.y() + v2.y()) * b};
    }

    return res;
}

Point kraemer_method(CGAL::Triangle_2<K> &T) {

    random_device rand;
    mt19937 generator(rand());
    uniform_real_distribution<double> distribution(0, 1);

    double a = distribution(generator);
    double b = distribution(generator);

    double q = abs(a - b);

    double s = q;
    double t = 0.5 * (a + b - q);
    double u = 1 - 0.5 * (q + a + b);

    return Point{s * T[0].x() + t * T[1].x() + u * T[2].x(), s * T[0].y() + t * T[1].y() + u * T[2].y()};
}

Point_2 trisample(CGAL::Triangle_2<K> &T) {

    random_device rand;
    mt19937 generator(rand());
    uniform_real_distribution<double> distribution(0, 1);

    double a = distribution(generator);
    double b = distribution(generator);

    double s1 = sqrt(a);

    return Point_2{T[0].x() * (1 - s1) + T[1].x() * (1 - b) * s1 + T[2].x() * b * s1,
                   T[0].y() * (1 - s1) + T[1].y() * (1 - b) * s1 + T[2].y() * b * s1};
}


std::vector<CGAL::Triangle_2<K>> triangulate(Polygon_2 &poly) {
    if (!poly.is_counterclockwise_oriented()) {
        poly.reverse_orientation();
    }
    Polygon_list partition_polys;

    CGAL::approx_convex_partition_2(poly.vertices_begin(), poly.vertices_end(), back_inserter(partition_polys));

    std::vector<CGAL::Triangle_2<K>> res;
    for (auto fac: partition_polys) {
        res.emplace_back(CGAL::Triangle_2<K>(fac.vertex(0), fac.vertex(1), fac.vertex(2)));
    }

    return res;
}


int main() {
    CGAL::IO::set_pretty_mode(std::cout);

    Polygon_2 polygon;
    polygon.push_back(Point_T_2(0, 0));
    polygon.push_back(Point_T_2(0, 3));
    polygon.push_back(Point_T_2(1, 1));
    polygon.push_back(Point_T_2(3, 0));

    std::vector<CGAL::Triangle_2<K>> triangles = triangulate(polygon);

    for (auto tri: triangles) {
        cout << tri << endl;
    }

    std::sort(triangles.begin(), triangles.end(),
              [](const auto &p1, const auto &p2) { return p1.area() < p2.area(); });
    std::vector<double> area_prefix_sum{triangles[0].area()};

    for (int i = 1; i < triangles.size(); i++) {
        area_prefix_sum.emplace_back(area_prefix_sum[i - 1] + triangles[i].area());
    }


    random_device rand;
    mt19937 generator(rand());
    uniform_real_distribution<double> distribution(0, area_prefix_sum[area_prefix_sum.size() - 1]);

    auto sample1 = distribution(generator);
    auto sample2 = distribution(generator);

    unsigned long triangle_index1 = triangle_index_search(sample1, area_prefix_sum);
    unsigned long triangle_index2 = triangle_index_search(sample2, area_prefix_sum);

    CGAL::Triangle_2<K> triangle1 = triangles[triangle_index1];
    CGAL::Triangle_2<K> triangle2 = triangles[triangle_index2];


    vector<Segment_2> p;

    for (auto it = polygon.vertex_pairs_begin(); it != polygon.vertex_pairs_end(); ++it) {
        auto x = Point_2(it->first.get().x(), it->first.get().y());
        auto y = Point_2(it->second.get().x(), it->second.get().y());

        p.emplace_back(Segment_2(x, y));
    }

    Point_2 point_sample = trisample(triangle1);
    Point_2 point_sample2 = trisample(triangle2);

    cout << "Sample: " << point_sample << endl;
    cout << "Sample: " << point_sample2 << endl;


    Arrangement_2 env;
    CGAL::insert_non_intersecting_curves(env, p.begin(), p.end());
    // find the face of the query point
    // (usually you may know that by other means)
    Arrangement_2::Face_const_handle *face;
    CGAL::Arr_naive_point_location<Arrangement_2> pl(env);
    CGAL::Arr_point_location_result<Arrangement_2>::Type obj = pl.locate(point_sample);
    // The query point locates in the interior of a face
    face = boost::get<Arrangement_2::Face_const_handle>(&obj);

    // compute non regularized visibility area
    // Define visibility object type that computes non-regularized visibility area
    Arrangement_2 non_regular_output;
    NSPV non_regular_visibility(env);
    non_regular_visibility.compute_visibility(point_sample, *face, non_regular_output);
    std::cout << "Non-regularized visibility region of q has "
              << non_regular_output.number_of_edges()
              << " edges:" << std::endl;
    for (Edge_const_iterator eit = non_regular_output.edges_begin(); eit != non_regular_output.edges_end(); ++eit)
        std::cout << "[" << eit->source()->point() << " -> " << eit->target()->point() << "]" << std::endl;

    CGAL::Polygon_2<Kernel> vis;


    for(auto var : non_regular_output.vertex_handles()){
        vis.push_back(var->point());
    }


    auto os = vis.oriented_side(point_sample2);

    cout << vis << endl;
    for (auto it = vis.edges_begin(); it != vis.edges_end(); ++it){
        cout << it->point(0) << "->" << it->point(1) << endl;
    }
    cout << point_sample << " : " << point_sample2 << endl;
    cout << nice(os) << endl;


    return 0;
}
