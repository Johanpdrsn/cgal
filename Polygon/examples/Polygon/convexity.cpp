#include <iostream>
#include <c++/v1/random>

#include <CGAL/partition_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_polygon_visibility_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_naive_point_location.h>

using namespace std;

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polygon_2<Kernel> Polygon_2;
typedef CGAL::Partition_traits_2<Kernel> Traits;
typedef Traits::Polygon_2 Polygon_T_2;
typedef std::list<Polygon_T_2> Polygon_list;


typedef Kernel::Point_2 Point_2;
typedef Kernel::Segment_2 Segment_2;
typedef CGAL::Arr_segment_traits_2<Kernel> Traits_2;
typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;

typedef CGAL::Simple_polygon_visibility_2<Arrangement_2, CGAL::Tag_false> NSPV;


unsigned long triangle_index_search(double sample, std::vector<double> &cdf) {
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

Point_2 kraemer_method(CGAL::Triangle_2<Kernel> &T) {

    random_device rand;
    mt19937 generator(rand());
    uniform_real_distribution<double> distribution(0, 1);

    double a = distribution(generator);
    double b = distribution(generator);

    double q = abs(a - b);

    double s = q;
    double t = 0.5 * (a + b - q);
    double u = 1 - 0.5 * (q + a + b);

    return Point_2{s * T[0].x() + t * T[1].x() + u * T[2].x(), s * T[0].y() + t * T[1].y() + u * T[2].y()};
}

Point_2 triSample(CGAL::Triangle_2<Kernel> &T) {

    random_device rand;
    mt19937 generator(rand());
    uniform_real_distribution<double> distribution(0, 1);

    double a = distribution(generator);
    double b = distribution(generator);

    double s1 = sqrt(a);

    return Point_2{T[0].x() * (1 - s1) + T[1].x() * (1 - b) * s1 + T[2].x() * b * s1,
                   T[0].y() * (1 - s1) + T[1].y() * (1 - b) * s1 + T[2].y() * b * s1};
}


std::vector<CGAL::Triangle_2<Kernel>> triangulate(Polygon_2 &poly) {
    if (!poly.is_counterclockwise_oriented()) {
        poly.reverse_orientation();
    }
    Polygon_list partition_polys;

    CGAL::approx_convex_partition_2(poly.vertices_begin(), poly.vertices_end(), back_inserter(partition_polys));

    std::vector<CGAL::Triangle_2<Kernel>> res;
    for (auto fac: partition_polys) {
        res.emplace_back(CGAL::Triangle_2<Kernel>(fac.vertex(0), fac.vertex(1), fac.vertex(2)));
    }

    return res;
}

std::vector<double> compute_prefix_sum(std::vector<CGAL::Triangle_2<Kernel>> &triangles) {
    std::vector<double> prefixSum{CGAL::to_double(triangles[0].area())};

    for (int i = 1; i < triangles.size(); i++)
        prefixSum[i] = prefixSum[i - 1] + CGAL::to_double(triangles[i].area());

    return prefixSum;
}


unsigned int two_point_visibility_sample(Polygon_2 &polygon, std::vector<CGAL::Triangle_2<Kernel>> &triangles,
                                         std::vector<double> &area_prefix_sum, mt19937 &generator,
                                         uniform_real_distribution<double> &distribution) {

    auto sample1 = distribution(generator);
    auto sample2 = distribution(generator);

    unsigned int triangle_index1 = triangle_index_search(sample1, area_prefix_sum);
    unsigned int triangle_index2 = triangle_index_search(sample2, area_prefix_sum);

    CGAL::Triangle_2<Kernel> triangle1 = triangles[triangle_index1];
    CGAL::Triangle_2<Kernel> triangle2 = triangles[triangle_index2];

    Point_2 point_sample1 = triSample(triangle1);
    Point_2 point_sample2 = triSample(triangle2);

    Segment_2 seg(point_sample1, point_sample2);

    for (auto it = polygon.edges_begin(); it != polygon.edges_end(); ++it) {
        if (CGAL::do_intersect(seg, Segment_2(it->point(0), it->point(1)))) {
            return 0;
        }
    }

    return 1;
}

double visibility_polygon_sample(Polygon_2 &polygon, std::vector<CGAL::Triangle_2<Kernel>> &triangles,
                                 std::vector<double> &area_prefix_sum, mt19937 &generator,
                                 uniform_real_distribution<double> &distribution) {


    auto sample1 = distribution(generator);

    unsigned int triangle_index = triangle_index_search(sample1, area_prefix_sum);

    CGAL::Triangle_2<Kernel> triangle1 = triangles[triangle_index];

    vector<Segment_2> p;

    for (auto it = polygon.edges_begin(); it != polygon.edges_end(); ++it) {
        p.emplace_back(Segment_2(it->point(0), it->point(1)));
    }

    Point_2 point_sample = kraemer_method(triangle1);

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

    CGAL::Polygon_2<Kernel> vis;
    for (auto var: non_regular_output.vertex_handles()) {
        vis.push_back(var->point());
    }

    return abs(CGAL::to_double(vis.area() / polygon.area()));
}


double simulate(int n, const function<double(Polygon_2 &, std::vector<CGAL::Triangle_2<Kernel>> &,
                                             std::vector<double> &, mt19937 &,
                                             uniform_real_distribution<double> &)> &vis_func) {
    Polygon_2 polygon;
    polygon.push_back(Point_2(0, 0));
    polygon.push_back(Point_2(-2, -1));
    polygon.push_back(Point_2(1, -1));
    polygon.push_back(Point_2(1, 2));

    std::vector<CGAL::Triangle_2<Kernel>> triangles = triangulate(polygon);
    std::vector<double> area_prefix_sum = compute_prefix_sum(triangles);

    random_device rand;
    mt19937 generator(rand());
    uniform_real_distribution<double> distribution(0, CGAL::to_double(polygon.area()));

    double sum = 0;

    for (int i = 0; i < n; i++) {
        sum += vis_func(polygon, triangles, area_prefix_sum, generator, distribution);
    }

    return sum / n;
}


int main() {
    int n = 10000;

    chrono::time_point<chrono::steady_clock, chrono::duration<double, milli>> t1;
    chrono::time_point<chrono::steady_clock, chrono::duration<double, milli>> t2;
    chrono::duration<double, std::milli> ms_double{};

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
