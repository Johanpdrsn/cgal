//
// Created by Johan Pedersen on 23/03/2022.
//

#include <CGAL/Convexity_measure_2.h>
#include <CGAL/Polygon_2.h>
#include <unistd.h>
#include <climits>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polygon_2<Kernel> Polygon_2;
typedef Kernel::Point_2 Point_2;
typedef CGAL::Arr_segment_traits_2<Kernel> Traits_2;
typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;

void read_polygon(const std::string &fileName, Polygon_2 &polygon) {
    std::ifstream in{fileName};

    if (!in.is_open()) {
        throw std::runtime_error("Could not open polygon file: " + fileName);
    }

    Kernel::FT x{};
    Kernel::FT y{};

    int n;
    in >> n;

    while (in >> x) {
        in >> y;
        polygon.push_back(Point_2{x, y});
    }
    in.close();
}

int main() {

    const int N = 100000;
    char dir[PATH_MAX];
    getcwd(dir, sizeof(dir));
    const std::string fileName{std::string(dir) + "/../100.line"};

    Polygon_2 polygon;
    polygon.push_back(Point_2{0.0, 0.0});
    polygon.push_back(Point_2{-2.0, -1.0});
    polygon.push_back(Point_2{1.0, -1.0});
    polygon.push_back(Point_2{1.0, 2.0});

    std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double, std::milli>> t1;
    std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double, std::milli>> t2;
    std::chrono::duration<double, std::milli> ms_double{};

    CGAL::Convexity_measure_2<Kernel> conv_pol{polygon};

    t1 = std::chrono::high_resolution_clock::now();
    Kernel::FT two_points = conv_pol.two_point_visibility_sample(N);
    t2 = std::chrono::high_resolution_clock::now();
    ms_double = t2 - t1;
    std::cout << "Two points prob for polygon: " << two_points << " in " << ms_double.count() << "ms" << std::endl;

    t1 = std::chrono::high_resolution_clock::now();
    Kernel::FT pol_vis = conv_pol.visibility_polygon_sample(N);
    t2 = std::chrono::high_resolution_clock::now();
    ms_double = t2 - t1;
    std::cout << "Polygon prob for polygon: " << pol_vis << " in " << ms_double.count() << "ms" << std::endl;


    Polygon_2 file_polygon;
    read_polygon(fileName, file_polygon);
    CGAL::Convexity_measure_2<Kernel> conv_path{file_polygon};
    const int M = N / 10;

    t1 = std::chrono::high_resolution_clock::now();
    two_points = conv_path.two_point_visibility_sample(M);
    t2 = std::chrono::high_resolution_clock::now();
    ms_double = t2 - t1;
    std::cout << "Two points prob for file: " << two_points << " in " << ms_double.count() << "ms" << std::endl;

    t1 = std::chrono::high_resolution_clock::now();
    pol_vis = conv_path.visibility_polygon_sample(M);
    t2 = std::chrono::high_resolution_clock::now();
    ms_double = t2 - t1;
    std::cout << "Polygon prob for file: " << pol_vis << " in " << ms_double.count() << "ms" << std::endl;


    return EXIT_SUCCESS;
}
