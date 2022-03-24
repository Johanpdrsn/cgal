//
// Created by Johan Pedersen on 23/03/2022.
//

#include <CGAL/Convexity_measure_2.h>
#include <CGAL/Polygon_2.h>

int main() {
    typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
    typedef CGAL::Polygon_2<Kernel> Polygon_2;
    typedef Kernel::Point_2 Point_2;

    const int N = 1000000;
    const std::string home = getenv("HOME");
    const std::string fileName{home+"/Documents/Thesis/cgal/Visibility_2/examples/Visibility_2/100.line"};

    Polygon_2 polygon;
    polygon.push_back(Point_2{0.0, 0.0});
    polygon.push_back(Point_2{-2.0, -1.0});
    polygon.push_back(Point_2{1.0, -1.0});
    polygon.push_back(Point_2{1.0, 2.0});

    CGAL::Convexity_measure_2 conv_pol{polygon};


    std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double, std::milli>> t1;
    std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double, std::milli>> t2;
    std::chrono::duration<double, std::milli> ms_double{};

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

        CGAL::Convexity_measure_2 conv_path{fileName};
        const double M = 1000;

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
