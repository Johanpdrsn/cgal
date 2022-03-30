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

void create_zigzag(Polygon_2 &polygon, int n) {

    std::deque<int> front(std::floor(static_cast<double>(n) / 2.0));
    std::iota(front.begin(), front.end(), std::ceil(static_cast<double>(n) / 2.0));


    for (auto a: front) {
        std::cout << a << std::endl;
    }

    for (int i = 0; !front.empty(); i++) {
        if (i % 2 == 0) {
            polygon.push_back(Point_2{i, 1});
        } else {
            polygon.push_back(Point_2{i, front.back()});
            front.pop_back();
        }
    }
    if (n % 2 != 0) polygon.push_back(Point_2{n - 1, 1});

    std::deque<int> back(std::floor(static_cast<double>(n) / 2.0));
    std::iota(back.begin(), back.end(), std::ceil(static_cast<double>(n) / 2.0) - 1);
    for (auto a: back) {
        std::cout << a << std::endl;
    }

    for (int i = n - 1; !back.empty(); i--) {

        if (i % 2 == 0) {
            polygon.push_back(Point_2{i, 0});
        } else {
            polygon.push_back(Point_2{i, back.front()});
            back.pop_front();
        }
    }
    polygon.push_back(Point_2{0, 0});

}


int main() {
    CGAL::IO::set_pretty_mode(std::cout);

    std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double, std::milli>> t1;
    std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double, std::milli>> t2;
    std::chrono::duration<double, std::milli> ms_double{};
    Kernel::FT two_points;
    Kernel::FT pol_vis;

    const int N = 100000;

    Polygon_2 polygon;
    create_zigzag(polygon, 20);

    
    CGAL::Convexity_measure_2<Kernel> conv_pol{polygon};

    t1 = std::chrono::high_resolution_clock::now();
    two_points = conv_pol.two_point_visibility_sample(N);
    t2 = std::chrono::high_resolution_clock::now();
    ms_double = t2 - t1;
    std::cout << "Two points prob for polygon: " << two_points << " in " << ms_double.count() << "ms" << std::endl;

    t1 = std::chrono::high_resolution_clock::now();
    pol_vis = conv_pol.visibility_polygon_sample(N);
    t2 = std::chrono::high_resolution_clock::now();
    ms_double = t2 - t1;
    std::cout << "Polygon prob for polygon: " << pol_vis << " in " << ms_double.count() << "ms" << std::endl;

//    char dir[PATH_MAX];
//    getcwd(dir, sizeof(dir));
//    const std::string fileName{std::string(dir) + "/../100.line"};
//
//    Polygon_2 file_polygon;
//    read_polygon(fileName, file_polygon);
//    CGAL::Convexity_measure_2<Kernel> conv_path{file_polygon};
//    const int M = 10000;
//
//    t1 = std::chrono::high_resolution_clock::now();
//    two_points = conv_path.two_point_visibility_sample(M);
//    t2 = std::chrono::high_resolution_clock::now();
//    ms_double = t2 - t1;
//    std::cout << "Two points prob for file: " << two_points << " in " << ms_double.count() << "ms" << std::endl;
//
//    t1 = std::chrono::high_resolution_clock::now();
//    pol_vis = conv_path.visibility_polygon_sample(M);
//    t2 = std::chrono::high_resolution_clock::now();
//    ms_double = t2 - t1;
//    std::cout << "Polygon prob for file: " << pol_vis << " in " << ms_double.count() << "ms" << std::endl;


    return EXIT_SUCCESS;
}
