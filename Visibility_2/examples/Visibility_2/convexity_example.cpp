//
// Created by Johan Pedersen on 23/03/2022.
//

#include <CGAL/Convexity_measure_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Cartesian.h>
#include <unistd.h>
#include <climits>
#include <iostream>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
//typedef CGAL::Cartesian<CGAL::Gmpq> Kernel;
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
    n+=2;
    polygon.push_back(Point_2{0,0});

    int i=1,k = 1,j=1;
    for (;i<n;i++,j+=i){
        if(i%2!=0){
            polygon.push_back(Point_2{j, k});

        } else{
            polygon.push_back(Point_2{j, -k});
            k++;
        }
    }

    j-=i;
    i--;
    auto t = k;
    if (n%2!=0)
        k-=2;
    else{
        k-=2;
        t++;
    }
    for (;i>0;j-=i,i--){
        if(i%2!=0){
            polygon.push_back(Point_2{j, t});
            t--;
        } else{
            polygon.push_back(Point_2{j, -k});
            k--;
        }
    }
    polygon.push_back(Point_2{0,1});
}


int main() {
    CGAL::IO::set_pretty_mode(std::cout);

    std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double, std::milli>> t1;
    std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double, std::milli>> t2;
    std::chrono::duration<double, std::milli> ms_double{};
    Kernel::FT two_points;
    Kernel::FT pol_vis;


    Polygon_2 polygon;
//    polygon.push_back(Point_2{0.0, 0.0});
//    polygon.push_back(Point_2{-2.0, -1.0});
//    polygon.push_back(Point_2{1.0, -1.0});
//    polygon.push_back(Point_2{1.0, 2.0});


    create_zigzag(polygon,10000);
//    std::cout << polygon << std::endl;
    std::cout << "CREATED POLYGON" << std::endl;
    CGAL::Convexity_measure_2<Kernel> conv_pol{polygon};
    const long N = 1000000000000;

    std::cout << "STARTING" << std::endl;
    t1 = std::chrono::high_resolution_clock::now();
    for(int i=1;i<=200;i++){
        two_points = conv_pol.two_point_visibility_sample(N);
    }
    t2 = std::chrono::high_resolution_clock::now();
    ms_double = t2 - t1;
//    std::cout << "Two points prob for polygon: " << CGAL::to_double(two_points) << " in " << ms_double.count() << "ms" << std::endl;
//
//    t1 = std::chrono::high_resolution_clock::now();
//    pol_vis = conv_pol.visibility_polygon_sample(N);
//    t2 = std::chrono::high_resolution_clock::now();
//    ms_double = t2 - t1;
//    std::cout << "Polygon prob for polygon: " << CGAL::to_double(pol_vis) << " in " << ms_double.count() << "ms" << std::endl;

//    char dir[PATH_MAX];
//    getcwd(dir, sizeof(dir));
//    const std::string fileName{std::string(dir) + "/../50.line"};
//
//    Polygon_2 file_polygon;
//    read_polygon(fileName, file_polygon);
//    CGAL::Convexity_measure_2<Kernel> conv_path{file_polygon};
//    const int M = 100000;
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
