//
// Created by Johan Pedersen on 25/04/2022.
//

#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Exact_convexity_measure.h>


typedef CGAL::Cartesian<double> Kernel;
typedef Kernel::Point_2 Point_2;
typedef CGAL::Polygon_2<Kernel> Polygon_2;

using namespace std;


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

//    polygon.push_back(Point_2{0, 0});
//    polygon.push_back(Point_2{2, -2.0});
//    polygon.push_back(Point_2{2, 0});
//    polygon.push_back(Point_2{1, 2});

    polygon.push_back(Point_2{0.0, 0.0});
    polygon.push_back(Point_2{-2.0, -1.0});
    polygon.push_back(Point_2{1.0, -1.0});
    polygon.push_back(Point_2{1.0, 2.0});


    auto test = Convexity_measure_exact_2<Kernel>(polygon);

    cout << "Result: " << test.generate_tree() << endl;

//    test.tree.prettyPrint();
//

    return 0;
}