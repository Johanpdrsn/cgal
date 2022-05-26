//
// Created by Johan Pedersen on 25/04/2022.
//

#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Exact_convexity_measure_2.h>


typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point_2;
typedef CGAL::Polygon_2<Kernel> Polygon_2;

using namespace std;


int main() {
    CGAL::IO::set_pretty_mode(std::cout);

    Polygon_2 polygon;

    polygon.push_back(Point_2{0, 0});
    polygon.push_back(Point_2{1, -2.0});
    polygon.push_back(Point_2{2, 0});
    polygon.push_back(Point_2{2, 2});




    auto test = Convexity_measure_exact_2<Kernel>(polygon);

    cout << "Result: " << test.measure_convexity() << endl;

//    test.tree.prettyPrint();
//

    return 0;
}