//
// Created by Johan Pedersen on 23/03/2022.
//

#ifndef VISIBILITY_2_EXAMPLES_CONVEXITY_MEASURE_2_H
#define VISIBILITY_2_EXAMPLES_CONVEXITY_MEASURE_2_H

#include <CGAL/license/Visibility_2.h>


#include <CGAL/Polygon_2.h>
#include <CGAL/Triangular_expansion_visibility_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

namespace CGAL {

    const int RECURSION_DEPTH = 10000;

    struct FaceInfo2 {
        FaceInfo2() = default;

        int nesting_level;

        bool in_domain() const {
            return nesting_level % 2 == 1;
        }
    };

    template<class Arrangement_2_, class RegularizationCategory = CGAL::Tag_true>
    class Convexity_measure_2 final {
    public:

        typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
        typedef CGAL::Polygon_2<Kernel> Polygon_2;
        typedef Kernel::Triangle_2 Triangle_2;
        typedef Kernel::Point_2 Point_2;
        typedef Kernel::Segment_2 Segment_2;

        explicit Convexity_measure_2(Polygon_2 &input_polygon) : polygon(input_polygon) {
            triangulate();
        }

        Kernel::FT two_point_visibility_sample(const int n) const {
            assert(n > 0);

            CGAL::Random_points_in_triangles_2<Point_2> point_generator{triangles};

            int sum = 0;
            Segment_2 seg;
            for (int i = 0; i < n; ++i) {
                seg = Segment_2{*point_generator++, *point_generator++};

                for (auto it = polygon.edges_begin(); it != polygon.edges_end(); ++it) {
                    if (CGAL::do_intersect(seg, it.make_value_type(CGAL::Tag_true()))) {
                        ++sum;
                        break;
                    }
                }
            }
            return (1.0 - ((sum / static_cast<double>(n))));
        }

        Kernel::FT visibility_polygon_sample(const int n) const {
            typedef Arrangement_2_ Arrangement_2;
            typedef CGAL::Arr_segment_traits_2<Kernel> Traits_2;
            typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;
            typedef CGAL::Triangular_expansion_visibility_2<Arrangement_2, RegularizationCategory> TEV;

            assert(n > 0);

            CGAL::Random_points_in_triangles_2<Point_2> point_generator{triangles};

            Arrangement_2 env;
            CGAL::insert(env, polygon.edges_begin(), polygon.edges_end());
            TEV regular_visibility{env};

            Arrangement_2::Face_const_handle *face;
            CGAL::Arr_naive_point_location<Arrangement_2> pl{env};
            CGAL::Arr_point_location_result<Arrangement_2>::Type obj;
            Arrangement_2 regular_output;

            Point_2 point_sample;

            Kernel::FT area = 0;
            Kernel::FT sum = 0;

            for (int i = 0; i < n; ++i) {
                // find the face of the query point
                point_sample = *point_generator++;
                obj = pl.locate(point_sample);
                // The query point locates in the interior of a face
                face = boost::get<Arrangement_2::Face_const_handle>(&obj);

                // compute regularized visibility_sample area
                regular_visibility.compute_visibility(point_sample, *face, regular_output);

                // Sadly the point type in the visibility result can't be used for area computation
                std::vector<Point_2> visible_points;
                std::transform(regular_output.vertices_begin(), regular_output.vertices_end(),
                               std::back_inserter(visible_points),
                               [](const auto &x) -> Point_2 { return x.point(); });

                //  Sum the area
                sum += CGAL::polygon_area_2(visible_points.begin(), visible_points.end(), Kernel());

                //Destructor for Lazy_nt is recursive, so we have to resolve at some interval smaller than ~25000
                if (i % RECURSION_DEPTH == 0) {
                    sum = sum.exact();
                }
            }
            // Normalize to the size of the polygon
            return sum / (polygon.area() * n);
        }

    private:
        typedef CGAL::Triangulation_vertex_base_2<Kernel> Vb;
        typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, Kernel> Fbb;
        typedef CGAL::Constrained_triangulation_face_base_2<Kernel, Fbb> Fb;
        typedef CGAL::Triangulation_data_structure_2<Vb, Fb> TDS;
        typedef CGAL::Exact_intersections_tag Itag;
        typedef CGAL::Constrained_triangulation_2<Kernel, TDS, Itag> CT;

        Polygon_2 polygon;
        std::vector<Triangle_2> triangles;

        static void mark_domains(CT &ct,
                                 CT::Face_handle start,
                                 int index,
                                 std::list<CT::Edge> &border) {
            if (start->info().nesting_level != -1) {
                return;
            }
            std::list<CT::Face_handle> queue;
            queue.push_back(start);
            while (!queue.empty()) {
                CT::Face_handle fh = queue.front();
                queue.pop_front();
                if (fh->info().nesting_level == -1) {
                    fh->info().nesting_level = index;
                    for (int i = 0; i < 3; ++i) {
                        CT::Edge e(fh, i);
                        CT::Face_handle n = fh->neighbor(i);
                        if (n->info().nesting_level == -1) {
                            if (ct.is_constrained(e)) border.push_back(e);
                            else queue.push_back(n);
                        }
                    }
                }
            }
        }

//explore set of facets connected with non-constrained edges,
//and attribute to each such set a nesting level.
//We start from facets incident to the infinite vertex, with a nesting
//level of 0. Then we recursively consider the non-explored facets incident
//to constrained edges bounding the former set and increase the nesting level by 1.
//Facets in the domain are those with an odd nesting level.
        static void mark_domains(CT &cdt) {
            for (CT::Face_handle f: cdt.all_face_handles()) {
                f->info().nesting_level = -1;
            }
            std::list<CT::Edge> border;
            mark_domains(cdt, cdt.infinite_face(), 0, border);
            while (!border.empty()) {
                CT::Edge e = border.front();
                border.pop_front();
                CT::Face_handle n = e.first->neighbor(e.second);
                if (n->info().nesting_level == -1) {
                    mark_domains(cdt, n, e.first->info().nesting_level + 1, border);
                }
            }
        }

        void triangulate() {
            CT ct;
            ct.insert_constraint(polygon.vertices_begin(), polygon.vertices_end(), true);
            mark_domains(ct);

            for (auto var: ct.finite_face_handles()) {
                if (var->info().in_domain()) {
                    triangles.emplace_back(ct.triangle(var));
                }
            }
        }
    };
}

#endif //VISIBILITY_2_EXAMPLES_CONVEXITY_MEASURE_2_H
