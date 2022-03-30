//
// Created by Johan Pedersen on 23/03/2022.
//

#ifndef CONVEXITY_MEASURE_2_VISIBILITY_2_H
#define CONVEXITY_MEASURE_2_VISIBILITY_2_H

#include <CGAL/license/Visibility_2.h>


#include <CGAL/Polygon_2.h>
#include <CGAL/Triangular_expansion_visibility_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

namespace CGAL {

    const int RECURSION_DEPTH = 10000;

    struct FaceInfo2 {
        FaceInfo2() = default;

        int nesting_level;

        bool in_domain() const {
            return nesting_level % 2 == 1;
        }
    };

    template<class Kernel>
    class Convexity_measure_2 final {
    public:

        typedef CGAL::Polygon_2<Kernel> Polygon_2;
        typedef typename Kernel::Point_2 Point_2;

        explicit Convexity_measure_2(Polygon_2 &input_polygon) : polygon(input_polygon) {
            triangulate();
        }

        typename Kernel::FT two_point_visibility_sample(const int n) const {
            typedef typename Kernel::Segment_2 Segment_2;

            assert(n > 0);

            CGAL::Random_points_in_triangles_2<Point_2> point_generator{triangles};

            int sum = 0;
            typename Kernel::Segment_2 seg;
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

        typename Kernel::FT visibility_polygon_sample(const int n) const {
            typedef CGAL::Arr_segment_traits_2<Kernel> Traits_2;
            typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;
            typedef CGAL::Triangular_expansion_visibility_2<Arrangement_2, CGAL::Tag_true> TEV;

            assert(n > 0);

            CGAL::Random_points_in_triangles_2<Point_2> point_generator{triangles};

            Arrangement_2 env;
            CGAL::insert(env, polygon.edges_begin(), polygon.edges_end());
            TEV regular_visibility{env};

            typename Arrangement_2::Face_const_handle *face;
            CGAL::Arr_naive_point_location<Arrangement_2> pl{env};
            typename CGAL::Arr_point_location_result<Arrangement_2>::Type obj;
            Arrangement_2 regular_output;

            Point_2 point_sample;

            typename Kernel::FT sum = 0;

            for (int i = 0; i < n; ++i) {
                // find the face of the query point
                point_sample = *point_generator++;
                obj = pl.locate(point_sample);
                // The query point locates in the interior of a face
                face = boost::get<typename Arrangement_2::Face_const_handle>(&obj);

                // compute regularized visibility_sample area
                regular_visibility.compute_visibility(point_sample, *face, regular_output);

                // Sadly the point type in the visibility result can't be used for area computation
                std::vector<Point_2> visible_points;
                std::transform(regular_output.vertices_begin(), regular_output.vertices_end(),
                               std::back_inserter(visible_points),
                               [](const auto &x) { return x.point(); });

                //  Sum the area
                sum += CGAL::polygon_area_2(visible_points.begin(), visible_points.end(), Kernel());


                //Destructor for Lazy_nt is recursive, so we have to resolve at some interval smaller than ~25000
                if (i % RECURSION_DEPTH == 0) {
                    sum = CGAL::exact(sum);
                }
            }
            // Normalize to the size of the polygon
            return sum / (CGAL::abs(polygon.area()) * n);
        }

    private:

        typedef CGAL::Triangulation_vertex_base_2<Kernel> Vb;
        typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, Kernel> Fbb;
        typedef CGAL::Constrained_triangulation_face_base_2<Kernel, Fbb> Fb;
        typedef CGAL::Triangulation_data_structure_2<Vb, Fb> TDS;
        typedef CGAL::Exact_intersections_tag Itag;
        typedef CGAL::Constrained_triangulation_2<Kernel, TDS, Itag> CT;
        typedef CGAL::Constrained_triangulation_plus_2<CT> CTP;

        typedef typename Kernel::Triangle_2 Triangle_2;

        Polygon_2 polygon;
        std::vector<Triangle_2> triangles;

        static void mark_domains(CTP &ctp,
                                 typename CTP::Face_handle start,
                                 int index,
                                 std::list<typename CTP::Edge> &border) {
            if (start->info().nesting_level != -1) {
                return;
            }
            std::list<typename CTP::Face_handle> queue;
            queue.push_back(start);
            while (!queue.empty()) {
                typename CTP::Face_handle fh = queue.front();
                queue.pop_front();
                if (fh->info().nesting_level == -1) {
                    fh->info().nesting_level = index;
                    for (int i = 0; i < 3; ++i) {
                        typename CTP::Edge e(fh, i);
                        typename CTP::Face_handle n = fh->neighbor(i);
                        if (n->info().nesting_level == -1) {
                            if (ctp.is_constrained(e)) border.push_back(e);
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
        static void mark_domains(CTP &ctp) {

            for (typename CTP::Face_handle f: ctp.all_face_handles()) {
                f->info().nesting_level = -1;
            }
            std::list<typename CTP::Edge> border;
            mark_domains(ctp, ctp.infinite_face(), 0, border);
            while (!border.empty()) {
                const typename CTP::Edge e = border.front();
                border.pop_front();
                typename CTP::Face_handle n = e.first->neighbor(e.second);
                if (n->info().nesting_level == -1) {
                    mark_domains(ctp, n, e.first->info().nesting_level + 1, border);
                }
            }
        }

        void triangulate() {
            CTP ctp;
            ctp.insert_constraint(polygon.vertices_begin(), polygon.vertices_end(), true);
            mark_domains(ctp);

            for (const auto &face: ctp.finite_face_handles()) {
                if (face->info().in_domain()) {
                    triangles.emplace_back(ctp.triangle(face));
                }
            }
        }
    };
}

#endif //CONVEXITY_MEASURE_2_VISIBILITY_2_H
