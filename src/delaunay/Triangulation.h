/****************************************************************************
 *                                                                          *
 *  Author : lukasz.iwaszkiewicz@gmail.com                                  *
 *  ~~~~~~~~                                                                *
 *  License : see COPYING file for details.                                 *
 *  ~~~~~~~~~                                                               *
 ****************************************************************************/

#ifndef DELAUNAYTRIANGULATION_H_
#define DELAUNAYTRIANGULATION_H_

#include <boost/polygon/voronoi.hpp>
#include <boost/math/special_functions/round.hpp>
#include <vector>
#include "Index.h"
#include "Triangle.h"
#include "Point.h"
#include "InputCollection.h"
#ifndef NDEBUG
#include <boost/timer/timer.hpp>
#endif

namespace Delaunay {

/**
 *
 */
template<typename T>
struct triangulation_voronoi_diagram_traits {
        typedef T coordinate_type;
        typedef boost::polygon::voronoi_cell<coordinate_type> cell_type;
        typedef boost::polygon::voronoi_vertex<coordinate_type> vertex_type;
        typedef boost::polygon::voronoi_edge<coordinate_type> edge_type;

        // This copes with degenerations.
        typedef struct {
                bool operator() (const vertex_type& v1, const vertex_type& v2) const
                {
                        return false;
                }
        } vertex_equality_predicate_type;
};

/**
 *
 */
typedef boost::polygon::voronoi_diagram<double, triangulation_voronoi_diagram_traits <double> > triangulation_voronoi_diagram;

/**
 *
 */
template <
        typename PointArg = Point,
        typename TriangleArg = Triangle,
        typename PointList = std::vector <PointArg>
>
class Triangulation {
public:

        typedef InputCollection <PointArg, TriangleArg, PointList> InputCollectionType;
        typedef Index <PointArg, TriangleArg, PointList> DelaunayIndexType;
        typedef TypeTraits <PointArg, TriangleArg, PointList> Traits;
        typedef typename Traits::PointType PointType;
        typedef typename Traits::TriangleType TriangleType;
        typedef typename Traits::PointListType PointListType;
        typedef typename Traits::PointListCollectionType PointListCollectionType;
        typedef typename Traits::PointTraitsType PointTraitsType;
        typedef typename Traits::CoordinateType CoordinateType;
        typedef typename Traits::EdgeType EdgeType;
        typedef typename Traits::TriangleTraitsType TriangleTraitsType;
        typedef typename Traits::IndexType IndexType;
        typedef typename Traits::TriangleEdgeType TriangleEdgeType;
        typedef typename Traits::TriangleEdgeList TriangleEdgeList;
        typedef typename Traits::TriangleEdgeVector TriangleEdgeVector;
        typedef typename Traits::TriangleVector TriangleVector;
        typedef typename Traits::TrianglePtrVector TrianglePtrVector;
        typedef typename Traits::TriangleIndex TriangleIndex;
        typedef typename Traits::IntersectionInfo IntersectionInfo;
        typedef triangulation_voronoi_diagram::vertex_type vertex_type;
        typedef triangulation_voronoi_diagram::edge_type edge_type;
        typedef triangulation_voronoi_diagram::cell_type cell_type;

        enum LineMode { LINES, LINE_LOOP };

        Triangulation () : index (input) {}

        void constructDelaunay (/*Geometry::LineString *crossing*/);

        void makeVoronoiDual ();
        void linkTriangles ();
        void findMissingConstraints (InputCollectionType const &allConstraints, TriangleEdgeList *missingConstraints, LineMode lineMode = LINE_LOOP);

        TriangleVector const &getTriangulation () const { return index.getTriangulation (); }

        void setPoints (PointListType const &p) { input.setPoints (p); }
        PointListType const &getPoints () const { return input.getPoints (); }

        void addConstraint (PointListType const &p) { input.addConstraint (p); }
        PointListCollectionType const &getConstraints () const { return input.getConstraints (); }

        size_t getDataSize () const { return input.size (); }

private:

        /**
         * Returns if diagonal (a, b) lays completely inside the polygon. It assumes two things:
         * - Points in polygon are stored in counter clockwise order.
         * - Checked diagonal do not intersects with boundary segments (defined by polygon itself)
         * in points other than a and b.
         * This second condition is easy to conform to, since we operate on diagonals generated
         * by triangulation algotihm.
         * \param a diagonal point a
         * \param ap point previous to a (in input).
         * \param an point next to a
         * \param b diagonal point b
         */
        bool diagonalInside (PointType const& a,
                             PointType const& ap,
                             PointType const& an,
                             PointType const& b) const;

        bool diagonalInside (TriangleEdgeType const &e) const;
        bool triangleInside (TriangleType const &t) const;

private:

        InputCollectionType input;
        DelaunayIndexType index;
        triangulation_voronoi_diagram vd;
};

/****************************************************************************/

template <
        typename PointArg,
        typename TriangleArg,
        typename PointList
>
void Triangulation<PointArg, TriangleArg, PointList>::constructDelaunay (/*Geometry::LineString *crossing*/)
{
        index.reserve ();
        makeVoronoiDual ();
        linkTriangles ();

        TriangleEdgeList missingConstraints;
        findMissingConstraints (input, &missingConstraints);

#if 1
        // 4. Add missing segments.
        for (typename TriangleEdgeList::const_iterator i = missingConstraints.begin (); i != missingConstraints.end (); ++i) {
                TriangleEdgeType const &missingConstraint = *i;

#if 0
                if (missingConstraint.a != 460) {
//                if (missingConstraint.a != 119) {
                        continue;
                }
                std::cerr << "Missing constraint : " << missingConstraint << std::endl;
#endif

                TriangleEdgeList crossingEdges;
                TrianglePtrVector crossingTriangles;

                // Paragraph 2.
                index.findCrossingEdges (missingConstraint, &crossingEdges, &crossingTriangles);

#if 0
                std::cerr << "Constraint " << missingConstraint << " crosses : " << crossingTriangles << std::endl;
#endif

#if 0
                for (typename TrianglePtrVector::const_iterator i = crossingTriangles.begin (); i != crossingTriangles.end (); ++i) {
                        Delaunay::Point const &a = input[Delaunay::a (**i)];
                        Delaunay::Point const &b = input[Delaunay::b (**i)];
                        Delaunay::Point const &c = input[Delaunay::c (**i)];

                        crossing->push_back (Geometry::makePoint (a.x, a.y));
                        crossing->push_back (Geometry::makePoint (b.x, b.y));
                        crossing->push_back (Geometry::makePoint (c.x, c.y));
                }
#endif

                // Paragraph 3.
                TriangleEdgeList newEdges;
                while (!crossingEdges.empty ()) {

                        // Paragraph 3.1
                        TriangleEdgeType currentCrossingEdge = crossingEdges.front ();
                        crossingEdges.pop_front ();

                        // Paragraph 3.2
                        if (currentCrossingEdge.covers (TriangleEdgeType (549, 8)))
                        if (!index.twoTrianglesConvex (currentCrossingEdge)) {
                                crossingEdges.push_back (currentCrossingEdge);

#if 0
                                std::cerr << "####> !CONVEX" << std::endl;
#endif

                                if (crossingEdges.size () == 1) {
                                        /*
                                         * TODO Jeśli jest tylko jedna przecinająca dany constraint, to jeśli convex,
                                         * to flip, a jeśli nie, to nie wiem, ale coś trzeba tu zrobić.
                                         */
                                        assert (0); // not implemented TODO - swoją drogą - to jest chyab niemożliwe (jeśli program dobrze działa)?!
                                }

                                continue;
                        }

#if 0
                        std::cerr << "####> +++CONVEX" << std::endl;
#endif
                        // Dwa przyległę trójkąty zawierające e tworzą czworobok wypukły.
                        TriangleEdgeType newDiagonal;

                        /*
                         * Tu trzeba
                         * - uaktualnić wierzchołki trójkątów.
                         * - uaktualnić ich zlinkowane trójkąty.
                         * - uaktualnic triangleIndex (potrzebny w findCrossing edges).
                         */
                        index.flip (currentCrossingEdge, &newDiagonal);

                        EdgeType a = index.triangleEdgeToEdge (missingConstraint);
                        EdgeType b = index.triangleEdgeToEdge (newDiagonal);

                        if (Delaunay::intersects (a, b)) {
                                crossingEdges.push_back (newDiagonal);
                        }
                        else {
                                newEdges.push_back (newDiagonal);
                        }
                }

                // 4. Make CDT from DT.
#if 0
                std::cerr << newEdges.size () << std::endl;
#endif
                for (typename TriangleEdgeList::iterator i = newEdges.begin (); i != newEdges.end (); ++i) {
                        TriangleEdgeType &newEdge = *i;

                        if (newEdge.covers (missingConstraint)) {
                                continue;
                        }

                        if (!index.twoTrianglesNotDelaunay (newEdge)) {
                                TriangleEdgeType newDiagonal;
                                index.flip (newEdge, &newDiagonal);

                                // TODO w piśmie napsali, żeby ją zachowac w kolekcji nowych krawędzi - ale po co?
                                newEdge = newDiagonal;
                        }
                }
        }

        // TODO z tego zrobić metodę i dać parametr czy CW czy CCW
        // 5. Remove superfluous triangles.  Input must be in COUNTER CLOCKWISE order.
        TriangleVector const &triangulation = index.getTriangulation ();

        for (typename TriangleVector::const_iterator i = triangulation.begin (), e = triangulation.end (); i != e; ++i) {
                TriangleType &triangle = const_cast <TriangleType &> (*i);

                if (a (triangle) == 0 && b (triangle) == 0 && c (triangle) == 0) {
                        continue;
                }

                if (!triangleInside (triangle)) {
                        a (triangle, 0);
                        b (triangle, 0);
                        c (triangle, 0);
                }
        }

        index.clean ();

#if 0
        printlog ("Triangulation time (derived from voronoi as its dual) : %f ms", t1.elapsed ().wall / 1000000.0);
        std::cerr << "CDT size : " << index.getNumTriangles () << " triangles." << std::endl;
//        std::cout << triangulation << std::endl;
#endif
#endif
}

/****************************************************************************/

template <
        typename PointArg,
        typename TriangleArg,
        typename PointList
>
void Triangulation <PointArg, TriangleArg, PointList>::makeVoronoiDual ()
{
#ifndef NDEBUG
        boost::timer::cpu_timer t0;
#endif

        // construct_voronoi (input.begin (), input.end (), &vd);
        boost::polygon::default_voronoi_builder builder;

        PointListCollectionType const &pointsCollection = input.getPointsCollection ();
        for (typename PointListCollectionType::const_iterator i = pointsCollection.begin (), e = pointsCollection.end (); i != e; ++i) {
                PointListType const *pointList = *i;
                boost::polygon::insert (pointList->begin (), pointList->end (), &builder);
        }

        builder.construct (&vd);

#ifndef NDEBUG
        std::cerr << "Voronoi diagram construction time : " << t0.elapsed ().wall / 1000000.0 << " ms" << std::endl;
        boost::timer::cpu_timer t1;
#endif

        // 1. Make triangles from voronoi.
        for (triangulation_voronoi_diagram::const_vertex_iterator it = vd.vertices ().begin (); it != vd.vertices ().end (); ++it) {
                vertex_type const &vertex = *it;
                edge_type const *edge = vertex.incident_edge ();

                IndexType triangleVertices[3];
                IndexType triangleCnt = index.getNumTriangles ();
                vertex.color (triangleCnt);

                size_t edgeCnt = 0;
                do {
                        if (!edge->is_primary ()) {
                                continue;
                        }

                        cell_type const *cell = edge->cell ();
                        size_t index = cell->source_index ();
                        triangleVertices[edgeCnt++] = index;
                        // Add 1 so that default 0 is invalid.
                        edge->color (triangleCnt + 1);
                        edge = edge->rot_next ();
                } while (edge != vertex.incident_edge () && edgeCnt < 3);

                if (edgeCnt == 3) {
                        TriangleType triangle;
                        a (triangle, triangleVertices[0]);
                        b (triangle, triangleVertices[1]);
                        c (triangle, triangleVertices[2]);
                        index.addTriangle (triangle);
                }
        }

        index.sortEdgeIndex ();

#if 0
        std::cerr << "Delaunay triangulation produced : " << index.getNumTriangles () << " triangles." << std::endl;
        std::cerr << triangulation << std::endl;
#endif
}

/****************************************************************************/

template <
        typename PointArg,
        typename TriangleArg,
        typename PointList
>
void Triangulation <PointArg, TriangleArg, PointList>::linkTriangles ()
{
        // 2. Link triangles.
        for (triangulation_voronoi_diagram::const_vertex_iterator it = vd.vertices ().begin (); it != vd.vertices ().end (); ++it) {
                vertex_type const &vertex = *it;
                TriangleType &triangle = index.getTriangle (vertex.color ());
                edge_type const *edge = vertex.incident_edge ();
                size_t edgeCnt = 0;

                do {
                        if (!edge->is_primary ()) {
                                continue;
                        }

                        uint32_t adjacentIndex = edge->twin ()->color ();

                        // Index with value 0 is invalid.
                        if (!adjacentIndex) {
                                edge = edge->rot_next ();
                                continue;
                        }

                        TriangleType &adjacentTriangle = index.getTriangle (adjacentIndex - 1);

                        // Sprawdzić którym bokiem się stykają i ustawić.
                        if (!hasVertex (adjacentTriangle, triangle.a)) {
                                index.setAdjacentTriangle(triangle, A, &adjacentTriangle);
                        }
                        else if (!hasVertex (adjacentTriangle, triangle.b)) {
                                index.setAdjacentTriangle(triangle, B, &adjacentTriangle);
                        }
                        else if (!hasVertex (adjacentTriangle, triangle.c)) {
                                index.setAdjacentTriangle(triangle, C, &adjacentTriangle);
                        }

                        ++edgeCnt;
                        edge = edge->rot_next ();
                } while (edge != vertex.incident_edge () && edgeCnt < 3);
        }

#if 0
        std::cout << triangulation << std::endl;
#endif
}

/****************************************************************************/

template <
        typename PointArg,
        typename TriangleArg,
        typename PointList
>
void Triangulation <PointArg, TriangleArg, PointList>::findMissingConstraints (InputCollectionType const &allConstraints, TriangleEdgeList *missingConstraints, LineMode lineMode)
{
        // 3. Find missing constraints. Update triangleVector (data structure for CDT).
        /*
         * TODO This is loop made for simple polygons (without holes). It is also possible to make
         * loop for discrete list of constraints (that are not linked).
         *
         * W tym kawałku chodzi o to, żeby znaleźć wszystkie constrainty. Szukamy w triangulacji
         * krawędzi o wierzchołkach [i, i+1]. Jeśli jakiegoś nie ma, to dodajemy go do listy
         * brakujących constraintów.
         */
        size_t inputSize = allConstraints.size ();

//        for (typename InputCollectionType::const_iterator it = allConstraints.begin (), e = allConstraints.end (); it != e; ++it) {
        for (size_t i = 0; i < inputSize; ++i) {
                size_t j = (i + 1) % inputSize;

                assert (index.getTriangleIndexSize () > i);
                TrianglePtrVector const &trianglesForPoint = index.getTrianglesForPoint (i);

#ifndef NDEBUG
                if (trianglesForPoint.empty ()) {
                        std::cerr << "UWAGA : trianglesForPoint.empty () : nie ma trójkątów stycznych do punktu o indeksie : " << i << ". TODO rozkminić czy to OK. i_max = " << inputSize -1  << std::endl;
                }
#endif

                assert (!trianglesForPoint.empty ());

                bool found = false;
                for (typename TrianglePtrVector::const_iterator k = trianglesForPoint.begin (); k != trianglesForPoint.end (); ++k) {
                        if (hasEdge (**k, TriangleEdgeType (i, j))) {
                                found = true;
                                break;
                        }
                }

                if (!found) {
#if 0
                        std::cerr << "Constraint (" << i << ", " << j << ") was **NOT** found in triangulation." << std::endl;
#endif
                        missingConstraints->push_back (TriangleEdgeType (i, j));
                }
        }
}

/****************************************************************************/

template <
        typename PointArg,
        typename TriangleArg,
        typename PointList
>
bool Triangulation <PointArg, TriangleArg, PointList>::diagonalInside (PointType const &a, PointType const &ap, PointType const &an, PointType const &b) const
{
        int apx = ap.x - a.x;
        int apy = ap.y - a.y;
        int anx = an.x - a.x;
        int any = an.y - a.y;
        int bx = b.x - a.x;
        int by = b.y - a.y;

        int apXan = apx * any - apy * anx;
        int apXb = apx * by - apy * bx;
        int bXan = bx * any - by * anx;

        return ((apXan >= 0 && apXb >= 0 && bXan >= 0) || (apXan < 0 && !(apXb < 0 && bXan < 0)));
}

/****************************************************************************/

template <
        typename PointArg,
        typename TriangleArg,
        typename PointList
>
bool Triangulation <PointArg, TriangleArg, PointList>::diagonalInside (TriangleEdgeType const &e) const
{
        // If e is one of constraints, we don't need to perform furher computations.
        if (std::abs (e.a - e.b) == 1) {
                return true;
        }

        size_t pointsSize = input.size ();
        PointType const &a = input[e.a];
        PointType const &ap = input[(e.a == 0) ? (pointsSize - 1) : (e.a - 1)];
        PointType const &an = input[(e.a == pointsSize - 1) ? (0) : (e.a + 1)];
        PointType const &b = input[e.b];

        return diagonalInside (a, ap, an, b);
}

/****************************************************************************/

template <
        typename PointArg,
        typename TriangleArg,
        typename PointList
>
bool Triangulation <PointArg, TriangleArg, PointList>::triangleInside (TriangleType const &t) const
{
        for (int i = 1; i <= 3; ++i) {
                if (!diagonalInside (getEdge (t, static_cast <SideEnum> (i)))) {
                        return false;
                }
        }

        return true;
}

} // namespace

#endif /* DELAUNAYTRIANGULATION_H_ */