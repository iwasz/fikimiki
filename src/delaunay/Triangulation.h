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

        Triangulation () : index (input) {}

        void constructDelaunay (/*Geometry::LineString *crossing*/);

        void makeVoronoiDual ();
        void linkTriangles ();
        void findMissingConstraints (size_t constraintNumber, TriangleEdgeList *missingConstraints, bool closed = true);
        void addMissingSegments (TriangleEdgeList const &missingSegments);
//        void removeSuperfulous (PointListType const &pointList, bool removeOutside = true);
        void removeSuperfulous ();

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

        bool diagonalInside (TriangleEdgeType const &e, PointListType const &pointList) const;
        bool triangleInside (TriangleType const &t, PointListType const &pointList) const;
        void establishSuperTriangle (PointListType const &p);

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
        findMissingConstraints (0, &missingConstraints);
        addMissingSegments (missingConstraints);

        missingConstraints.clear ();
        findMissingConstraints (1, &missingConstraints);
        addMissingSegments (missingConstraints);

//        removeSuperfulous ();

        // TODO : Konstrukcja super-rectangle.
        // Usuwanie zewnętrznych.
        // Usuwanie wewnętrznych (dziur).

        // Numer pierwszego punktu super-rectangle
        IndexType first = input.getPoints().size ();
        // Pierwszy trójkat, który na pewno jest nazewnatrz.
        TriangleType *outside = index.getTrianglesForPoint(first).front ();
        TriangleType *inside = 0;

        /*
         * - Tworzyć na poczatku sper-rectangle.
         * -- Znaleźć bounding box zewnętrznego ringu.
         * -- Rozepchnąć ten bounding box o jakąś wartość i dodać punkty do wejściowych punktów.
         *
         * - Zaimplementować depth first search (DFS).
         * -- Dla początkowego trójkąta w triangulacji algorytm DFS ma biec dalej po trójkątach i dla kazdego z nich uruchomic jakąś funkcję.
         * -- W tej funkcji będziemy uzywać kolejnych funkcji, które będą klasyfikować trójkąt i oznaczac go jako do usunięcia lub nie.
         */


        while (true) {

                bool deadEnd = true;
                for (unsigned int i = 1; i <= 3; ++i) {
                        TriangleEdgeType e = getEdge (*outside, static_cast <SideEnum> (i));
                        TriangleType *t = index.getAdjacentTriangle(*outside, static_cast <SideEnum> (i));

                        std::cerr << "out = [" << a (*outside) << "," << b (*outside) << "," << c (*outside) << "], i = " << i << ", t=" << t << ", abse = " << std::abs (e.a - e.b) << std::endl;

                        if (!t) {
                                continue;
                        }

                        if (std::abs (e.a - e.b) == 1) {
                                inside = t;
                                continue;
                        }

                        IndexType ta = a (*t);
                        IndexType tb = b (*t);
                        IndexType tc = c (*t);

                        if (ta == 0 && tb == 0 && tc == 0) {
                                std::cerr << "t-zero" << std::endl;
                                continue;
                        }

                        a (*outside, 0);
                        b (*outside, 0);
                        c (*outside, 0);

                        outside = t;
                        deadEnd = false;
                        break;
                }

                if (deadEnd) {
                        break;
                }
        }

        index.clean ();

#if 0
        printlog ("Triangulation time (derived from voronoi as its dual) : %f ms", t1.elapsed ().wall / 1000000.0);
        std::cerr << "CDT size : " << index.getNumTriangles () << " triangles." << std::endl;
//        std::cout << triangulation << std::endl;
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
void Triangulation <PointArg, TriangleArg, PointList>::findMissingConstraints (size_t constraintNumber, TriangleEdgeList *missingConstraints, bool closed)
{
        // 3. Find missing constraints. Update triangleVector (data structure for CDT).
        /*
         * W tym kawałku chodzi o to, żeby znaleźć wszystkie constrainty. Szukamy w triangulacji
         * krawędzi o wierzchołkach [i, i+1]. Jeśli jakiegoś nie ma, to dodajemy go do listy
         * brakujących constraintów.
         */
        PointListType const &allConstraints = input.getConstraint (constraintNumber);
        size_t inputSize = allConstraints.size ();
        size_t constraintOffset = input.getConstraintOffset (constraintNumber);

        for (size_t io = 0; io < inputSize - (size_t)(!closed); ++io) {
                size_t i = io + constraintOffset;
                size_t j = ((io + 1) % inputSize) + constraintOffset;

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
void Triangulation <PointArg, TriangleArg, PointList>::addMissingSegments (TriangleEdgeList const &missingSegments)
{
        // 4. Add missing segments.
        for (typename TriangleEdgeList::const_iterator i = missingSegments.begin (); i != missingSegments.end (); ++i) {
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

                                // TODO w paperze napsali, żeby ją zachowac w kolekcji nowych krawędzi - ale po co?
                                newEdge = newDiagonal;
                        }
                }
        }
}

/****************************************************************************/

//template <
//        typename PointArg,
//        typename TriangleArg,
//        typename PointList
//>
//void Triangulation <PointArg, TriangleArg, PointList>::removeSuperfulous (PointListType const &pointList, bool removeOutside)
//{
//        // 5. Remove superfluous triangles.  Input must be in COUNTER CLOCKWISE order.
//        TriangleVector const &triangulation = index.getTriangulation ();
//
//        for (typename TriangleVector::const_iterator i = triangulation.begin (), e = triangulation.end (); i != e; ++i) {
//                TriangleType &triangle = const_cast <TriangleType &> (*i);
//
//                IndexType ta = a (triangle);
//                IndexType tb = b (triangle);
//                IndexType tc = c (triangle);
//
//                // TODO
//                if (ta < 10 || ta > 13 || tb < 10 || tb > 13 || tc < 10 || tc > 13) {
//                        continue;
//                }
//
//                if (a (triangle) == 0 && b (triangle) == 0 && c (triangle) == 0) {
//                        continue;
//                }
//
//                if (removeOutside ^ triangleInside (triangle, pointList)) {
//                        a (triangle, 0);
//                        b (triangle, 0);
//                        c (triangle, 0);
//                }
//        }
//}


template <
        typename PointArg,
        typename TriangleArg,
        typename PointList
>
void Triangulation <PointArg, TriangleArg, PointList>::removeSuperfulous ()
{
        // 5. Remove superfluous triangles.  Input must be in COUNTER CLOCKWISE order.
        TriangleVector const &triangulation = index.getTriangulation ();

        for (typename TriangleVector::const_iterator i = triangulation.begin (), e = triangulation.end (); i != e; ++i) {
                TriangleType &triangle = const_cast <TriangleType &> (*i);

                IndexType ta = a (triangle);
                IndexType tb = b (triangle);
                IndexType tc = c (triangle);

                // Sprawdź czy wierzchołki trójkąta naeżą do jednego konstraintu.
                size_t ca = input.getConstraintForIndex (ta);
                size_t cb = input.getConstraintForIndex (tb);
                size_t cc = input.getConstraintForIndex (tc);

                // Jesli nie, to continue.
                if (ca != cb || cb != cc) {
                        continue;
                }

                if (ta == 0 && tb == 0 && tc == 0) {
                        continue;
                }

                // Jeśli tak, to pobierz referencję do tego konstraintu.
                PointListType const &pointList = input.getConstraint (ca);
                bool removeOutside = (ca == 0);

                if (removeOutside ^ triangleInside (triangle, pointList)) {
                        a (triangle, 0);
                        b (triangle, 0);
                        c (triangle, 0);
                }
        }

        // Clear triangles marked as to-delete.
        index.clean ();
}

/****************************************************************************/

template <
        typename PointArg,
        typename TriangleArg,
        typename PointList
>
bool Triangulation <PointArg, TriangleArg, PointList>::diagonalInside (PointType const &a, PointType const &ap, PointType const &an, PointType const &b) const
{
        int apx = getX (ap) - getX (a);
        int apy = getY (ap) - getY (a);
        int anx = getX (an) - getX (a);
        int any = getY (an) - getY (a);
        int bx = getX (b) - getX (a);
        int by = getY (b) - getY (a);

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
bool Triangulation <PointArg, TriangleArg, PointList>::diagonalInside (TriangleEdgeType const &e, PointListType const &pointList) const
{
        // If e is one of constraints, we don't need to perform furher computations.
        if (std::abs (e.a - e.b) == 1) {
                return true;
        }

        size_t pointsSize = pointList.size ();
        PointType const &a = pointList[e.a];
        PointType const &ap = pointList[(e.a == 0) ? (pointsSize - 1) : (e.a - 1)];
        PointType const &an = pointList[(e.a == pointsSize - 1) ? (0) : (e.a + 1)];
        PointType const &b = pointList[e.b];

        return diagonalInside (a, ap, an, b);
}

/****************************************************************************/

template <
        typename PointArg,
        typename TriangleArg,
        typename PointList
>
bool Triangulation <PointArg, TriangleArg, PointList>::triangleInside (TriangleType const &t, PointListType const &pointList) const
{
        for (int i = 1; i <= 3; ++i) {
                if (!diagonalInside (getEdge (t, static_cast <SideEnum> (i)), pointList)) {
                        return false;
                }
        }

        return true;
}

/****************************************************************************/

template <
        typename PointArg,
        typename TriangleArg,
        typename PointList
>
void Triangulation <PointArg, TriangleArg, PointList>::establishSuperTriangle (PointListType const &p)
{
        std::vector <CoordinateType> coordinate;
        coordinate.reserve (p.size ());

        // TODO nie kopiowac tych współrzędznych, tylko liniowe wyszukiwanie zrobić. Przecież zajmie tyle samo czasu, a nie zajmie w ogóle pamięci.
        for (typename PointListType::const_iterator i = p.begin (), e = p.end (); i != e; ++i) {
                coordinate.push_back (getX (*i));
        }

        std::sort (coordinate.begin (), coordinate.end ());

        CoordinateType xmin = coordinate.front ();
        CoordinateType xmax = coordinate.back ();
        CoordinateType dx = (xmax - xmin) / 4.0;
        xmin -= dx;
        xmax += dx;

        int j = 0;
        for (typename PointListType::const_iterator i = p.begin (), e = p.end (); i != e; ++i, ++j) {
                coordinate[j] = getY (*i);
        }

        std::sort (coordinate.begin (), coordinate.end ());

        CoordinateType ymin = coordinate.front ();
        CoordinateType ymax = coordinate.back ();
        CoordinateType dy = (ymax - ymin) / 4.0;
        ymin -= dy;
        ymax += dy;

#if 0
        std::cerr << xmin << ", " << xmax << ", " << ymin << ", " << ymax << std::endl;
#endif

        PointType pt;
        setX (pt, xmin);
        setY (pt, ymin);
        superTriangle.push_back (pt);

        setX (pt, xmax);
        setY (pt, ymin);
        superTriangle.push_back (pt);

        setX (pt, xmax);
        setY (pt, ymax);
        superTriangle.push_back (pt);

        setX (pt, xmin);
        setY (pt, ymax);
        superTriangle.push_back (pt);

        pointsCollection.push_back (&superTriangle);
        size_ += superTriangle.size ();
}

} // namespace

#endif /* DELAUNAYTRIANGULATION_H_ */
