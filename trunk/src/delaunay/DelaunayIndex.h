/****************************************************************************
 *                                                                          *
 *  Author : lukasz.iwaszkiewicz@gmail.com                                  *
 *  ~~~~~~~~                                                                *
 *  License : see COPYING file for details.                                 *
 *  ~~~~~~~~~                                                               *
 ****************************************************************************/

#ifndef DELAUNAYINDEX_H_
#define DELAUNAYINDEX_H_

#include "DelaunayTriangle.h"
#include "DelaunayPoint.h"
#include "DelaunayEdge.h"
#include "DelaunayTriangleEdge.h"
#include "DelaunayTraits.h"
#include "InputCollection.h"
#include <list>
#include <boost/tuple/tuple.hpp>
#include <map>

namespace Delaunay {

/*
 * This is mutable, R/W index, thus most of ots methods are non const. It is meant to be a private
 * implementation detail not accessible to user.
 */
template <
        typename PointArg = Point,
        typename TriangleArg = Triangle,
        template<typename, typename> class PointList = std::vector,
        template<typename, typename> class ConstraintList = std::vector,
        template<typename> class PointAlloc = std::allocator,
        template<typename> class ConstraintAlloc = std::allocator
>
class DelaunayIndex {
public:

        typedef PointArg PointType;
        typedef TriangleArg TriangleType;
        typedef PointList <PointArg, PointAlloc <PointArg> > PointListType;
        typedef ConstraintList <PointListType, ConstraintAlloc <PointListType> > ConstraintListType;
        typedef InputCollection <PointArg, TriangleArg, PointList, ConstraintList, PointAlloc, ConstraintAlloc> InputCollectionType;
        typedef PointTraits<PointType> PointTraitsType;
        typedef typename PointTraitsType::CoordinateType CoordinateType;
        typedef Edge <PointType> EdgeType;
        typedef TriangleTraits <TriangleType> TriangleTraitsType;
        typedef typename TriangleTraitsType::IndexType IndexType;
        typedef TriangleEdge<TriangleType> TriangleEdgeType;
        typedef std::list <TriangleEdgeType> TriangleEdgeList;
        typedef std::vector <TriangleEdgeType> TriangleEdgeVector;
        typedef std::vector <TriangleType> TriangleVector;
        typedef std::vector <TriangleType *> TrianglePtrVector;
        typedef std::vector <TrianglePtrVector> TriangleIndex;
        typedef boost::tuple <int, SideEnum, SideEnum> IntersectionInfo;

        /**
         *
         */
//        class HalfEdge {
//        public:
//
//                HalfEdge (IndexType b_ = 0, TriangleType *t_ = 0) : vertexB (b_), triangle (t_), twin (0), next (0) {}
//
//                IndexType getVertexB () const { return vertexB; }
//                IndexType getVertexC (IndexType a) const
//                {
//                        SideEnum cSide = otherThan (getVertexSide (*triangle, a), getVertexSide (*triangle, vertexB));
//                        return getVertex (*triangle, cSide);
//                }
//
//                TriangleType const *getTriangle () const { return triangle; }
//                TriangleType *getTriangle () { return triangle; }
//
//                HalfEdge const *getTwin () const { return twin; }
//                HalfEdge *getTwin () { return twin; }
//
//                HalfEdge const *getNext () const { return next; }
//                HalfEdge *getNext () { return next; }
//
//        private:
//
//                friend class DelaunayIndex;
//                IndexType vertexB;
//                TriangleType *triangle;
//                HalfEdge *twin;
//                HalfEdge *next;
//        };

//        /*
//         * http://www.lafstern.org/matt/col1.pdf : Why you shouldn't use set (and what you should use instead) by Matt Austern
//         */
//        typedef std::vector <HalfEdge> HalfEdgeVector;
//        typedef std::list <HalfEdge> HalfEdgeList;
//
//        struct HalfEdgeNode {
//                HalfEdgeNode () : first (0), last (0) {}
//
//                HalfEdge *first;
//                HalfEdge *last;
//                HalfEdgeList all;
//        };

//        typedef std::vector <HalfEdgeNode> HalfEdgeIndex;

        DelaunayIndex (InputCollectionType const &i) : input (i)
        {
        }

        void reserve ()
        {
                triangleIndex.resize (input.size ());
                // TODO ***KRYTYCZNE*** ustawić tu tyle ile ma być. Da się to wyliczyć na początku!?
                triangulation.reserve (input.size () * 10);
//                adjacentTrianglesIndex.reserve (triangulation.size ());
        }

        TriangleVector const &getTriangulation () const { return triangulation; }
        size_t getNumTriangles () const { return triangulation.size (); }
        void addTriangle (IndexType index, TriangleType *triangle);
        void addTriangle (TriangleType &triangle);
        TriangleType const &getTriangle (size_t i) const { return triangulation.at (i); }
        TriangleType &getTriangle (size_t i) { return triangulation.at (i); }

        size_t getTriangleIndexSize () const { return triangleIndex.size (); }
        void clean ();

        void setVertex (TriangleType &t, SideEnum s, IndexType v);

        TrianglePtrVector &getTrianglesForPoint (IndexType i) { return triangleIndex[i]; }
//        std::pair <TriangleType const *, TriangleType const *> getTrianglesForEdge (TriangleEdge const &e) const;
//        std::pair <TriangleType *, TriangleType *> getTrianglesForEdge (TriangleEdge const &e);
//
//        HalfEdge *findEdge (IndexType a, IndexType b);
//        HalfEdge const *findEdge (IndexType a, IndexType b) const { return const_cast <DelaunayIndex *> (this)->findEdge (a, b); }
//
//        HalfEdge *findEdge (IndexType a, IndexType b, IndexType c);
//        HalfEdge const *findEdge (IndexType a, IndexType b, IndexType c) const { return const_cast <DelaunayIndex *> (this)->findEdge (a, b, c); }

        typedef std::pair <TriangleType const *, TriangleType const *> ConstTrianglePair;
        typedef std::pair <TriangleType *, TriangleType *> TrianglePair;

        ConstTrianglePair getTrianglesForEdge (TriangleEdgeType const &e) const { return const_cast <DelaunayIndex *> (this)->getTrianglesForEdge (e); }
        TrianglePair getTrianglesForEdge (TriangleEdgeType const &e);

/****************************************************************************/

        TriangleType *getAdjacentTriangle (TriangleType const &triangle, SideEnum side);
        void setAdjacentTriangle (TriangleType &t, SideEnum s, TriangleType *a);

        /**
         * Input : this->triangulation (performed earlier), edge to check for. Edge must have endpoints
         * in this->input.
         * Output : list of edges which intersects with edge.
         * Non const because it returns R-W data by output parameter 'crossingTriangles'.
         */
        void findCrossingEdges (TriangleEdgeType const &edge, TriangleEdgeList *crossingEdges, TrianglePtrVector *crossingTriangles);

        /**
         * First tuple element : number of intersections (0, 1 or 2).
         * Second tuple element : number of first edge (if any) which intersects with edge e.
         * Third tuple element : number of second edge (if any) which intersects with edge e.
         * - 1 : c-b
         * - 2 : c-a
         * - 3 : b-a
         */
        IntersectionInfo intersects (TriangleType const &t, EdgeType const &e) const;

        /**
         * Are two adjacent triangles form quadrilateral which is convex?
         */
        bool twoTrianglesConvex (TriangleEdgeType const &e/*, TriangleType const &a, TriangleType const &b*/) const;
        bool twoTrianglesNotDelaunay (TriangleEdgeType const &e) const;

        /**
         * Based on Cline and Renka
         */
        bool pointInCircumcircle (TriangleType const &triangle, IndexType point) const;

        /**
         * Perform a flip, and return new diagonal. Triangle index.
         */
        void flip (TriangleEdgeType const &oldDiagonal, TriangleEdgeType *newDiagonal);

        /**
         * Index based edge to coordinate based edge.
         */
        EdgeType triangleEdgeToEdge (TriangleEdgeType const &e) const
        {
                return EdgeType (input[e.a], input[e.b]);
        }

        /*
         *
         */
        void sortEdgeIndex ();

private:

        /*
         * Sort vertices of triangle in CCW order.
         */
        void sortTriangle (TriangleType &triangle);

        /*
         *
         */
//        class HalfEdgeCompare {
//        public:
//                HalfEdgeCompare (IndexType a_) : a (a_), b (0), c (0), mode (REGULAR) {}
//                HalfEdgeCompare (IndexType a_, IndexType b_) : a (a_), b (b_), c (0), mode (CUSTOM_AB) {}
//                HalfEdgeCompare (IndexType a_, IndexType b_, IndexType c_) : a (a_), b (b_), c (c_), mode (CUSTOM_ABC) {}
//
//                bool operator () (HalfEdge const &e1, HalfEdge const &e2) const
//                {
//                        switch (mode) {
//                        default:
//                        case REGULAR:
//                                return opRegular (e1, e2);
//
//                        case CUSTOM_AB:
//                                return opCustomAB (e1, e2);
//
//                        case CUSTOM_ABC:
//                                return opCustomABC (e1, e2);
//                        }
//                }
//
//                bool opRegular (HalfEdge const &e1, HalfEdge const &e2) const
//                {
//                        if (e1.getVertexB () == e2.getVertexB ()) {
//                                IndexType c1 = e1.getVertexC (a);
//                                IndexType c2 = e2.getVertexC (a);
//                                return c1 < c2;
//                        }
//
//                        return e1.getVertexB () < e2.getVertexB ();
//                }
//
//                bool opCustomAB (HalfEdge const &e1, HalfEdge const &) const
//                {
//                        return e1.getVertexB () < b;
//                }
//
//                bool opCustomABC (HalfEdge const &e1, HalfEdge const &) const
//                {
//                        if (e1.getVertexB () == b) {
//                                IndexType c1 = e1.getVertexC (a);
//                                return c1 < c;
//                        }
//
//                        return e1.getVertexB () < b;
//                }
//
//        private:
//                IndexType a;
//                IndexType b;
//                IndexType c;
//                enum Mode { REGULAR, CUSTOM_AB, CUSTOM_ABC } mode;
//        };

        /*
         *
         */
//        class PHalfEdgeCompare {
//        public:
//                PHalfEdgeCompare (IndexType a_, IndexType b_, IndexType c_) : impl (a_, b_, c_) {}
//                bool operator () (HalfEdge const *e1, HalfEdge const *e2) const { return impl.operator () (*e1, *e2); }
//        private:
//                HalfEdgeCompare impl;
//        };

        /*
         * Works only if triangles are CCW ordered.
         */
        struct TriangleCompare {
                TriangleCompare (IndexType a_) : a (a_), mode (REGULAR) {}
                TriangleCompare (IndexType a_, IndexType b_) : a (a_), b (b_), mode (CUSTOM_AB) {}

                bool operator () (TriangleType const *t1, TriangleType const *t2) const
                {
                        if (mode == REGULAR) {
                                return opRegular (t1, t2);
                        }
                        else {
                                return opCustomAB (t1, t2);
                        }
                }

                bool opRegular (TriangleType const *t1, TriangleType const *t2) const
                {
                        SideEnum a1s = getVertexSide (*t1, a);
                        IndexType b1 = getVertex (*t1, static_cast <SideEnum> ((a1s + 1) % 3));

                        SideEnum a2s = getVertexSide (*t2, a);
                        IndexType b2 = getVertex (*t2, static_cast <SideEnum> ((a2s + 1) % 3));

                        if (b1 == b2) {
                                IndexType c1 = getVertex (*t1, static_cast <SideEnum> ((a1s + 2) % 3));
                                IndexType c2 = getVertex (*t2, static_cast <SideEnum> ((a2s + 2) % 3));
                                return c1 < c2;
                        }

                        return b1 < b2;
                }

                bool opCustomAB (TriangleType const *t1, TriangleType const *) const
                {
                        SideEnum t1s = getVertexSide (*t1, a);
                        IndexType n1 = getVertex (*t1, static_cast <SideEnum> ((t1s + 1) % 3));
                        return n1 < b;
                }

                IndexType a;
                IndexType b;
                enum Mode { REGULAR, CUSTOM_AB } mode;
        };

        struct TraingleRemovePredicate {
                bool operator () (TriangleType const &t) { return !a (t) && !b (t) && !c (t); }
        };

//        void addHalfEdge (IndexType a, IndexType b, TriangleType *t);

//        /**
//         * Sort edges, so thier triangles will be in order. Edges will be stroed in clockwise order.
//         */
//        void topologicalSort (HalfEdgeNode &sortedBC, IndexType a);

private:

        // Input points.
        InputCollectionType const &input;
        // Adjacency list - points -> triangles.
        TriangleIndex triangleIndex;
        // Triangles in triangulation in random order.
        TriangleVector triangulation;

//        typedef TriangleType *AdjacentTriangles[3];

        struct AdjacentTriangles {
                AdjacentTriangles () : tA (0), tB (0), tC (0) {}
                TriangleType *tA, *tB, *tC;
        };
        typedef std::map <TriangleType const *, AdjacentTriangles> AdjacentTrianglesIndex;

        // 3 (or less) adjacent triangles for every single triangle.
        AdjacentTrianglesIndex adjacentTrianglesIndex;

};

/****************************************************************************/

template <
        typename PointArg,
        typename TriangleArg,
        template<typename, typename> class PointList,
        template<typename, typename> class ConstraintList,
        template<typename> class PointAlloc,
        template<typename> class ConstraintAlloc
>
typename DelaunayIndex <PointArg, TriangleArg, PointList,ConstraintList, PointAlloc, ConstraintAlloc>::TriangleType *
DelaunayIndex <PointArg, TriangleArg, PointList,ConstraintList, PointAlloc, ConstraintAlloc>::getAdjacentTriangle (TriangleType const &t, SideEnum side)
{
        AdjacentTriangles &at = adjacentTrianglesIndex[&t];

        switch (side) {
        case A:
                return at.tA;
        case B:
                return at.tB;
        case C:
                return at.tC;
        default:
                return 0;
        }

}

/****************************************************************************/

template <
        typename PointArg,
        typename TriangleArg,
        template<typename, typename> class PointList,
        template<typename, typename> class ConstraintList,
        template<typename> class PointAlloc,
        template<typename> class ConstraintAlloc
>
void DelaunayIndex <PointArg, TriangleArg, PointList,ConstraintList, PointAlloc, ConstraintAlloc>::setAdjacentTriangle (TriangleType &t, SideEnum s, TriangleType *a)
{
        AdjacentTriangles &at = adjacentTrianglesIndex[&t];

        switch (s) {
        case A:
                at.tA = a;
                break;
        case B:
                at.tB = a;
                break;
        case C:
                at.tC = a;
                break;
        default:
                break;
        }
}

/****************************************************************************/

template <
        typename PointArg,
        typename TriangleArg,
        template<typename, typename> class PointList,
        template<typename, typename> class ConstraintList,
        template<typename> class PointAlloc,
        template<typename> class ConstraintAlloc
>
typename DelaunayIndex <PointArg, TriangleArg, PointList,ConstraintList, PointAlloc, ConstraintAlloc>::IntersectionInfo
DelaunayIndex <PointArg, TriangleArg, PointList,ConstraintList, PointAlloc, ConstraintAlloc>::intersects (TriangleType const &t, EdgeType const &e) const
{
        IntersectionInfo ret;

        int cnt = 0;
        EdgeType edge = triangleEdgeToEdge (TriangleEdgeType (c (t), b (t)));

        if (Delaunay::intersects (e, edge)) {
                ret.get<1> () = A;
                ++cnt;
        }

        edge.a = input[t.c];
        edge.b = input[t.a];
        if (Delaunay::intersects (e, edge)) {
                if (cnt) {
                        ret.get<2> () = B;
                }
                else {
                        ret.get<1> () = B;
                }
                ++cnt;
        }

        edge.a = input[t.b];
        edge.b = input[t.a];
        if (Delaunay::intersects (e, edge)) {
                if (cnt) {
                        ret.get<2> () = C;
                }
                else {
                        ret.get<1> () = C;
                }
                ++cnt;
        }

        ret.get<0> () = cnt;
        return ret;
}

/****************************************************************************/

template <
        typename PointArg,
        typename TriangleArg,
        template<typename, typename> class PointList,
        template<typename, typename> class ConstraintList,
        template<typename> class PointAlloc,
        template<typename> class ConstraintAlloc
>
void DelaunayIndex <PointArg, TriangleArg, PointList,ConstraintList, PointAlloc, ConstraintAlloc>::findCrossingEdges (TriangleEdgeType const &edge, TriangleEdgeList *crossingEdges, TrianglePtrVector *crossingTriangles)
{
        TrianglePtrVector const &incidentTriangles = triangleIndex[edge.a];
        EdgeType e = triangleEdgeToEdge (edge);
        TriangleType *start = NULL;
        IntersectionInfo intersections;

        for (typename TrianglePtrVector::const_iterator k = incidentTriangles.begin (); k != incidentTriangles.end (); ++k) {
                intersections = intersects (**k, e);
                if (intersections.get<0> ()) {
                        start = *k;
                        break;
                }
        }

#if 0
        if (start) {
                std::cerr << "Missing constraint : (" << edge.first << "," << edge.second << "), first triangle intersecinting this constraint : " << *start << std::endl;
        }
#endif

        assert (start);

        if (crossingTriangles) {
                crossingTriangles->push_back (start);
        }

        SideEnum commonEdgeNumber = intersections.get<1> ();
        TriangleType *next = start;

        while (true) {
                TriangleEdgeType commonEdge = getEdge (*next, commonEdgeNumber);
                next = getAdjacentTriangle (*next, commonEdgeNumber);
                commonEdgeNumber = getEdgeSide (*next, commonEdge);

                // Musi być, bo gdyby nie było, to by znaczyło, że constraint wychodzi poza zbiór punktów (ma jeden koniec gdzieś w powietrzu).
                assert (next);

                if (crossingTriangles) {
//                        crossingEdges->push_back (crossingEdge);
                        crossingEdges->push_back (commonEdge);
                        crossingTriangles->push_back (next);
                }

                if (hasVertex (*next, edge.b)) {
                        break;
                }

                intersections = intersects (*next, e);

                // Musi się przecinać w 2 punktach, bo inaczej by wyszło z funkcji.
                assert (intersections.get<0> () == 2);

                // Eliminate commonEdge from equation - we know about it already. Find new commonEdge.
                if (intersections.get<1> () == commonEdgeNumber) {
                        commonEdgeNumber = intersections.get<2> ();
                }
                else if (intersections.get<2> () == commonEdgeNumber) {
                        commonEdgeNumber = intersections.get<1> ();
                }
        }

}

/****************************************************************************/

template <
        typename PointArg,
        typename TriangleArg,
        template<typename, typename> class PointList,
        template<typename, typename> class ConstraintList,
        template<typename> class PointAlloc,
        template<typename> class ConstraintAlloc
>
bool DelaunayIndex <PointArg, TriangleArg, PointList,ConstraintList, PointAlloc, ConstraintAlloc>::twoTrianglesConvex (TriangleEdgeType const &firstDiagonal/*, TriangleType const &a, TriangleType const &b*/) const
{
        ConstTrianglePair pair = getTrianglesForEdge (firstDiagonal);
        TriangleType const *a = pair.first;
        TriangleType const *b = pair.second;
        if (!a || !b) {
                std::cerr << a << ", " << b << ", " << firstDiagonal << std::endl;
        }
        assert (a && b);

        SideEnum aSide = getEdgeSide (*a, firstDiagonal);
        SideEnum bSide = getEdgeSide (*b, firstDiagonal);

        TriangleEdgeType secondDiagonal = TriangleEdgeType (getVertex (*a, aSide), getVertex (*b, bSide));
        EdgeType e1 = triangleEdgeToEdge (firstDiagonal);
        EdgeType e2 = triangleEdgeToEdge (secondDiagonal);
        return Delaunay::intersects (e1, e2);
}

/****************************************************************************/

template <
        typename PointArg,
        typename TriangleArg,
        template<typename, typename> class PointList,
        template<typename, typename> class ConstraintList,
        template<typename> class PointAlloc,
        template<typename> class ConstraintAlloc
>
bool DelaunayIndex <PointArg, TriangleArg, PointList,ConstraintList, PointAlloc, ConstraintAlloc>::pointInCircumcircle (TriangleType const &triangle, IndexType point) const
{
        PointType const &ta = input[a (triangle)];
        PointType const &tb = input[b (triangle)];
        PointType const &tc = input[c (triangle)];
        PointType const &tp = input[point];

        double cosa = (ta.x - tc.x) * (tb.x - tc.x) + (ta.y - tc.y) * (tb.y - tc.y);
        double cosb = (tb.x - tp.x) * (ta.x - tp.x) + (tb.y - tp.y) * (ta.y - tp.y);

        if (cosa >= 0 && cosb >= 0) {
                return false;
        }

        if (cosa < 0 && cosb < 0) {
                return true;
        }

        double sinab = ((ta.x - tc.x) * (tb.y - tc.y) - (tb.x - tc.x) * (ta.y - tc.y)) * cosb + ((tb.x - tp.x) * (ta.y - tp.y) - (ta.x - tp.x) * (tb.y - tp.y)) * cosa;

        if (sinab < 0) {
                return true;
        }

        return false;
}

/****************************************************************************/

template <
        typename PointArg,
        typename TriangleArg,
        template<typename, typename> class PointList,
        template<typename, typename> class ConstraintList,
        template<typename> class PointAlloc,
        template<typename> class ConstraintAlloc
>
bool DelaunayIndex <PointArg, TriangleArg, PointList,ConstraintList, PointAlloc, ConstraintAlloc>::twoTrianglesNotDelaunay (TriangleEdgeType const &firstDiagonal) const
{
        ConstTrianglePair pair = getTrianglesForEdge (firstDiagonal);
        TriangleType const *a = pair.first;
        TriangleType const *b = pair.second;

        SideEnum aSide = getEdgeSide (*a, firstDiagonal);
        SideEnum bSide = getEdgeSide (*b, firstDiagonal);

        return pointInCircumcircle (*a, getVertex (*b, bSide)) || pointInCircumcircle (*b, getVertex (*a, aSide));
}

/****************************************************************************/

template <
        typename PointArg,
        typename TriangleArg,
        template<typename, typename> class PointList,
        template<typename, typename> class ConstraintList,
        template<typename> class PointAlloc,
        template<typename> class ConstraintAlloc
>
void DelaunayIndex <PointArg, TriangleArg, PointList,ConstraintList, PointAlloc, ConstraintAlloc>::addTriangle (IndexType index, TriangleType *triangle)
{
        assert (getTriangleIndexSize () > index);
        triangleIndex[index].push_back (triangle);
}

/****************************************************************************/

template <
        typename PointArg,
        typename TriangleArg,
        template<typename, typename> class PointList,
        template<typename, typename> class ConstraintList,
        template<typename> class PointAlloc,
        template<typename> class ConstraintAlloc
>
void DelaunayIndex <PointArg, TriangleArg, PointList,ConstraintList, PointAlloc, ConstraintAlloc>::addTriangle (TriangleType &triangle)
{
        sortTriangle (triangle);
        triangulation.push_back (triangle);

        // Update triangle index.
        TriangleType *t = &triangulation.back ();

        addTriangle (a (*t), t);
        addTriangle (b (*t), t);
        addTriangle (c (*t), t);
}

/****************************************************************************/

template <
        typename PointArg,
        typename TriangleArg,
        template<typename, typename> class PointList,
        template<typename, typename> class ConstraintList,
        template<typename> class PointAlloc,
        template<typename> class ConstraintAlloc
>
void DelaunayIndex <PointArg, TriangleArg, PointList,ConstraintList, PointAlloc, ConstraintAlloc>::sortTriangle (TriangleType &triangle)
{
        IndexType t1 = a (triangle);
        IndexType t2 = b (triangle);
        IndexType t3 = c (triangle);

        PointType const &ta = input[t1];
        PointType const &tb = input[t2];
        PointType const &tc = input[t3];

        double det = getX (ta) * getY (tb) + getX (tb) * getY (tc) + getX (tc) * getY (ta) - getY (tb) * getX (tc) - getY (tc) * getX (ta) - getY (ta) * getX (tb);

        if (det > 0) { // already CCW
                return;
        }

        // Not CCW - swap
        b (triangle, t3);
        c (triangle, t2);
}

/****************************************************************************/

template <
        typename PointArg,
        typename TriangleArg,
        template<typename, typename> class PointList,
        template<typename, typename> class ConstraintList,
        template<typename> class PointAlloc,
        template<typename> class ConstraintAlloc
>
void DelaunayIndex <PointArg, TriangleArg, PointList,ConstraintList, PointAlloc, ConstraintAlloc>::sortEdgeIndex ()
{
        IndexType a = 0;
        for (typename TriangleIndex::iterator i = triangleIndex.begin (), e = triangleIndex.end (); i != e; ++i, ++a) {
                TrianglePtrVector &trianglesForIndex = *i;
                std::sort (trianglesForIndex.begin (), trianglesForIndex.end (), TriangleCompare (a));
        }

//        std::cerr << "--------TRIANGLE INDEX--------" << std::endl;
//        std::cerr << triangleIndex << std::endl;
//        std::cerr << "--------EDGE INDEX--------" << std::endl;
//
//        int cnt = 0;
//        for (typename HalfEdgeIndex::const_iterator i = edgeIndex.begin (); i != edgeIndex.end (); ++i, ++cnt) {
//                HalfEdgeVector const &edgesForIndex = *i;
//
//                std::cerr << cnt << ". ";
//
//                for (typename HalfEdgeVector::const_iterator j = edgesForIndex.begin (); j != edgesForIndex.end (); ++j) {
//                        HalfEdge const &edge = *j;
//                        std::cerr << edge.getVertexB () << ":" << *edge.getTriangle () << " | ";
//                }
//
//                std::cerr << std::endl;
//        }
}

/****************************************************************************/

#if 0
template <typename Input, typename Traits>
void DelaunayIndex<Input, Traits>::topologicalSort (HalfEdgeVector &input, IndexType a)
{
#if 0
        std::cerr << "#### " << a << std::endl;
#endif
        size_t initialSize = input.size ();
        HalfEdgeList lex (input.begin (), input.end ());
        lex.sort (HalfEdgeCompare (a));

        // Find first (and the only, if any) pair of consecutive HalfEdges whose b-vertices don't match.
        // Size is always even.
        typename HalfEdgeList::iterator i = lex.begin ();
        for (; i != lex.end (); ++i, ++i) {
                typename HalfEdgeList::iterator j = i;
                ++j;
                if (i->getVertexB () != j->getVertexB ()) {
                        break;
                }
        }

        /*
         *  If found - that means that this "fan" of triangles is not closed, and we must start
         *  sorting from one of its ends. If fan is closed (i.e. there is no gap between triangles),
         *  we can start from any point.
         */
        if (i == lex.end ()) {
                i = lex.begin ();
        }
#if 0
        for (typename HalfEdgeList::const_iterator j = lex.begin (); j != lex.end (); ++j) {
                HalfEdge const &edge = *j;
                std::cerr << edge.getVertexB () << ":" << edge.getVertexC (a) << " | ";
        }
        std::cerr << std::endl;
#endif
        input.clear ();

        HalfEdge edge = *i;
        input.push_back (edge);
#if 0
        std::cerr << edge.getVertexB() << "," << edge.getVertexC(a) << std::endl;
#endif
        lex.erase (i);

        bool directionDown = true;
        while (true) {
                typename HalfEdgeList::iterator j = std::lower_bound (lex.begin (), lex.end (), HalfEdge (), HalfEdgeCompare (a, edge.getVertexC (a), edge.getVertexB ()));
                edge = *j;
                input.push_back (edge);
#if 0
                std::cerr << edge.getVertexB() << "," << edge.getVertexC(a) << std::endl;
#endif
                if (input.size () >= initialSize) {
                        break;
                }

                typename HalfEdgeList::iterator down = j;
                ++down;
                typename HalfEdgeList::iterator up = (j != lex.begin ()) ? (j) : (lex.end ());
                --up;
                lex.erase (j);

#if 0
                std::cerr << "down : " << down->getVertexB () << ", up : " << up->getVertexB () << std::endl;
#endif

                retry:
                // Znajdz następny (up, lub down)
                if (directionDown) {
                        if (down == lex.end () || down->getVertexB () != edge.getVertexB ()) {
                                directionDown = false;
                                goto retry;
                        }

                        edge = *down;
                        lex.erase (down);
                }
                else {
                        if (up == lex.end () || up->getVertexB () != edge.getVertexB ()) {
                                directionDown = true;
                                goto retry;
                        }

                        edge = *up;
                        lex.erase (up);
                }

                input.push_back (edge);
#if 0
                std::cerr << edge.getVertexB() << "," << edge.getVertexC(a) << std::endl;
#endif
        }
#if 0
        for (typename HalfEdgeVector::const_iterator j = input.begin (); j != input.end (); ++j) {
                HalfEdge const &edge = *j;
                std::cerr << edge.getVertexB () << ":" << edge.getVertexC (a) << " | ";
        }
        std::cerr << std::endl;
#endif
}
#endif

#if 0
/**
 * TODO pozbyć się kolekcji sorted i sortowac inplace w lex. Da się chyba.
 * TODO Powiązać HalfEdge w druga stronę (nie wiem po co, ale mają wskaxniki, to powiązać).
 * TODO Zamknąć kółko. Next osrtatniego powinien wskazywac na następny za pierwszym jesli pierwszy i ostatni mają to samo B.
 */
template <typename Input, typename Traits>
void DelaunayIndex<Input, Traits>::topologicalSort (HalfEdgeNode &node, IndexType a)
{
        HalfEdgeList &all = node.all;
        size_t initialSize = all.size ();
        all.sort (HalfEdgeCompare (a));

        typedef std::list <HalfEdge *> PHalfEdgeList;
        typedef std::vector <HalfEdge *> PHalfEdgeVector;
        PHalfEdgeList lex;

        for (typename HalfEdgeList::iterator i = all.begin (), e = all.end (); i != e; ++i) {
                lex.push_back (&*i);
        }

        PHalfEdgeVector sorted;
        sorted.reserve (initialSize);

        // Find first (and the only, if any) pair of consecutive HalfEdges whose b-vertices don't match.
        // Size is always even.
        typename PHalfEdgeList::iterator i = lex.begin ();
        for (; i != lex.end (); ++i, ++i) {
                typename PHalfEdgeList::iterator j = i;
                ++j;
                if ((*i)->getVertexB () != (*j)->getVertexB ()) {
                        break;
                }
        }

        /*
         *  If found - that means that this "fan" of triangles is not closed, and we must start
         *  sorting from one of its ends. If fan is closed (i.e. there is no gap between triangles),
         *  we can start from any point.
         */
        if (i == lex.end ()) {
                i = lex.begin ();
        }

        HalfEdge *edge = *i;
        sorted.push_back (edge);
        lex.erase (i);

        bool directionDown = true;
        while (true) {
                typename PHalfEdgeList::iterator j = std::lower_bound (lex.begin (), lex.end (), (HalfEdge *)0, PHalfEdgeCompare (a, edge->getVertexC (a), edge->getVertexB ()));
                edge = *j;
                sorted.push_back (edge);

                if (sorted.size () >= initialSize) {
                        break;
                }

                typename PHalfEdgeList::iterator down = j;
                ++down;
                typename PHalfEdgeList::iterator up = (j != lex.begin ()) ? (j) : (lex.end ());
                --up;
                lex.erase (j);

                retry:
                // Znajdz następny (up, lub down)
                if (directionDown) {
                        if (down == lex.end () || (*down)->getVertexB () != edge->getVertexB ()) {
                                directionDown = false;
                                goto retry;
                        }

                        edge = *down;
                        lex.erase (down);
                }
                else {
                        if (up == lex.end () || (*up)->getVertexB () != edge->getVertexB ()) {
                                directionDown = true;
                                goto retry;
                        }

                        edge = *up;
                        lex.erase (up);
                }

                sorted.push_back (edge);
        }

#if 1
        std::ostringstream o1, o2;
        std::cerr << "SORTED for " << a << " : ";
        for (typename PHalfEdgeVector::const_iterator j = sorted.begin (); j != sorted.end (); ++j) {
                HalfEdge const *edge = *j;
                std::cerr << edge->getVertexB () << ":" << edge->getVertexC (a) << " | ";
                o1 << edge->getVertexB () << ":" << edge->getVertexC (a) << " | ";
        }
        std::cerr << std::endl;
#endif

        // 2. Poustawiaj wskaźniki.
        for (typename PHalfEdgeVector::iterator i = sorted.begin (), e = sorted.end (); i != e;) {
                HalfEdge *current = *i;
                HalfEdge *next = (i + 1 != e) ? (*(i + 1)) : (0);

                if (!next) {
                        break;
                }

                // Zwiększanie iteratora
                bool match = current->getVertexB () == next->getVertexB ();
                i += (match) ? (2) : (1);

                // Ustawienie nextów:
                if (i != e) {
                        current->next = *i;
                }

                // Ustawienie twinów:
                if (match) {
                        current->twin = next;
                        next->twin = current;
                }
        }

        node.first = sorted.front ();
        node.last = sorted.back ();

        if (node.first->getVertexB () == node.last->getVertexB ()) {
                node.first->twin = node.last;
                node.last->twin = node.first;
        }

#if 0
        {
                std::cerr << "DATA   : ";
                // first
                HalfEdge *edge = &sorted.front ();

                do {
                        std::cerr << edge->getVertexB () << ":" << edge->getVertexC (a) << " | ";
                        o2 << edge->getVertexB () << ":" << edge->getVertexC (a) << " | ";

                        if (edge->twin) {
                                std::cerr << edge->twin->getVertexB () << ":" << edge->twin->getVertexC (a) << " | ";
                                o2 << edge->twin->getVertexB () << ":" << edge->twin->getVertexC (a) << " | ";
                        }
                }
                while ((edge = edge->next));
                std::cerr << std::endl;

                if (std::string (o1.str ()) !=  std::string (o2.str ())) {
                        std::cerr << "AAAAAAAAAAAA" << std::endl;
                        exit (0);
                }
        }
#endif
}
#endif

/*
 *      top
 *       /\
 *      /  \                /|\
 *      ----     ->        / | \ right
 *      \  /          left \ | /
 *       \/                 \|/
 *     bottom
 *
 */
template <
        typename PointArg,
        typename TriangleArg,
        template<typename, typename> class PointList,
        template<typename, typename> class ConstraintList,
        template<typename> class PointAlloc,
        template<typename> class ConstraintAlloc
>
void DelaunayIndex <PointArg, TriangleArg, PointList,ConstraintList, PointAlloc, ConstraintAlloc>::flip (TriangleEdgeType const &oldDiagonal, TriangleEdgeType *newDiagonal)
{
        TrianglePair pair = getTrianglesForEdge (oldDiagonal);
        TriangleType *top = pair.first;
        TriangleType *bottom = pair.second;
#if 0
        std::cerr << "oldDiagonal : " << oldDiagonal << ", top : ";
        if (top) { std::cerr << *top; } else { std::cerr << "NULL"; }
        std::cerr << ", bottom : ";
        if (bottom) { std::cerr << *bottom; } else { std::cerr << "NULL"; }
        std::cerr << ", ptr : " << top << ", " << bottom << std::endl;
#endif
        if (!bottom) {
                *newDiagonal = oldDiagonal;
                std::cerr << "!s - return" << std::endl;
                return;
        }

        SideEnum topLeftVertexSide = getVertexSide (*top, oldDiagonal.a);
        SideEnum topRightVertexSide = getVertexSide (*top, oldDiagonal.b);
        SideEnum topTopVertexSide = otherThan (topLeftVertexSide, topRightVertexSide);
        IndexType topTopVertex = getVertex (*top, topTopVertexSide);
#if 0
        std::cerr << "foa : " << (int)topLeftVertexSide << ", fob : " << (int)topRightVertexSide << ", fc : " << (int)topTopVertexSide << ", fcIndex : " << topTopVertex << std::endl;
#endif
        SideEnum bottomLeftVertexSide = getVertexSide (*bottom, oldDiagonal.a);
        SideEnum bottomRightVertexSide = getVertexSide (*bottom, oldDiagonal.b);
        SideEnum bottomBottomVertexSide = otherThan (bottomRightVertexSide, bottomLeftVertexSide);
        IndexType bottomBottomVertex = getVertex (*bottom, bottomBottomVertexSide);
#if 0
        std::cerr << "soa : " << (int)bottomRightVertexSide << ", sob : " << (int)bottomLeftVertexSide << ", sc : " << (int)bottomBottomVertexSide << ", scIndex : " << bottomBottomVertex << std::endl;
#endif

        // Store before modification.
        TriangleType *tmp = getAdjacentTriangle (*top, topLeftVertexSide);

        // top becomes left.
        this->setVertex (*top, topRightVertexSide, bottomBottomVertex);
        setAdjacentTriangle (*top, topTopVertexSide, getAdjacentTriangle (*bottom, bottomRightVertexSide));
        setAdjacentTriangle (*top, topLeftVertexSide, bottom);

        // bottom becomes right.
        this->setVertex (*bottom, bottomLeftVertexSide, topTopVertex);
        setAdjacentTriangle (*bottom, bottomBottomVertexSide, tmp);
        setAdjacentTriangle (*bottom, bottomRightVertexSide, top);

        newDiagonal->a = topTopVertex;
        newDiagonal->b = bottomBottomVertex;
#if 0
        std::cerr << "newDiagonal : " << *newDiagonal << ", top : " << *top << ", bottom : " << *bottom << std::endl;
#endif
}

/****************************************************************************/

template <
        typename PointArg,
        typename TriangleArg,
        template<typename, typename> class PointList,
        template<typename, typename> class ConstraintList,
        template<typename> class PointAlloc,
        template<typename> class ConstraintAlloc
>
typename DelaunayIndex <PointArg, TriangleArg, PointList,ConstraintList, PointAlloc, ConstraintAlloc>::TrianglePair
DelaunayIndex <PointArg, TriangleArg, PointList,ConstraintList, PointAlloc, ConstraintAlloc>::getTrianglesForEdge (TriangleEdgeType const &e)
{
#if 0
        TrianglePtrVector &trianglesForIndex = triangleIndex[e.a];
        typename TrianglePtrVector::iterator i = std::lower_bound (trianglesForIndex.begin (), trianglesForIndex.end (), static_cast <TriangleType*> (0), TriangleCompare (e.a, e.b));

        if (i == trianglesForIndex.end ()) {
                return TrianglePair ();
        }

        TrianglePair foundTriangles;
        foundTriangles.first = *i;
        SideEnum aSide = getVertexSide (*foundTriangles.first, e.a);
        SideEnum adj = static_cast <SideEnum> ((aSide - 1) % 3);
        foundTriangles.second = getAdjacentTriangle (*foundTriangles.first, adj);
        return foundTriangles;
#else
        TrianglePtrVector const &triaglesA = triangleIndex[e.a];
        TrianglePair foundTriangles;
#if 0
        std::cerr << "Edge : " << e << std::endl;
#endif
        for (typename TrianglePtrVector::const_iterator i = triaglesA.begin (); i != triaglesA.end (); ++i) {
                TriangleType *t = *i;
#if 0
                std::cerr << "getTrForE : " << *t;
#endif
                SideEnum s = getVertexSide (*t, e.a);
                TriangleEdgeType me = getEdge (*t, s);

                if (me.a == e.b || me.b == e.b) {
                        if (!foundTriangles.first) {
                                foundTriangles.first = t;
#if 0
                                std::cerr << " +++a";
#endif
                        }
                        else if (!foundTriangles.second) {
                                foundTriangles.second = t;
#if 0
                                std::cerr << " +++b"  << std::endl;
#endif
                                return foundTriangles;
                        }
                }
#if 0
                std::cerr  << std::endl;
#endif
        }

        return foundTriangles;
#endif
}

/****************************************************************************/

template <
        typename PointArg,
        typename TriangleArg,
        template<typename, typename> class PointList,
        template<typename, typename> class ConstraintList,
        template<typename> class PointAlloc,
        template<typename> class ConstraintAlloc
>
void DelaunayIndex <PointArg, TriangleArg, PointList,ConstraintList, PointAlloc, ConstraintAlloc>::setVertex (TriangleType &t, SideEnum s, IndexType v)
{
        IndexType current = Delaunay::getVertex (t, s);
        TrianglePtrVector &triangles = triangleIndex[current];
        triangles.erase (std::remove (triangles.begin (), triangles.end (), &t), triangles.end ());

        Delaunay::setVertex (t, s, v);
        TrianglePtrVector &trianglesV = triangleIndex[v];

        // insert zachowujący porządek sortowania.
        typename TrianglePtrVector::iterator i = std::lower_bound (trianglesV.begin (), trianglesV.end (), &t, TriangleCompare (v));
        trianglesV.insert (i, &t);
//        trianglesV.push_back (&t);
}

/****************************************************************************/

template <
        typename PointArg,
        typename TriangleArg,
        template<typename, typename> class PointList,
        template<typename, typename> class ConstraintList,
        template<typename> class PointAlloc,
        template<typename> class ConstraintAlloc
>
void DelaunayIndex <PointArg, TriangleArg, PointList,ConstraintList, PointAlloc, ConstraintAlloc>::clean ()
{
        triangulation.erase (std::remove_if (triangulation.begin (), triangulation.end (), TraingleRemovePredicate ()), triangulation.end ());
}

/****************************************************************************/

#ifndef NDEBUG
std::ostream &operator<< (std::ostream &o, std::vector <Triangle> const &e)
{
        typedef std::vector <Triangle> TriangleVector;

        size_t cnt = 0;
        for (typename TriangleVector::const_iterator i = e.begin (); i != e.end (); ++i, ++cnt) {
                o << cnt << ". " << *i << " A=";

//                if (i->tA) {
//                        o << *(i->tA);
//                }
//                else {
//                        o << "NULL";
//                }
//
//                o << " B=";
//
//                if (i->tB) {
//                        o << *(i->tB);
//                }
//                else {
//                        o << "NULL";
//                }
//
//                o << " C=";
//
//                if (i->tC) {
//                        o << *(i->tC);
//                }
//                else {
//                        o << "NULL";
//                }

                o << "\n";
        }

        return o;
}

std::ostream &operator<< (std::ostream &o, std::vector <Triangle *> const &e)
{
        size_t cnt = 0;
        for (std::vector <Triangle *>::const_iterator i = e.begin (); i != e.end (); ++i, ++cnt) {
                o << **i;

                if (i + 1 != e.end ()) {
                        o << " | ";
                }
        }

        return o;
}

std::ostream &operator<< (std::ostream &o, std::vector <std::vector <Triangle const *> > const &e)
{
        typedef std::vector <Triangle const *> TrianglePtrVector;
        typedef std::vector <TrianglePtrVector> TriangleIndex;

        size_t cnt = 0;
        for (TriangleIndex::const_iterator i = e.begin (); i != e.end (); ++i, ++cnt) {
                TrianglePtrVector const &trianglesForPoint = *i;

                for (TrianglePtrVector::const_iterator k = trianglesForPoint.begin (); k != trianglesForPoint.end (); ++k) {
                        o << cnt << " : " << **k;

                        if (k + 1 != trianglesForPoint.end ()) {
                                o << " | ";
                        }
                }

                o << "\n";
        }

        return o;
}
#endif

} // namespace

#endif /* DELAUNAYINDEX_H_ */
