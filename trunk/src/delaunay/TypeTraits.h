/****************************************************************************
 *                                                                          *
 *  Author : lukasz.iwaszkiewicz@gmail.com                                  *
 *  ~~~~~~~~                                                                *
 *  License : see COPYING file for details.                                 *
 *  ~~~~~~~~~                                                               *
 ****************************************************************************/

#ifndef FIKIMIKI_TYPETRAITS_H_
#define FIKIMIKI_TYPETRAITS_H_

#include <vector>
#include <list>
#include <boost/tuple/tuple.hpp>

namespace Delaunay {

template <
        typename PointArg = Point,
        typename TriangleArg = Triangle,
        typename PointList = std::vector <PointArg>
>
struct TypeTraits {

        typedef PointArg PointType;
        typedef TriangleArg TriangleType;
        typedef PointList PointListType;
        typedef std::vector <PointListType const *> ConstraintListType;
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

};

} // namespace

#endif /* TYPETRAITS_H_ */
