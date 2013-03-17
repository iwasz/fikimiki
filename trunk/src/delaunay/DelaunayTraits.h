/****************************************************************************
 *                                                                          *
 *  Author : lukasz.iwaszkiewicz@gmail.com                                  *
 *  ~~~~~~~~                                                                *
 *  License : see COPYING file for details.                                 *
 *  ~~~~~~~~~                                                               *
 ****************************************************************************/

#ifndef DELAUNAYTRIANGULATIONTRAITS_H_
#define DELAUNAYTRIANGULATIONTRAITS_H_

#include "DelaunayTriangle.h"
#include "DelaunayPoint.h"

namespace Delaunay {

/**
 *
 */
template <typename P = Point, typename T = Triangle>
struct DelaunayTriangulationTraits {
        typedef P PointType;
        typedef T TriangleType;

        typedef TriangleTraits <TriangleType> TriangleTraitsType;
        typedef typename TriangleTraitsType::IndexType IndexType;

};

} // namespace Delaunay

#endif /* DELAUNAYTRIANGULATIONTRAITS_H_ */
