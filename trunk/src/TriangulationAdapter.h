/****************************************************************************
 *                                                                          *
 *  Author : lukasz.iwaszkiewicz@gmail.com                                  *
 *  ~~~~~~~~                                                                *
 *  License : see COPYING file for details.                                 *
 *  ~~~~~~~~~                                                               *
 ****************************************************************************/

#ifndef TRIANGULATIONADAPTER_H_
#define TRIANGULATIONADAPTER_H_

#include "delaunay/Triangulation.h"
#include "Point.h"
#include "Triangle.h"

namespace boost {
namespace polygon {

template <>
struct geometry_concept<Point> {
        typedef point_concept type;
};

template<>
struct point_traits<Point> {
//        typedef Delaunay::PointTraits <Delaunay::Point> DelaunayPointTraitsType;
        typedef int coordinate_type;

        static inline coordinate_type get (const Point& point, orientation_2d orient)
        {
                // TODO This 1000 factor is arbitrary - consider other options.
                return (orient == HORIZONTAL) ? boost::math::iround (point.X * 1000) : boost::math::iround (point.Y * 1000);
        }
};

}
}

namespace Delaunay {

template <>
struct TriangleTraits <typename ::Triangle> {
        typedef unsigned int IndexType;

        /**
         * Gets one of triangle vertices.
         */
        static inline IndexType get (typename ::Triangle const &triangle, Delaunay::SideEnum side)
        {
                switch (side) {
                case Delaunay::A:
                        return triangle.a;
                case Delaunay::B:
                        return triangle.b;
                case Delaunay::C:
                        return triangle.c;
                default:
                        return 0;
                }
        }
};

template <>
struct TriangleMutableTraits <typename ::Triangle> {
        typedef typename TriangleTraits <typename ::Triangle>::IndexType IndexType;

        /**
         * Sets one of triangle vertices.
         */
        static inline void set (typename ::Triangle &triangle, Delaunay::SideEnum side, IndexType value)
        {
                switch (side) {
                case Delaunay::A:
                        triangle.a = value;
                        break;
                case Delaunay::B:
                        triangle.b = value;
                        break;
                case Delaunay::C:
                        triangle.c = value;
                        break;
                default:
                        break;
                }
        }

        /**
         * Zero initialized triangle.
         */
        static inline typename ::Triangle construct (IndexType a = 0, IndexType b = 0, IndexType c = 0)
        {
                typename ::Triangle t;
                t.a = a;
                t.b = b;
                t.c = c;
                return t;
        }
};


template<>
struct PointTraits <typename ::Point> {
        typedef float CoordinateType;
        typedef int IntCoordinateType;

        static inline CoordinateType get (const typename ::Point& point, PointCoordinate orient)
        {
                return ((orient == X) ? (point.X) : (point.Y));
        }
};

template<>
struct PointMutableTraits <typename ::Point> {
        typedef typename PointTraits <typename ::Point>::CoordinateType CoordinateType;

        static inline void set (typename ::Point& point, PointCoordinate orient, CoordinateType c)
        {
                if (orient == X) {
                        point.X = c;
                }
                else {
                        point.Y = c;
                }
        }


        /**
         * Zero initialized point.
         */
        static inline typename ::Point construct (CoordinateType x = 0, CoordinateType y = 0)
        {
                typename ::Point p;
                p.X = x;
                p.Y = y;
                return p;
        }
};

} // namespace Delaunay

typedef Delaunay::Triangulation <Point, Triangle> MyTriagulation;
typedef MyTriagulation::TriangleVector TriangleVector;

#endif /* TRIANGULATIONADAPTER_H_ */
