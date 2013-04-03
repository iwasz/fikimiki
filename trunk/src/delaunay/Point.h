/****************************************************************************
 *                                                                          *
 *  Author : lukasz.iwaszkiewicz@gmail.com                                  *
 *  ~~~~~~~~                                                                *
 *  License : see COPYING file for details.                                 *
 *  ~~~~~~~~~                                                               *
 ****************************************************************************/

#ifndef DELAUNAYPOINT_H_
#define DELAUNAYPOINT_H_

#include <boost/polygon/isotropy.hpp>
#include <boost/polygon/point_concept.hpp>
#include <boost/math/special_functions/round.hpp>

namespace Delaunay {

/**
 *
 */
enum PointCoordinate { X = 1, Y = 2 };

/**
 * Defaul Point class. It is assumed that user has his own implementation, thus this class
 * is provided only for saske of completeness.
 */
struct Point {
        typedef float CoordinateType;

        Point (CoordinateType x_ = 0, CoordinateType y_ = 0) : x (x_), y (y_)  {}
        CoordinateType x, y;

        inline CoordinateType get (PointCoordinate s) const
        {
                return ((s == X) ? (x) : (y));
        }

        inline void set (PointCoordinate s, CoordinateType v)
        {
                if (s == X) {
                        x = v;
                }
                else {
                        y = v;
                }
        }
};


template<typename T>
struct PointTraits {
        typedef typename T::CoordinateType CoordinateType;
        typedef int IntCoordinateType;

        static inline CoordinateType get (const T& point, PointCoordinate orient)
        {
                return point.get (orient);
        }
};

template<typename T>
struct PointMutableTraits {
        typedef typename PointTraits <T>::CoordinateType CoordinateType;

        static inline void set (T& point, PointCoordinate orient, CoordinateType c)
        {
                return point.set (orient, c);
        }


        /**
         * Zero initialized point.
         */
        static inline T construct (CoordinateType x = 0, CoordinateType y = 0)
        {
                return T (x, y);
        }
};

/**
 * Helper accessor.
 */
template <typename T>
typename PointTraits <T>::CoordinateType getX (T const &point)
{
        return PointTraits <T>::get (point, X);
}

/**
 * Helper accessor.
 */
template <typename T>
void setX (T const &point, typename PointTraits <T>::CoordinateType p)
{
        PointMutableTraits <T>::set (point, X, p);
}

/**
 * Helper accessor.
 */
template <typename T>
typename PointTraits <T>::CoordinateType getY (T const &point)
{
        return PointTraits <T>::get (point, Y);
}

/**
 * Helper accessor.
 */
template <typename T>
void setY (T const &point, typename PointTraits <T>::CoordinateType p)
{
        PointMutableTraits <T>::set (point, Y, p);
}

//} // namespace Helper
} // namespace Delaunay

namespace boost {
namespace polygon {

template <>
struct geometry_concept<Delaunay::Point> {
        typedef point_concept type;
};

template<>
struct point_traits<Delaunay::Point> {
        typedef Delaunay::PointTraits <Delaunay::Point> DelaunayPointTraitsType;
        typedef DelaunayPointTraitsType::IntCoordinateType coordinate_type;

        /*
         * (ok 2*n) Sprawdzić jak często to się wykonuje i skąd. Jak tylko w kroku inicjacji, to OK.
         * Dał bym tutaj też mnożnik, żeby uniknąć sytuacji, kiedy dwa zmiennoprzecinkowe punkty
         * wejściowe, które są bardzo blisko siebie zostaną przez poniższy get zwróceone jako
         * ten sam punkt. Można albo mnożyć przez stała (np 1000), albo znaleźć najmniejszą różnicę
         * mięczy dwoma współrzednymi w danych wejściowych i przeskalować je odpowiednio.
         */
        static inline coordinate_type get (const Delaunay::Point& point, orientation_2d orient)
        {
                return (orient == HORIZONTAL) ? boost::math::iround (Delaunay::getX (point)) : boost::math::iround (Delaunay::getY (point));
        }
};

}
}

#endif /* DELAUNAYPOINT_H_ */
