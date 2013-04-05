/****************************************************************************
 *                                                                          *
 *  Author : lukasz.iwaszkiewicz@gmail.com                                  *
 *  ~~~~~~~~                                                                *
 *  License : see COPYING file for details.                                 *
 *  ~~~~~~~~~                                                               *
 ****************************************************************************/

#ifndef MY_POINT_H_
#define MY_POINT_H_

#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>
#include <vector>

struct Point {
        float X;
        float Y;
};

BOOST_STATIC_ASSERT (boost::has_trivial_assign <Point>::value);
BOOST_STATIC_ASSERT (boost::has_trivial_copy <Point>::value);
BOOST_STATIC_ASSERT (boost::is_pod <Point>::value);

typedef std::vector <Point> PointVector;

#endif /* POINT_H_ */
