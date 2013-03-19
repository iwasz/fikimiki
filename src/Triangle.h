/****************************************************************************
 *                                                                          *
 *  Author : lukasz.iwaszkiewicz@gmail.com                                  *
 *  ~~~~~~~~                                                                *
 *  License : see COPYING file for details.                                 *
 *  ~~~~~~~~~                                                               *
 ****************************************************************************/

#ifndef MY_TRIANGLE_H_
#define MY_TRIANGLE_H_

#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>

struct Triangle {

        unsigned int a;
        unsigned int b;
        unsigned int c;
};

BOOST_STATIC_ASSERT (boost::has_trivial_assign <Point>::value);
BOOST_STATIC_ASSERT (boost::has_trivial_copy <Point>::value);
BOOST_STATIC_ASSERT (boost::is_pod <Point>::value);

#endif /* TRIANGLE_H_ */
