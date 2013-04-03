/****************************************************************************
 *                                                                          *
 *  Author : lukasz.iwaszkiewicz@gmail.com                                  *
 *  ~~~~~~~~                                                                *
 *  License : see COPYING file for details.                                 *
 *  ~~~~~~~~~                                                               *
 ****************************************************************************/

#ifndef DELAUNAYEDGE_H_
#define DELAUNAYEDGE_H_

#include "Point.h"

namespace Delaunay {

/**
 *
 */
template <typename P>
struct Edge {
        typedef P PointType;
        Edge (PointType const a_, PointType const b_) : a (a_), b (b_) {}
        PointType a, b;
};

template <typename P>
float getArea (P const &a, P const &b, P const &c)
{
        return (getX (b) - getX (a)) * (getY (c) - getY (a)) - (getX (c) - getX (a)) * (getY (b) - getY (a));
}

template <typename P>
bool isLeft (P const &a, P const &b, P const &c)
{
        return getArea (a, b, c) > 0;
}

template <typename P>
bool isCollinear (P const &a, P const &b, P const &c)
{
        // TODO To musi pracować na INTACH z voronoia! Nie może być tego epsilona!
        return fabs (getArea (a, b, c)) < 0.00001;
}

/**
 * Computational geometry in C - second edition.
 */
template <typename P>
bool intersects (Edge <P> const &a, Edge <P> const &b)
{
        if (isCollinear (a.a, a.b, b.a) ||
            isCollinear (a.a, a.b, b.b) ||
            isCollinear (b.a, b.b, a.a) ||
            isCollinear (b.a, b.b, a.b)) {
                return false;
        }

        return (isLeft (a.a, a.b, b.a) ^ isLeft (a.a, a.b, b.b)) && (isLeft (b.a, b.b, a.a) ^ isLeft (b.a, b.b, a.b));
}

} // namespace

#endif /* DELAUNAYEDGE_H_ */
