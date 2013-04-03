/****************************************************************************
 *                                                                          *
 *  Author : lukasz.iwaszkiewicz@gmail.com                                  *
 *  ~~~~~~~~                                                                *
 *  License : see COPYING file for details.                                 *
 *  ~~~~~~~~~                                                               *
 ****************************************************************************/

#ifndef INPUTCOLLECTION_H_
#define INPUTCOLLECTION_H_

#include <vector>
#include "TypeTraits.h"

namespace Delaunay {

/**
 * Nie zadziała w przypadku, kiedy ktoś w miedzy czasie doda punkty do input lub do constraint.
 */
template <
        typename PointArg,
        typename TriangleArg,
        typename PointList
>
class InputCollection {
public:

        InputCollection () : points (0), size_ (0) {}
        typedef TypeTraits <PointArg, TriangleArg, PointList> Traits;
        typedef typename Traits::PointType PointType;
        typedef typename Traits::TriangleType TriangleType;
        typedef typename Traits::PointListType PointListType;
        typedef typename Traits::ConstraintListType ConstraintListType;

        void setPoints (PointListType const &p);
        PointListType const &getPoints () const { return *points; }

        void addConstraint (PointListType const &p);
        ConstraintListType const &getConstraints () const { return constraints; }

        PointType &operator[] (size_t i) { return points->operator[] (i); }
        PointType const &operator[] (size_t i) const { return points->operator[] (i); }

        size_t size () const { return size_; }

private:

        PointListType const *points;
        ConstraintListType constraints;
        size_t size_;
        std::vector <int> bounds;

};

template <
        typename PointArg,
        typename TriangleArg,
        typename PointList
>
void InputCollection <PointArg, TriangleArg, PointList>::setPoints (PointListType const &p)
{
        points = &p;
        size_ = points->size ();
        bounds.push_back (0);
}

template <
        typename PointArg,
        typename TriangleArg,
        typename PointList
>
void InputCollection <PointArg, TriangleArg, PointList>::addConstraint (PointListType const &p)
{
        constraints.push_back (&p);
        size_ += p.size ();
        bounds.push_back ();
}

}


#endif /* INPUTCOLLECTION_H_ */
