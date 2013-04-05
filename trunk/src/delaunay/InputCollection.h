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
#include <algorithm>
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

        InputCollection () : /*points (0),*/ size_ (0) {}
        typedef TypeTraits <PointArg, TriangleArg, PointList> Traits;
        typedef typename Traits::PointType PointType;
        typedef typename Traits::TriangleType TriangleType;
        typedef typename Traits::PointListType PointListType;
        typedef typename Traits::PointListCollectionType PointListCollectionType;

        void setPoints (PointListType const &p);
        PointListType const &getPoints () const { return *pointsCollection.front (); }

        void addConstraint (PointListType const &p);
        PointListType const &getConstraint (size_t i) const { return *pointsCollection[i]; }
        size_t getConstraintOffset (size_t i) const;

        PointListCollectionType const &getPointsCollection () const { return pointsCollection; }

        PointType const &operator[] (size_t i) const;
        size_t size () const { return size_; }

private:

//        PointListType const *points;
        PointListCollectionType pointsCollection;
        size_t size_;
        typedef std::vector <size_t> BoundsVector;
        typedef BoundsVector::const_iterator BoundsIterator;
        BoundsVector bounds;

};

/****************************************************************************/

template <
        typename PointArg,
        typename TriangleArg,
        typename PointList
>
void InputCollection <PointArg, TriangleArg, PointList>::setPoints (PointListType const &p)
{
        assert (pointsCollection.empty ());
        pointsCollection.push_back (&p);
        size_ = p.size ();
//        bounds.push_back (size_ - 1);
}

/****************************************************************************/

template <
        typename PointArg,
        typename TriangleArg,
        typename PointList
>
void InputCollection <PointArg, TriangleArg, PointList>::addConstraint (PointListType const &p)
{
        assert (!pointsCollection.empty ());
        pointsCollection.push_back (&p);
        bounds.push_back (size_ - 1);
        size_ += p.size ();
}

/****************************************************************************/

template <
        typename PointArg,
        typename TriangleArg,
        typename PointList
>
typename InputCollection <PointArg, TriangleArg, PointList>::PointType const &
InputCollection <PointArg, TriangleArg, PointList>::operator[] (size_t i) const
{
        PointListType const *points = 0;
        size_t j = i;

        if (pointsCollection.size () == 1) {
                points = pointsCollection.front ();
        }
        else {

                // pierwszy który nie jest mniejszy (czyli równy większy).
                BoundsIterator b = std::lower_bound (bounds.begin (), bounds.end (), i);

                if (b == bounds.begin ()) {
                        points = pointsCollection.front ();
                }
                else if (b == bounds.end ()) {
                        points = pointsCollection.back ();
                        j -= bounds.back () + 1; // One before the last size.
                }
                else {
                        points = pointsCollection[(b - bounds.begin ()) - 1];
                        return pointsCollection[(b - bounds.begin ()) - 1]->operator[] (i - 10);
                }
        }

        return points->operator[] (j);
}

/****************************************************************************/

template <
        typename PointArg,
        typename TriangleArg,
        typename PointList
>
size_t InputCollection <PointArg, TriangleArg, PointList>::getConstraintOffset (size_t i) const
{
        if (!i) {
                return 0;
        }

        return bounds[i - 1] + 1;
}

}

#endif /* INPUTCOLLECTION_H_ */
