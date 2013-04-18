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

        InputCollection () : size_ (0) {}
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
        size_t getConstraintForIndex (size_t index) const;

        PointListCollectionType const &getPointsCollection () const { return pointsCollection; }

        PointType const &operator[] (size_t i) const;
        size_t size () const { return size_; }

private:

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

/****************************************************************************/

template <
        typename PointArg,
        typename TriangleArg,
        typename PointList
>
size_t InputCollection <PointArg, TriangleArg, PointList>::getConstraintForIndex (size_t i) const
{
        BoundsIterator b = std::lower_bound (bounds.begin (), bounds.end (), i);

        if (b == bounds.begin ()) {
                return 0;
        }
        else if (b == bounds.end ()) {
                return pointsCollection.size () - 1;
        }
        else {
                return (b - bounds.begin ()) - 1;
        }
}

/****************************************************************************/

template <
        typename PointArg,
        typename TriangleArg,
        typename PointList
>
void InputCollection <PointArg, TriangleArg, PointList>::establishSuperTriangle (PointListType const &p)
{
        std::vector <CoordinateType> coordinate;
        coordinate.reserve (p.size ());

        for (typename PointListType::const_iterator i = p.begin (), e = p.end (); i != e; ++i) {
                coordinate.push_back (getX (*i));
        }

        std::sort (coordinate.begin (), coordinate.end ());

        CoordinateType xmin = coordinate.front ();
        CoordinateType xmax = coordinate.back ();
        CoordinateType dx = (xmax - xmin) / 4.0;
        xmin -= dx;
        xmax += dx;

        int j = 0;
        for (typename PointListType::const_iterator i = p.begin (), e = p.end (); i != e; ++i, ++j) {
                coordinate[j] = getY (*i);
        }

        std::sort (coordinate.begin (), coordinate.end ());

        CoordinateType ymin = coordinate.front ();
        CoordinateType ymax = coordinate.back ();
        CoordinateType dy = (ymax - ymin) / 4.0;
        ymin -= dy;
        ymax += dy;

#if 0
        std::cerr << xmin << ", " << xmax << ", " << ymin << ", " << ymax << std::endl;
#endif

        PointType pt;
        setX (pt, xmin);
        setY (pt, ymin);
        superTriangle.push_back (pt);

        setX (pt, xmax);
        setY (pt, ymin);
        superTriangle.push_back (pt);

        setX (pt, xmax);
        setY (pt, ymax);
        superTriangle.push_back (pt);

        setX (pt, xmin);
        setY (pt, ymax);
        superTriangle.push_back (pt);

        pointsCollection.push_back (&superTriangle);
        size_ += superTriangle.size ();
}

} // namespace

#endif /* INPUTCOLLECTION_H_ */
