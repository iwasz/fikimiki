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

namespace Delaunay {

template <
        typename PointArg,
        typename TriangleArg,
        template<typename, typename> class PointList = std::vector,
        template<typename, typename> class ConstraintList = std::vector,
        template<typename> class PointAlloc = std::allocator,
        template<typename> class ConstraintAlloc = std::allocator
>
struct InputCollection {

        typedef PointArg PointType;
        typedef TriangleArg TriangleType;
        typedef PointList <PointArg, PointAlloc <PointArg> > PointListType;
        typedef ConstraintList <PointListType, ConstraintAlloc <PointListType> > ConstraintListType;

        void setPoints (PointListType const &p);
        PointListType &getPoints ();

        void addConstraint (PointListType const &p);
        ConstraintListType &getConstraints ();

        PointType &operator[] (size_t i);
        PointType const &operator[] (size_t i) const;
};

}


#endif /* INPUTCOLLECTION_H_ */
