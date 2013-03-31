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
class InputCollection {
public:

        typedef PointArg PointType;
        typedef TriangleArg TriangleType;
        typedef PointList <PointArg, PointAlloc <PointArg> > PointListType;
        typedef ConstraintList <PointListType const *, ConstraintAlloc <PointListType const *> > ConstraintListType;

        void setPoints (PointListType const &p) { points = &p; }
        PointListType const &getPoints () const { return *points; }

        void addConstraint (PointListType const &p) { constraints.push_back (&p); }
        ConstraintListType const &getConstraints () const { return constraints; }

        PointType &operator[] (size_t i);
        PointType const &operator[] (size_t i) const;

private:

        PointListType const *points;
        ConstraintListType constraints;

};

}


#endif /* INPUTCOLLECTION_H_ */
