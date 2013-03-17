/****************************************************************************
 *                                                                          *
 *  Author : lukasz.iwaszkiewicz@gmail.com                                  *
 *  ~~~~~~~~                                                                *
 *  License : see COPYING file for details.                                 *
 *  ~~~~~~~~~                                                               *
 ****************************************************************************/

#ifndef DELAUNAYTRIANGLE_H_
#define DELAUNAYTRIANGLE_H_

#ifndef NDEBUG
#include <iosfwd>
#endif

namespace Delaunay {

/**
 * Numbers triangle edges and vertices. Edge A (first edge) id at the opposite
 * side than vertex A.
 */
enum SideEnum { N = 0, A = 1, B = 2, C = 3 };

/**
 * Default triangle type (index based).
 */
struct Triangle {

        typedef unsigned int IndexType;

        Triangle (IndexType x = 0, IndexType y = 0, IndexType z = 0) : a (x), b (y), c (z) {}

        IndexType get (SideEnum s) const
        {
                switch (s) {
                case A:
                        return a;
                case B:
                        return b;
                case C:
                        return c;
                default:
                        return 0;
                }
        }

        void set (SideEnum s, IndexType v)
        {
                switch (s) {
                case A:
                        a = v;
                        break;
                case B:
                        b = v;
                        break;
                case C:
                        c = v;
                        break;
                default:
                        break;
                }
        }

        // Endpoints of a triangle.
        IndexType a, b, c;
};

/**
 *
 */
template<typename T>
struct TriangleTraits {
        typedef typename T::IndexType IndexType;

        /**
         * Gets one of triangle vertices.
         */
        static inline IndexType get (T const &triangle, SideEnum side)
        {
                return triangle.get (side);
        }
};

/**
 *
 */
template<typename T>
struct TriangleMutableTraits {
        typedef typename TriangleTraits <T>::IndexType IndexType;

        /**
         * Sets one of triangle vertices.
         */
        static inline void set (T &triangle, SideEnum side, IndexType value)
        {
                return triangle.set (side, value);
        }

        /**
         * Zero initialized triangle.
         */
        static inline T construct (IndexType a = 0, IndexType b = 0, IndexType c = 0)
        {
                return T (a, b, c);
        }
};

template <typename T>
typename TriangleTraits <T>::IndexType a (T const &triangle)
{
        return TriangleTraits <T>::get (triangle, A);
}

template <typename T>
void a (T &triangle, typename TriangleTraits <T>::IndexType i)
{
        TriangleMutableTraits <T>::set (triangle, A, i);
}

template <typename T>
typename TriangleTraits <T>::IndexType b (T const &triangle)
{
        return TriangleTraits <T>::get (triangle, B);
}

template <typename T>
void b (T &triangle, typename TriangleTraits <T>::IndexType i)
{
        TriangleMutableTraits <T>::set (triangle, B, i);
}

template <typename T>
typename TriangleTraits <T>::IndexType c (T const &triangle)
{
        return TriangleTraits <T>::get (triangle, C);
}

template <typename T>
void c (T &triangle, typename TriangleTraits <T>::IndexType i)
{
        TriangleMutableTraits <T>::set (triangle, C, i);
}

template <typename T>
bool hasVertex (T const &t, typename TriangleTraits <T>::IndexType p)
{
        return p == t.a || p == t.b || p == t.c;
}

template <typename T>
typename TriangleTraits <T>::IndexType getVertex (T const &t, SideEnum side)
{
        return TriangleTraits <T>::get (t, side);
}

template <typename T>
void setVertex (T &t, SideEnum side, typename TriangleTraits <T>::IndexType i)
{
        TriangleMutableTraits <T>::set (t, side, i);
}

template <typename T>
SideEnum getVertexSide (T const &t, typename TriangleTraits <T>::IndexType i)
{
        return static_cast <SideEnum> ((i == a (t)) | ((i == b (t)) << 1) | (i == c (t)) | ((i == c (t)) << 1));
}

inline
SideEnum otherThan (SideEnum a, SideEnum b)
{
      return static_cast <SideEnum> (6 - static_cast <int> (a) - static_cast <int> (b));
}

inline
std::pair <SideEnum, SideEnum> otherThan (SideEnum s)
{
        switch (s) {
        case A:
                return std::make_pair (B, C);
        case B:
                return std::make_pair (A, C);
        case C:
                return std::make_pair (A, B);
        default:
                return std::make_pair (N, N);
        }
}

#ifndef NDEBUG
std::ostream &operator<< (std::ostream &o, Triangle const &t)
{
        o << a (t) << "," << b (t) << "," << c (t);
        return o;
}
#endif

} // namespace Delaunay

#endif /* DELAUNAYTRIANGLE_H_ */
