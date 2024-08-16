#ifndef UTILSVEC_H
#define UTILSVEC_H

// #include "../../thirdparty/martinus/svector/svector.h"
#include "../../thirdparty/gharveymn/small_vector/small_vector.hpp"


namespace mp {

/// Typedef small vector
// template <class T, size_t N>
// using SmallVec = ankerl::svector<T, N>;

/// Typedef small vector
template <class T, size_t N>
using SmallVec = gch::small_vector<T, N>;

/// Typedef small vector, default size of 64 bytes
template <class T>
using SmallVecDefSz = gch::small_vector<T>;

/// Grow vector capacity by a factor if needed.
/// @return reference to the element at \a i.
/// @note better preallocate, or call in the reverse order of indexes.
template <class Vec>
auto AutoExpand(Vec& vec, typename Vec::size_type i) {
  if (vec.size()<=i) {
    if (vec.capacity()<=i)
      vec.reserve(((i+1)*13)/10);
    vec.resize(i+1);
  }
  return vec[i];
}

}  // namespace mp


#endif // UTILSVEC_H
