#ifndef UTILSVEC_H
#define UTILSVEC_H


namespace mp {

/// Grow vector capacity by a factor if needed.
/// @return reference to the element at \a i.
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
