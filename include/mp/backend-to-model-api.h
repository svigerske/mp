#ifndef BACKENDTOMODELAPI_H
#define BACKENDTOMODELAPI_H

#include <cassert>

namespace mp {

/// A wrapper for CommonInfo providing some functionality
/// to connect Backend to ModelAPI.
/// @param CommonInfo: common information
///   (e.g., env/model handlers) stored in both
///   Backend and ModelAPI
template <class CommonInfo>
class Backend2ModelAPIConnector :
    public CommonInfo {
public:
  // Retrieve the other object
  CommonInfo* get_other() { assert(other_); return other_; }
  /// Set the other object.
  /// This is standardly called from CreateModelManagerWithModelAPI()
  void set_other(CommonInfo* other) { other_ = other; }

  /// Copy shared info from the other object.
  /// Typically after initializing solver env/model
  void copy_common_info_from_other() { (CommonInfo&)(*this) = *get_other(); }
  /// Copy into.
  /// Typically after initializing solver env/model
  void copy_common_info_to_other() { *get_other() = (CommonInfo&)(*this); }

protected:
  /// Call this in Backend constructor to create
  /// a ModelManager with ModelAPI.
  /// For recompilation speed, define the template
  /// in <solver>modelapi.cc.
  template <class MMInfo>
  static void CreateModelManager(MMInfo& );

private:
  CommonInfo* other_ = nullptr;
};


} // namespace mp

#endif // BACKENDTOMODELAPI_H
