#ifndef VALUE_PRESOLVE_LINK_H
#define VALUE_PRESOLVE_LINK_H

#include <array>
#include <vector>
#include <deque>
#include <algorithm>
#include <mp/arrayref.h>

#include "mp/common.h"
#include "valcvt-node.h"


namespace mp {

namespace pre {

/// Index range in a link
using LinkIndexRange = IndexRange;

/// Declare ValuePresolver
class ValuePresolver;


/// ValuePresolveLink interface
///
/// A link is a node of the structural conversion graph.
/// It contains an array of concrete nodes ('link entries'),
/// each describes a transformation from
/// some source vars+cons+objs values into/from some target ones.
/// All concrete nodes are of the same type.
class BasicLink {
public:
  /// Constructor
  BasicLink(ValuePresolver& pre) : value_presolver_(pre) { }

  /// Virtual destructor
  virtual ~BasicLink() = default;

  /// Type name
  virtual const char* GetTypeName() const = 0;

  /// The below pre- / postsolve methods
  /// work on a range of link entries.
  /// Postsolves should usually loop the range backwards
#undef PRESOLVE_KIND
#define PRESOLVE_KIND(name, ValType) \
  virtual void Presolve ## name (LinkIndexRange, \
    ValueResolution ) = 0; \
  virtual void Postsolve ## name(LinkIndexRange, \
    ValueResolution ) = 0;

  LIST_PRESOLVE_METHODS


  /// Typedef item container
  using ItemRangeList = std::vector<NodeRange>;

  /// Container of source and target items
  /// for a certain link entry. Used for graph export
  struct EntryItems {
    ItemRangeList src_items_, dest_items_;
  };

  /// Get source/target nodes for a given link entry.
  /// This is used for graph export
  virtual void ExportEntryItems(EntryItems& ei, int i) const = 0;


protected:
  /// Add a range of link entries to the Presolver's list.
  /// Every derived link should call either of the next two methods
  /// whenever adding a link entry.
  /// Version 1: add single link
  void RegisterLinkIndex(int i)
  { RegisterLinkIndexRange( {i, i+1} ); }

  /// Version 2: add a range of links
  void RegisterLinkIndexRange(LinkIndexRange );


private:
  ValuePresolver& value_presolver_;
};


/// Link range: range of conversion specifiers of certain type.
/// The link object is specified as well
struct LinkRange {
  /// Try & extend the range.
  /// @return true iff extension worked,
  /// otherwise the caller would have to record the new range
  /// separately
  bool TryExtendBy(LinkRange br) {
    if (&b_ == &br.b_) {                  // same link?
      if (ir_.end_ == br.ir_.beg_) {      // and consecutive ranges?
        ir_.end_ = br.ir_.end_;
        return true;
      }
    }
    return false;
  }

  /// Reference to BasicLink
  BasicLink& b_;

  /// Link index range
  LinkIndexRange ir_;
};


/// A specific link: each entry just copies a range of values.
/// Useful to transfer values between NL and internal model,
/// internal model and solver, and in conversions preserving
/// all values.
/// Feature: copies values without comparing with existing values,
/// in contrast to other links using Set() method.
/// This is useful for objectives which use both CopyLink
/// and One2ManyLink's, to retrieve objective values
/// and not max them out on duals, etc.
class CopyLink : public BasicLink {
public:
  /// Constructor
  CopyLink(ValuePresolver& pre) : BasicLink(pre) { }

  /// Type name
  const char* GetTypeName() const override { return "CopyLink"; }

  /// Single link entry,
  /// stores src + dest ranges
  using LinkEntry = std::pair<NodeRange, NodeRange>;

  /// Collection of entries
  using CollectionOfEntries = std::deque<LinkEntry>;

  /// Add entry.
  /// Instead of a new entry, tries to extend the last one
  /// if exists
  void AddEntry(LinkEntry be) {
    if (entries_.empty() ||
        !entries_.back().first.ExtendableBy(be.first) ||
        !entries_.back().second.ExtendableBy(be.second)) {
      entries_.push_back(be);             // Add new entry
      RegisterLinkIndex(entries_.size()-1);
    } else {                              // Extend last entry
      entries_.back().first.ExtendBy(be.first);
      entries_.back().second.ExtendBy(be.second);
    }
  }

  /// Get source/target nodes for a given link entry.
  /// This is used for graph export
  void ExportEntryItems(EntryItems& ei, int i) const override {
    const auto& en = entries_.at(i);
    ei.src_items_.clear();
    ei.src_items_.push_back(en.first);
    ei.dest_items_.clear();
    ei.dest_items_.push_back(en.second);
  }

  /// Copy everything, MaxAmongNon0 should not apply
#undef PRESOLVE_KIND
#define PRESOLVE_KIND(name, ValType) \
  void Presolve ## name (LinkIndexRange ir, ValueResolution ) override \
  { CopySrcDest <ValType>(ir); } \
  void Postsolve ## name(LinkIndexRange ir, ValueResolution ) override \
  { CopyDestSrc <ValType>(ir); }

  LIST_PRESOLVE_METHODS

protected:
  /// Copy src -> dest for index range ir
  template <class T>
  void CopySrcDest(LinkIndexRange ir) {
    for (int i=ir.beg_; i!=ir.end_; ++i) {
      const auto& br = entries_[i];
      Copy<T>(br.first, br.second);
    }
  }
  /// Copy src <- dest for index range ir. Loop backwards
  template <class T>
  void CopyDestSrc(LinkIndexRange ir) {
    for (int i=ir.end_; (i--)!=ir.beg_; ) {
      const auto& br = entries_[i];
      Copy<T>(br.second, br.first);
    }
  }


private:
  CollectionOfEntries entries_;
};


/// A very general link:
/// each entry copies values from/to
/// *every* position of the source range
/// into/from *every* position of the target range.
/// Is specialized by One2Many and Many2One, and
/// is not recommended to be used directly.
/// In the cases when such a Many2Many connection
/// exists, consider splitting via a (dummy) intermediate node
/// and linking via Many2One + One2Many.
/// Example: Unary Encoding.
class Many2ManyLink : public BasicLink {
public:
  /// Constructor
  Many2ManyLink(ValuePresolver& pre) : BasicLink(pre) { }

  /// Type name
  const char* GetTypeName() const override { return "Many2ManyLink"; }

  /// Single link entry,
  /// stores src + dest ranges.
  using LinkEntry = std::pair<NodeRange, NodeRange>;

  /// Collection of entries
  using CollectionOfEntries = std::deque<LinkEntry>;

  /// Add entry.
  /// Instead of a new entry, tries to extend the last one
  /// if exists
  void AddEntry(LinkEntry be) {
    if (entries_.empty() ||
        !(
          (entries_.back().first==be.first   // same sources
           && entries_.back().second.TryExtendBy(be.second))
          ||
          (entries_.back().second==be.second // same targets
           && entries_.back().first.TryExtendBy(be.first))
          )) {
      entries_.push_back(be);             // Add new entry
      if_consider_for_last_link_.push_back(false);
      RegisterLinkIndex(entries_.size()-1);
    }
  }

  /// Mark entry as contining the last link of a group.
  /// needed for solution pre-postsolve where
  /// we actually need the last link of a group,
  /// corr. to the main constraint.
  /// Example: a nonlinear con is redefined into a bunch
  /// of smaller constraints, but only the main
  /// resulting con receives/sends the dual value.
  void MarkLastEntryAsLastInGroup() {
    assert(if_consider_for_last_link_.size() == entries_.size());
    assert(entries_.size());
    if_consider_for_last_link_.back() = true;
  }

  /// Get source/target nodes for a given link entry.
  /// This is used for graph export
  void ExportEntryItems(EntryItems& ei, int i) const override {
    const auto& en = entries_.at(i);
    ei.src_items_.clear();
    ei.src_items_.push_back(en.first);
    ei.dest_items_.clear();
    ei.dest_items_.push_back(en.second);
  }

  /// All pre- / postsolves just take max from non-0
#undef PRESOLVE_KIND
#define PRESOLVE_KIND(name, ValType) \
  void Presolve ## name (LinkIndexRange ir, ValueResolution vr) override \
  { DistributeFromSrc2Dest <ValType>(ir, vr); } \
  void Postsolve ## name(LinkIndexRange ir, ValueResolution vr) override \
  { CollectFromDest2Src <ValType>(ir, vr); }

  LIST_PRESOLVE_METHODS

protected:
  /// Distribute values of type T from nr1 to nr2
  template <class T>
  void Distr(NodeRange nr1, NodeRange nr2, ValueResolution vr) {
    auto ir1 = nr1.GetIndexRange();
    auto ir2 = nr2.GetIndexRange();
    for (auto i0=ir1.beg_; i0!=ir1.end_; ++i0) {
      // Need reference here for reference counting
      // in PresolveNames():
      const auto& val = nr1.GetValueNode()->
          GetVal<T>(i0);
      if (ValueResolution::ValResAll == vr) {
        for (auto i=ir2.beg_; i!=ir2.end_; ++i)
          nr2.GetValueNode()->SetVal(i, val);
      } else {
        assert(ValueResolution::ValResLast == vr);
        nr2.GetValueNode()->SetVal(ir2.end_ - 1, val);
      }
    }
  }

  /// Collect values of type T from nr2 to nr1
  template <class T>
  void Collect(NodeRange nr1, NodeRange nr2, ValueResolution vr) {
    auto ir1 = nr1.GetIndexRange();
    auto ir2 = nr2.GetIndexRange();
    auto& vec2 = nr2.GetValueNode()->GetValVec<T>();
    for (auto i0=ir1.beg_; i0!=ir1.end_; ++i0) {
      if (ValueResolution::ValResAll == vr) {
        for (auto i=ir2.beg_; i!=ir2.end_; ++i) {
          nr1.GetValueNode()->SetVal(i0, vec2.at(i));
        }
      } else {
        assert(ValueResolution::ValResLast == vr);
        nr1.GetValueNode()->SetVal(i0, vec2.at(ir2.end_ - 1));
      }
    }
  }

  /// Src -> dest for the entries index range ir
  template <class T>
  void DistributeFromSrc2Dest(LinkIndexRange ir, ValueResolution vr) {
    for (int i=ir.beg_; i!=ir.end_; ++i) {
      if (ValueResolution::ValResAll == vr
          || if_consider_for_last_link_[i]) {
        const auto& br = entries_[i];
        Distr<T>(br.first, br.second, vr);
      }
    }
  }

  /// Collect src <- dest for index range ir. Loop backwards
  template <class T>
  void CollectFromDest2Src(LinkIndexRange ir, ValueResolution vr) {
    for (int i=ir.end_; (i--)!=ir.beg_; ) {
      if (ValueResolution::ValResAll == vr
          || if_consider_for_last_link_[i]) {
        const auto& br = entries_[i];
        Collect<T>(br.first, br.second, vr);
      }
    }
  }


private:
  CollectionOfEntries entries_;
  std::vector<bool> if_consider_for_last_link_;
};


/// A specialization of Many2ManyLink:
/// each entry just copies from 1 source value
/// into a range of values.
/// Useful to transfer values into expression subtree,
/// in conversions, e.g., from 1 var / constraint / objective
/// into several new ones
class One2ManyLink : public Many2ManyLink {
public:
  /// Constructor
  One2ManyLink(ValuePresolver& pre) : Many2ManyLink(pre) { }

  /// Type name
  const char* GetTypeName() const override { return "One2ManyLink"; }

  /// Add entry.
  void AddEntry(LinkEntry be) {
    assert(be.first.IsSingleIndex());
    Many2ManyLink::AddEntry(be);
  }
};


/// A specialization of Many2ManyLink:
/// each entry just copies from a range of values
/// into 1 target value.
/// Useful for globalizing conversions,
/// to transfer values from several items into a new item
class Many2OneLink : public Many2ManyLink {
public:
  /// Constructor
  Many2OneLink(ValuePresolver& pre) : Many2ManyLink(pre) { }

  /// Type name
  const char* GetTypeName() const override { return "Many2OneLink"; }

  /// Add entry.
  void AddEntry(LinkEntry be) {
    assert(be.second.IsSingleIndex());
    Many2ManyLink::AddEntry(be);
  }
};


/// A stub for links which process each entry individually,
/// in contrast to Presolve...(LinkIndexRange ...).
/// Some of such links could be optimized to store
/// a range of transformations in each entry if that helps
/// (but this could be an overoptimization?)
///
/// Usage: a derived class should define some methods
/// Presolve...(const LinkEntry& ) and Postsolve...(const LinkEntry&).
///
/// We could allow default methods as follows:
/// Those methods which are not defined, just copy values
/// (which might be correct in _some_ cases).
/// Then, need a "default copy" method.
///
/// Using CRTP:
/// https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
template <class Impl, class LinkEntry>
class BasicIndivEntryLink : public BasicLink {
public:
  /// Constructor
  BasicIndivEntryLink(ValuePresolver& pre) :
    BasicLink(pre) { }

  /// Add entry
  void AddEntry(LinkEntry be) {
    entries_.push_back(be);             // Add physically
    RegisterLinkIndex(entries_.size()-1);
  }

  /// Get source/target nodes for a given link entry.
  /// This is used for graph export
  void ExportEntryItems(EntryItems& ei, int i) const override {
    MPCD( FillEntryItems( ei, entries_.at(i) ) );
  }


  /// Pre- / postsolve loops over link entries
  /// and calls the derived class' method for each.
#undef PRESOLVE_KIND
#define PRESOLVE_KIND(name, ValType) \
  void Presolve ## name (LinkIndexRange ir, ValueResolution) override { \
    for (int i=ir.beg_; i!=ir.end_; ++i) \
      MPD( Presolve ## name ## Entry(entries_.at(i)) ); } \
  void Postsolve ## name(LinkIndexRange ir, ValueResolution) override { \
    for (int i=ir.end_; i--!=ir.beg_; ) \
      MPD( Postsolve ## name ## Entry(entries_.at(i)) ); }

  LIST_PRESOLVE_METHODS


private:
  std::deque<LinkEntry> entries_;
};


/// A static indiv entry link has a fixed number (\a NNodes)
/// of ValueNodes and indexes (\a NIndexes) into them.
/// Generally \a NNodes==\a NIndexes.
/// \a NSources < \a NIndexes is the number of source nodes/indexes
/// (coming first in each LinkEntry, this is used for export).
template <class Impl, int NNodes, int NIndexes, int NSources>
class BasicStaticIndivEntryLink :
    public BasicIndivEntryLink<Impl, std::array<int, NIndexes> > {
public:
  /// Base class
  using Base = BasicIndivEntryLink<Impl, std::array<int, NIndexes> >;

  /// Typedef NodeList: list of ValueNode pointers.
  /// A ValueNode is a node of the (structural) transformation graph
  /// and can store arrays of values for a list of model items
  using NodeList = std::array<ValueNode*, NNodes>;

  /// Typedef: LinkEntry is just an array if node indexes
  using LinkEntry = std::array<int, NIndexes>;

  /// Constructor
  BasicStaticIndivEntryLink(ValuePresolver& pre, const NodeList& ndl) :
    Base(pre), ndl_(ndl) { }

  /// Get source/target nodes for a given link entry.
  /// This is used for graph export
  void FillEntryItems(typename Base::EntryItems& ei,
                      const LinkEntry& en) const {
    ei.src_items_.clear();
    int i=0;
    for ( ; i<NSources; ++i)
      ei.src_items_.push_back({ndl_.at(i), en.at(i)});
    ei.dest_items_.clear();
    for ( ; i<NIndexes; ++i)
      ei.dest_items_.push_back({ndl_.at(i), en.at(i)});
  }


protected:
  /// Access whole node at specific index, const
  const ValueNode& GetNode(size_t i) const { return *ndl_.at(i); }
  /// Access whole node at specific index
  ValueNode& GetNode(size_t i) { return *ndl_.at(i); }

  /// Access int value at the node \a pos at the index
  /// stored in \a be[pos]
  ///
  /// @param be: a LinkEntry
  /// @param pos: node number from 0..NNodes-1
  int GetInt(const LinkEntry& be, size_t pos) const
  { return GetNode(pos).GetInt(be[pos]); }
  /// SetInt
  void SetInt(const LinkEntry& be, size_t pos, int i)
  { GetNode(pos).SetInt(be[pos], i); }

  /// Access double value at the node \a pos at the index
  /// stored in \a be[pos]
  ///
  /// @param be: a LinkEntry
  /// @param pos: node number from 0..NNodes-1
  double GetDbl(const LinkEntry& be, size_t pos) const
  { return GetNode(pos).GetDbl(be[pos]); }
  /// SetDbl
  void SetDbl(const LinkEntry& be, size_t pos, double i)
  { GetNode(pos).SetDbl(be[pos], i); }

  /// Access string value at the node \a pos at the index
  /// stored in \a be[pos]
  ///
  /// @param be: a LinkEntry
  /// @param pos: node number from 0..NNodes-1
  const VCString& GetStr(const LinkEntry& be, size_t pos) const
  { return GetNode(pos).GetStr(be[pos]); }
  /// SetStr
  void SetStr(const LinkEntry& be, size_t pos, VCString v)
  { GetNode(pos).SetStr(be[pos], std::move(v)); }

private:
  NodeList ndl_;
};


/// A class providing autolinking in the destructor (RAII).
/// Usage: create an object, and all conversions in its scope
/// will be autolinked to the source node range (unless
/// autolinking is switched off during the execution).
template <class ModelConverter>
class AutoLinkScope {
public:
  /// Constructor.
  /// @param cvt: the ModelConverter doing the autolinking
  /// @param src: the source node for the links
  AutoLinkScope(ModelConverter& cvt, NodeRange src) : cvt_(cvt) {
    assert(src.IsSingleIndex());
    cvt_.SetAutoLinkSource(src);
    assert(cvt_.GetAutoLinkTargets().empty());
  }

  /// Destructor
  ~AutoLinkScope() {
    const auto& targets = cvt_.GetAutoLinkTargets();
    if (auto sz = targets.size()) {
      assert(cvt_.DoingAutoLinking());
      if (1==sz && targets.front().IsSingleIndex()) {
        cvt_.GetCopyLink().AddEntry(   // use CopyLink for 1:1
              { cvt_.GetAutoLinkSource(), targets.front() } );
      } else {
        for (const auto& t: targets) {
          assert(t.IsValid());
          cvt_.GetOne2ManyLink().AddEntry(   // use One2ManyLink
                { cvt_.GetAutoLinkSource(), t } );
        }
        cvt_.GetOne2ManyLink().MarkLastEntryAsLastInGroup();
      }
    }
    cvt_.TurnOffAutoLinking();
  }

private:
  ModelConverter& cvt_;
};

} // namespace pre

} // namespace mp

#endif // VALUE_PRESOLVE_LINK_H
