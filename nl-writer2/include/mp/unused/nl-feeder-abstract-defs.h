#ifndef NLFEEDERABSTRACTDEFS_H
#define NLFEEDERABSTRACTDEFS_H


namespace mp {

/// SparseVectorWriter<> interface
template <class T>
class BasicSparseVectorWriter {
public:
  /// Virtual destruct
  virtual ~BasicSparseVectorWriter() { }
  /// Number of elements left to write
  virtual int NLeft() const = 0;
  /// Write entry
  virtual void Write(int index, T value) = 0;
};

/// Typedef sparse int vec writer
using BasicSparseIntVecWriter = BasicSparseVectorWriter<int>;

/// Typedef sparse double vec writer
using BasicSparseDblVecWriter = BasicSparseVectorWriter<double>;


/// SparseVectorWriterFactory interface.
/// An implementation can be used once to produce a SparseVectorWriter<>.
template <class T>
class BasicSparseVectorWriterFactory {
public:
  /// Virtual destruct
  virtual ~BasicSparseVectorWriterFactory() { }
  /// Make the SparseVectorWriter.
  /// @note std::size_t without include:
  ///   https://stackoverflow.com/questions/36594569/which-header-should-i-include-for-size-t.
  virtual BasicSparseVectorWriter<T>&
  MakeVectorWriter(decltype(sizeof(int)) nnz) = 0;
};

/// Typedef sparse int vec writer factory
using BasicSparseIntVecWrtFactory = BasicSparseVectorWriterFactory<int>;

/// Typedef sparse double vec writer factory
using BasicSparseDblVecWrtFactory = BasicSparseVectorWriterFactory<double>;


/// NLExprWriter interface.
/// @note Store ExprArgWriter's as local objects.
template <class Expr>
class BasicNLExprWriter {
public:
  /// Virtual destruct
  virtual ~BasicNLExprWriter() { }

  /// Write the next arg as Feeder's native expression.
  /// This recursively calls Feeder::FeedExpr().
  virtual void EPut(Expr e) = 0;

  /** Write the next arg as 'variable reference'.
     *  0 <= index < num_vars is a solver variable;
     *  index >= num_vars is a defined variable. */
  virtual void VPut(int v, const char* descr="") = 0;

  /** Write numeric constant expression. */
  virtual void NPut(double x) = 0;

  /** Write string constant expression. */
  virtual void StrPut(const char* ) = 0;


  /// typedef ExprArgWriter
  /// @todo return temporary
  using ExprArgWriter = int; // BasicNLExprWriter;

  /** Write the next arg as function call expression. */
  virtual ExprArgWriter FuncPut(
      int index, int nArgs, const char* descr="") = 0;

  /// Write the next arg as AMPL opcode for a unary op.
  /// @return 1-arg writer.
  virtual ExprArgWriter OPut1(int opcode, const char* descr="") = 0;
  /// Write AMPL opcode for a binary op.
  virtual ExprArgWriter OPut2(int opcode, const char* descr="") = 0;
  /// Write AMPL opcode for a 3-arg op.
  virtual ExprArgWriter OPut3(int opcode, const char* descr="") = 0;

  /// Write AMPL opcode for an iterated op (min, exists, sum, etc).
  /// For a piecewise-linear expression, \a nArgs should be
  /// 2*(N slopes) and the arguments are:
  /// break points, slopes, argument variable.
  virtual ExprArgWriter OPutN(
      int opcode, int nArgs, const char* descr="") = 0;

  /// Shortcut: OPut1( struct Opcode )
  template <class Opcode> ExprArgWriter OPut1(Opcode oc)
  { return OPut1(oc.code, oc.name); }
  /// Shortcut: OPut2( struct Opcode )
  template <class Opcode> ExprArgWriter OPut2(Opcode oc)
  { return OPut2(oc.code, oc.name); }
  /// Shortcut: OPut3( struct Opcode )
  template <class Opcode> ExprArgWriter OPut3(Opcode oc)
  { return OPut3(oc.code, oc.name); }
  /// Shortcut: OPutN( struct Opcode, int nArgs )
  template <class Opcode> ExprArgWriter OPutN(
      Opcode oc, int nArgs)
  { return OPutN(oc.code, nArgs, oc.name); }
};


/** Interface: writer of a defined variable. */
template <class Expr>
class BasicNLDefVarWriter {
public:
  /// Virtual destruct
  virtual ~BasicNLDefVarWriter() { }

  /// Write entries c*var[v] to the linear part
  /// of the defining expression.
  /// All nnz entries should be written before
  /// the nonlinear expression.
  /// This method can be called only once.
  virtual BasicSparseDblVecWriter& GetLinExprWriter() = 0;

  /// Retrieve the nonlinear expression writer.
  /// Should be used after the linear expression.
  virtual BasicNLExprWriter<Expr>& GetExprWriter() = 0;
};


/** Interface: BasicNLDefVarWriterFactory.
 *  Write a sequence of defined variables,
 *  for example, all such used in a certain constraint. */
template <class Expr>
class BasicNLDefVarWriterFactory {
public:
  /// Virtual destruct.
  virtual ~BasicNLDefVarWriterFactory() { }

  /// Start writing a defined variable.
  ///
  /// @param index: defined variable index, used to
  /// reference it in subsequent expression graphs.
  /// Thus, the index should be >= num_vars,
  ///   < num_vars+num_common_exprs.
  /// Providing the index explicitly, because classical NL
  /// likes special order of defined variables.
  ///
  /// @param nnz: number of nonzeros in the linear part.
  ///
  /// @param descr: meta-information - what is this variable,
  /// for example, "nl(t[2])".
  /// Providing it here because it's not
  /// included in ColNames().
  ///
  /// @return A callback object writing
  /// a single defined variable.
  virtual BasicNLDefVarWriter<Expr>& StartDefVar(
      int index, int nnz, const char* descr="") = 0;
};


/** Interface: write \a num_vars variable bounds
   *  (all except defined variables). */
class BasicNLVarBndWriter {
public:
  /// virtual destruct
  virtual ~BasicNLVarBndWriter() { }
  /// Write range for the next variable.
  virtual void WriteLbUb(double lb, double ub) = 0;
};


/** Interface: write \a num_algebraic_cons constraint bounds. */
class BasicNLConBndWriter {
public:
  /// Virtual destruct
  virtual ~BasicNLConBndWriter() { }
  /// Write range/complementarity for the next constraint.
  virtual void WriteAlgConRange(
      double L, double U,
      int k, int cvar) = 0;
};


/// Interface: int suffix writer
using BasicNLSuffixIntWriter = BasicSparseIntVecWriter;

/// Interface: double-valued suffix writer
using BasicNLSuffixDblWriter = BasicSparseDblVecWriter;

/** Interface: suffixes. */
class BasicNLSuffixWriterFactory {
public:
  /// Virtual destruct
  virtual ~BasicNLSuffixWriterFactory() { }
  /// Start writing an int-valued suffix.
  virtual BasicNLSuffixIntWriter& StartIntSuffix(
      const char* name, int kind, int nnz) = 0;
  /// Start writing a dbl-valued suffix.
  virtual BasicNLSuffixDblWriter& StartDblSuffix(
      const char* name, int kind, int nnz) = 0;
};


/** Interface: PL-SOS constraints. */
class BasicPLSOSWriter {
public:
  /// Virtual destruct.
  virtual ~BasicPLSOSWriter() { }
  /// .sos for variables
  virtual BasicNLSuffixIntWriter& StartSOSVars(int nnz) = 0;
  /// .sos for constraints
  virtual BasicNLSuffixIntWriter& StartSOSCons(int nnz) = 0;
  /// .sosref for variables
  virtual BasicNLSuffixDblWriter& StartSOSREFVars(int nnz) = 0;
};


/** Interface: write num_vars+num_rand_vars-1
   *  column sizes */
class BasicNLColSizeWriter {
public:
  /// Virtual destruct
  virtual ~BasicNLColSizeWriter() { }

  /// Write next col's size
  virtual void Write(int s) = 0;
};


/// Interface: string file writer
class BasicNLNameFileWriter {
public:
  /// Virtual destruct
  virtual ~BasicNLNameFileWriter() { }

  /// operator bool
  virtual operator bool() const = 0;

  /// Write 1 string
  virtual void Write(const char* name) = 0;

  /// Write 2 strings (name, comment)
  virtual void Write(const char* name, const char* cmt) = 0;

  /// Write 1 string, 1 double
  virtual void Write(const char* name, double num) = 0;
};

}  // namespace mp

#endif // NLFEEDERABSTRACTDEFS_H
