#include "mp/flat/redef/MIP/converter_mip.h"
#include "mp/flat/model_api_connect.h"

#include "mp2nlmodelapi.h"


namespace mp {

/// Defining the function in ...-modelapi-connect.cc
/// for recompilation speed
std::unique_ptr<BasicModelManager>
CreateMP2NLModelMgr(MP2NLCommon& cc, Env& e,
										 pre::BasicValuePresolver*& pPre) {
	return CreateModelMgrWithFlatConverter<
      MP2NLModelAPI, MIPFlatConverter >(cc, e, pPre);
}

}  // namespace mp
