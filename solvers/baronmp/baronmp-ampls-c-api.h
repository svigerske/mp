#ifndef BARONMPAMPLSCAPI_H
#define BARONMPAMPLSCAPI_H
/*
 * C API for MP/Baronmp
 */

//#include "baronmp.h"

#include "mp/ampls-c-api.h"

/*
 * Below are Baronmp-specific AMPLS API functions.
 * They complement the 'public' AMPLS API defined in ampls-c-api.h.
 */

/// Initialize AMPLS baronmp.

/// @param slv_opt: a string of solver options
/// (normally provided in the <solver>_options string).
/// Can be NULL.
/// @return pointer to struct AMPLS_MP_Solver to be populated.
void*  AMPLSOpenBaronmp(const char* slv_opt, CCallbacks cb);

/// Shut down solver instance
void AMPLSCloseBaronmp(AMPLS_MP_Solver* slv);

/// Extract the Baronmp model handle
void* GetBaronmpmodel(AMPLS_MP_Solver* slv);


#endif // BARONMPAMPLSCAPI_H


