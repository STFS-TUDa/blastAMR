#ifndef ABORT_HANDLE_H
#define ABORT_HANDLE_H

#include <csetjmp>
#include <csignal>
#include <cstdlib>
#include <functional>
#include "fvCFD.H"


// This is to make sure OpenFOAM doesn't bail on us at any point
// Typically, some functions with abort even if you set throwExceptions; which is not cool
// also, code might run into deadlocks or infinite loops, timeout will kill try to kill the process
// This section protects against these scenatios
// ********************************************************************************************** //
using namespace Foam;

// Declared in abortHandle.C
extern jmp_buf amr_env;
void onSigabrt(int signum);
void tryAndCatchAbortingCode(std::function<void(void)> func);

// ********************************************************************************************** //

#endif // ABORT_HANDLE_H
