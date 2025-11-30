#include "abortHandle.H"
#include "catch2/catch_test_macros.hpp"

jmp_buf amr_env;

void onSigabrt(int signum)
{
    signal(signum, SIG_DFL);
    longjmp(amr_env, 1);
}

void tryAndCatchAbortingCode(std::function<void(void)> func)
{
    FatalError.dontThrowExceptions();
    if (setjmp(amr_env) == 0)
    {
        signal(SIGUSR1, &onSigabrt);
        signal(SIGUSR2, &onSigabrt);
        signal(SIGABRT, &onSigabrt);
        signal(SIGTERM, &onSigabrt);
        signal(SIGQUIT, &onSigabrt);
        func();
        signal(SIGQUIT, SIG_DFL);
        signal(SIGUSR1, SIG_DFL);
        signal(SIGUSR2, SIG_DFL);
        signal(SIGABRT, SIG_DFL);
        signal(SIGTERM, SIG_DFL);
    }
    else
    {
        Pout<< "Either this code tried to abort or there was"
            " an attempt to terminate it (e.g. with a timeout) on " <<
            Pstream::myProcNo() << "..." << endl;
        bool abortedOrTerminated = true;
        REQUIRE(abortedOrTerminated == false);
    }
}
