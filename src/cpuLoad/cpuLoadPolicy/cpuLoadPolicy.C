/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 Tyler Voskuilen
     \\/     M anipulation  |
-------------------------------------------------------------------------------
21-05-2020 Synthetik Applied Technologies: |    Modified original
                            dynamicRefineBalanceBlastFvMesh class
                            to be more appilcable to compressible flows.
                            Improved compatibility with snappyHexMesh.
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "cpuLoadPolicy.H"
#include "addToRunTimeSelectionTable.H"
#include "mpi.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cpuLoadPolicy, 0);
    addToRunTimeSelectionTable(loadPolicy, cpuLoadPolicy, dictionary);

    __attribute__((constructor))
    void onCPULoadLib() {
        WarningInFunction
            << "'libamrCPULoad.so' library loaded. Some MPI calls are overrided just by loading this library!"
            << nl << tab << "This is the case even if you don't use cpuLoad for load-balancing..."
            << nl << tab << "Of course this is possible only when this library is loaded before real MPI libs."
            << nl << endl;
    }
}

// MPI call overrides
WRAP_MPI_FUNCTION(
    Foam::MPITypeLevel::P2PSend,
    int,
    Send,
    (const void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm),
    buf, count, datatype, dest, tag, comm
)
WRAP_MPI_FUNCTION(
    Foam::MPITypeLevel::P2PSend,
    int,
    Isend,
    (const void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request* request),
    buf, count, datatype, dest, tag, comm, request
)
WRAP_MPI_FUNCTION(
    Foam::MPITypeLevel::P2PSend,
    int,
    Bsend,
    (const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm),
    buf, count, datatype, dest, tag, comm
)
WRAP_MPI_FUNCTION(
    Foam::MPITypeLevel::P2PSend,
    int,
    Ssend,
    (const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm),
    buf, count, datatype, dest, tag, comm
)
WRAP_MPI_FUNCTION(
    Foam::MPITypeLevel::P2PSend,
    int,
    Issend,
    (const void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request* request),
    buf, count, datatype, dest, tag, comm, request
)
WRAP_MPI_FUNCTION(
    Foam::MPITypeLevel::P2PRecv,
    int,
    Recv,
    (void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status* status),
    buf, count, datatype, source, tag, comm, status
)
WRAP_MPI_FUNCTION(
    Foam::MPITypeLevel::P2PRecv,
    int,
    Irecv,
    (void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request* request),
    buf, count, datatype, source, tag, comm, request
)
WRAP_MPI_FUNCTION(
    Foam::MPITypeLevel::Collective,
    int,
    Allgather,
    (const void* sendbuf, int sendcount, MPI_Datatype sendtype,
     void* recvbuf, int recvcount, MPI_Datatype recvtype,
     MPI_Comm comm),
    sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm
)
WRAP_MPI_FUNCTION(
    Foam::MPITypeLevel::Collective,
    int,
    Allreduce,
    (const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype,
     MPI_Op op, MPI_Comm comm),
    sendbuf, recvbuf, count, datatype, op, comm
)
WRAP_MPI_FUNCTION(
    Foam::MPITypeLevel::Collective,
    int,
    Bcast,
    (void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm),
    buffer, count, datatype, root, comm
)
WRAP_MPI_FUNCTION(
    Foam::MPITypeLevel::Collective,
    int,
    Barrier,
    (MPI_Comm comm),
    comm
)
WRAP_MPI_FUNCTION(
    Foam::MPITypeLevel::Other,
    int,
    Wait,
    (MPI_Request *request, MPI_Status *status),
    request, status
)
WRAP_MPI_FUNCTION(
    Foam::MPITypeLevel::Other,
    int,
    Waitall,
    (int count, MPI_Request *array_of_requests, MPI_Status *array_of_statuses),
    count, array_of_requests, array_of_statuses
)
WRAP_MPI_FUNCTION(
    Foam::MPITypeLevel::Other,
    int,
    Waitany,
    (int count, MPI_Request *array_of_requests, int *index, MPI_Status *status),
    count, array_of_requests, index, status
)
WRAP_MPI_FUNCTION(
    Foam::MPITypeLevel::Other,
    int,
    Waitsome,
    (int incount, MPI_Request *array_of_requests, int *outcount, int *array_of_indices, MPI_Status *array_of_statuses),
    incount, array_of_requests, outcount, array_of_indices, array_of_statuses
)
WRAP_MPI_FUNCTION(
    Foam::MPITypeLevel::Other,
    int,
    Test,
    (MPI_Request *request, int *flag, MPI_Status *status),
    request, flag, status
)
WRAP_MPI_FUNCTION(
    Foam::MPITypeLevel::Other,
    int,
    Testall,
    (int count, MPI_Request *array_of_requests, int *flag, MPI_Status *array_of_statuses),
    count, array_of_requests, flag, array_of_statuses
)
WRAP_MPI_FUNCTION(
    Foam::MPITypeLevel::Other,
    int,
    Testany,
    (int count, MPI_Request *array_of_requests, int *index, int *flag, MPI_Status *status),
    count, array_of_requests, index, flag, status
)
WRAP_MPI_FUNCTION(
    Foam::MPITypeLevel::Other,
    int,
    Testsome,
    (int incount, MPI_Request *array_of_requests, int *outcount, int *array_of_indices, MPI_Status *array_of_statuses),
    incount, array_of_requests, outcount, array_of_indices, array_of_statuses
)
WRAP_MPI_FUNCTION(
    Foam::MPITypeLevel::Other,
    int,
    Accumulate,
    (const void* origin_addr, int origin_count, MPI_Datatype origin_datatype,
     int target_rank, MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype,
     MPI_Op op, MPI_Win win),
    origin_addr, origin_count, origin_datatype, target_rank, target_disp,
    target_count, target_datatype, op, win
)
WRAP_MPI_FUNCTION(
    Foam::MPITypeLevel::Other,
    int,
    Get,
    (void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
     int target_rank, MPI_Aint target_disp, int target_count,
     MPI_Datatype target_datatype, MPI_Win win),
    origin_addr, origin_count, origin_datatype, target_rank, target_disp,
    target_count, target_datatype, win
)
WRAP_MPI_FUNCTION(
    Foam::MPITypeLevel::Other,
    int,
    Put,
    (const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
     int target_rank, MPI_Aint target_disp, int target_count,
     MPI_Datatype target_datatype, MPI_Win win),
    origin_addr, origin_count, origin_datatype, target_rank, target_disp,
    target_count, target_datatype, win
)
WRAP_MPI_FUNCTION(
    Foam::MPITypeLevel::Other,
    int,
    Win_lock,
    (int lock_type, int rank, int assert, MPI_Win win),
    lock_type, rank, assert, win
)
WRAP_MPI_FUNCTION(
    Foam::MPITypeLevel::Other,
    int,
    Win_unlock,
    (int rank, MPI_Win win),
    rank, win
)
WRAP_MPI_FUNCTION(
    Foam::MPITypeLevel::Other,
    int,
    Win_flush,
    (int rank, MPI_Win win),
    rank, win
)
WRAP_MPI_FUNCTION(
    Foam::MPITypeLevel::Other,
    int,
    Win_flush_all,
    (MPI_Win win),
    win
)
WRAP_MPI_FUNCTION(
    Foam::MPITypeLevel::Other,
    int,
    Win_flush_local,
    (int rank, MPI_Win win),
    rank, win
)
WRAP_MPI_FUNCTION(
    Foam::MPITypeLevel::Other,
    int,
    Win_flush_local_all,
    (MPI_Win win),
    win
)



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace
{
    // Helper to parse time unit from string
    Foam::TimeUnit parseTimeUnit(const Foam::word& unitStr)
    {
        if (unitStr == "nano" || unitStr == "nanoseconds" || unitStr == "ns")
        {
            return Foam::TimeUnit::Nano;
        }
        else if (unitStr == "milli" || unitStr == "milliseconds" || unitStr == "ms")
        {
            return Foam::TimeUnit::Milli;
        }
        // Default to microseconds
        return Foam::TimeUnit::Micro;
    }
}

Foam::cpuLoadPolicy::cpuLoadPolicy
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    loadPolicy(mesh, dict),
    maxCycleLength_(dict.lookupOrDefault("maxLBCycleLength", 5*readLabel(dict.lookup("refineInterval")))),
    isActive_(true),
    timeUnit_(parseTimeUnit(dict.getOrDefault<word>("timeUnit", "micro")))
{
    // Set the global time unit for output formatting
    profilerTimeUnit() = timeUnit_;

    Info<< "    timeUnit: " << timeUnitSuffix(timeUnit_) << endl;

    mpiCommsStats.reset();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cpuLoadPolicy::~cpuLoadPolicy()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::cpuLoadPolicy::canBalance()
{
    if (isActive_ && mesh_.time().timeIndex() > 5) {
        // check that MPI calls where intercepted
        // it's highly unlikely that after 5 iterations, no measuremants were picked up!
        if (
            returnReduce(mpiCommsStats.nP2PSends, sumOp<int>()) == 0 &&
            returnReduce(mpiCommsStats.nP2PRecvs, sumOp<int>()) == 0 &&
            returnReduce(mpiCommsStats.nCollectives, sumOp<int>()) == 0 &&
            returnReduce(mpiCommsStats.nOthers, sumOp<int>()) == 0
        ) {
            FatalErrorInFunction
                << "Seems like MPI call interception is not working. Make sure to preload the intercepting libraries:"
                << nl << nl << "LD_PRELOAD=\"$FOAM_USER_LIBBIN/libamrLoadPolicies.so $FOAM_USER_LIBBIN/libamrCPULoad.so\" <your-solver-command>"
                << nl << nl<< "Library loading order is important, libamrCPULoad.so must load before MPI libs."
                << nl << "So, no point in continuing..."
                << abort(FatalError);
        }
    }
    // Stat a new measuring cycle when LB is viable, bound by maxCycleLength_ of
    // time steps
    bool newCycle = !returnReduce(stillSameCycle(), orOp<bool>());
    if (mpiCommsStats.nTSInCycle > maxCycleLength_) newCycle = true;
    auto now = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<Duration>(now-mpiCommsStats.cycleStart).count();
    if (elapsed == 0) return false;
    if (newCycle) myLoad_ = 0;
    myLoad_ += elapsed
        - mpiCommsStats.p2pSendTime - mpiCommsStats.p2pRecvTime
        - mpiCommsStats.collectiveTime - mpiCommsStats.otherTime;
    scalar idealLoad = returnReduce(myLoad_, sumOp<scalar>()) / scalar(Pstream::nProcs());
    scalar maxImbalance = returnReduce(mag(myLoad_ - idealLoad) / idealLoad, maxOp<scalar>());
    Pout<< "Maximum imbalance found = " << 100*maxImbalance << " %" << endl;
    Pout<< "Current processor load stats:" << nl << mpiCommsStats << endl;
    if (maxImbalance < allowedImbalance_)
    {
        return false;
    }
    // start new cylce only if rebalancing...
    mpiCommsStats.reset();
    myLoadHistory_.set(mesh_.time().timeIndex(), myLoad_);
    //Pout << "Last history datapoint: " << myLoadHistory_.toc() << tab << myLoadHistory_[myLoadHistory_.toc().last()] << endl;
    if (!myLoadHistory_.empty()) {
        Pout<< "Current processor load computed from Time index = "
        << myLoadHistory_.toc().last()
        << " to Time index = " << mesh_.time().timeIndex() << "s"
        << nl;
    }
    return true;
}

Foam::scalarField Foam::cpuLoadPolicy::cellWeights() {
    label idealLoad = returnReduce(myLoad_, sumOp<scalar>()) / scalar(Pstream::nProcs());
    scalar imbalance = mag(myLoad_ - idealLoad) / idealLoad;
    if (myLoad_ > idealLoad) imbalance = 1.0/imbalance;
    return scalarField(mesh_.nCells(), myLoad_ != 0 ? imbalance : 1.0);
}

bool Foam::cpuLoadPolicy::willBeBeneficial
(
    const labelList distribution
) {
    // WARN: hard to know beforehand which cells will be beneficial
    // to move; it's a very dynamic property anyway, so no attempt
    // to make any guesses
    return true;
}

// ************************************************************************* //
