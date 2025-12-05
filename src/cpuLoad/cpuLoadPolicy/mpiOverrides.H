/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 Tyler Voskuilen
     \\/     M anipulation  |
-------------------------------------------------------------------------------
21-05-2020 Synthetik Applied Technologies: |   Modified original
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

#ifndef mpiOverrides_H
#define mpiOverrides_H

#include <mpi.h>
#include <array>
#include <string_view>
#include <chrono>
#include <mutex>
#include "Ostream.H"
#include "Pstream.H"

// This file attempts to override original MPI calls to profile processor loads
// for load-balancing related tasks
//
// To generate a list of MPI calls that can potentially be overriden in a codebase.
// These are the MPI function calls appearing in the code base and having PMPI counterparts
// ```bash
// tree-sitter query --scope source.cpp <query_path>
//             $(find <codebase_path> \( -name "*.H" -o -name "*.C" \) ! -path "*/lnInclude/*") |
//             sed -n 's/.*\(`\(MPI_\w*\)\).*/\2/p' |
//             sort | uniq |
//             while read -r mpi_func; do nm -D <path_to_libmpi.so> | grep "P${mpi_func}$"; done
// ```
// The query file should look like:
// ```query
// (
//   (call_expression
//     function: (identifier) @func
//     arguments: (argument_list))
//   (#match? @func "^MPI_")
// )
// ```

namespace Foam
{

enum class TimeUnit {
    Nano,
    Micro,
    Milli
};

// Runtime-configurable time unit (default: milliseconds)
// Can be set via dictionary in cpuLoadPolicy
inline TimeUnit& profilerTimeUnit() {
    static TimeUnit unit = TimeUnit::Milli;
    return unit;
}

// Always measure in nanoseconds internally, convert on output
using Duration = std::chrono::nanoseconds;

// Conversion factor from nanoseconds to the selected time unit
inline long long timeUnitDivisor(TimeUnit unit) {
    switch (unit) {
        case TimeUnit::Milli: return 1000000LL;
        case TimeUnit::Micro: return 1000LL;
        case TimeUnit::Nano:  return 1LL;
        default:              return 1000LL;
    }
}

inline const char* timeUnitSuffix(TimeUnit unit) {
    switch (unit) {
        case TimeUnit::Milli: return "ms";
        case TimeUnit::Micro: return "us";
        case TimeUnit::Nano:  return "ns";
        default:              return "us";
    }
}

inline const char* timeUnitSuffix() {
    return timeUnitSuffix(profilerTimeUnit());
}

struct MPICommsStats {
    long long p2pSendTime = 0;  // total P2P time inside send ops
    long long p2pRecvTime = 0;  // total P2P time inside receive ops
    long long collectiveTime = 0; // total time inside collective comms ops
    long long otherTime = 0;
    int nP2PSends = 0; // number of calls to P2P sends
    int nP2PRecvs = 0; // number of calls to P2P receives
    int nCollectives = 0; // number of calls to collective routines
    int nOthers = 0;

    // when latest cycle of LB has started
    std::chrono::high_resolution_clock::time_point cycleStart = std::chrono::high_resolution_clock::now();
    int nTSInCycle = 0;

    static constexpr std::array<std::string_view, 28> overridenMPICalls{
        "MPI_Send",
        "MPI_Bsend",
        "MPI_Isend",
        "MPI_Ssend",
        "MPI_Issend",
        "MPI_Recv",
        "MPI_Irecv",
        "MPI_Allgather",
        "MPI_Allreduce",
        "MPI_Bcast",
        "MPI_Barrier",
        "MPI_Wait",
        "MPI_Waitall",
        "MPI_Waitany",
        "MPI_Waitsome",
        "MPI_Test",
        "MPI_Testall",
        "MPI_Testany",
        "MPI_Testsome",
        "MPI_Accumulate",
        "MPI_Get",
        "MPI_Put",
        "MPI_Win_lock",
        "MPI_Win_unlock",
        "MPI_Win_flush",
        "MPI_Win_flush_all",
        "MPI_Win_flush_local",
        "MPI_Win_flush_local_all",
    };

    void reset() {
        p2pSendTime = 0;
        p2pRecvTime = 0;
        collectiveTime = 0;
        otherTime = 0;
        nP2PSends = 0;
        nP2PRecvs = 0;
        nCollectives = 0;
        nOthers = 0;
        nTSInCycle = 0;
        cycleStart = std::chrono::high_resolution_clock::now();
    }

    MPICommsStats operator+(const MPICommsStats& other) const {
        MPICommsStats result;
        result.p2pSendTime     = this->p2pSendTime + other.p2pSendTime;
        result.p2pRecvTime     = this->p2pRecvTime + other.p2pRecvTime;
        result.collectiveTime  = this->collectiveTime + other.collectiveTime;
        result.otherTime       = this->otherTime + other.otherTime;
        result.nP2PSends       = this->nP2PSends + other.nP2PSends;
        result.nP2PRecvs       = this->nP2PRecvs + other.nP2PRecvs;
        result.nCollectives    = this->nCollectives + other.nCollectives;
        result.nOthers         = this->nOthers + other.nOthers;
        return result;
    }

    Ostream& operator<<(Ostream& os) const {
        const char* unit = timeUnitSuffix();
        long long divisor = timeUnitDivisor(profilerTimeUnit());
        os << "P2P Send Time: " << label(p2pSendTime/divisor) << " " << unit << ", "
           << "P2P Recv Time: " << label(p2pRecvTime/divisor) << " " << unit << ", "
           << "Collective Time: " << label(collectiveTime/divisor) << " " << unit << ", "
           << "Other Time: " << label(otherTime/divisor) << " " << unit << "," << nl
           << "nP2PSends: " << nP2PSends << ", "
           << "nP2PRecvs: " << nP2PRecvs << ", "
           << "nCollectives: " << nCollectives << ", "
           << "nOthers: " << nOthers;
        return os;
    }
};

inline Ostream& operator<<(Ostream& os, const MPICommsStats& stats) {
    return stats.operator<<(os);
}

// Use inline to ensure single definition across translation units (C++17)
inline MPICommsStats mpiCommsStats;
inline std::mutex mpiCommsStatsMutex;

enum class MPITypeLevel { P2PSend, P2PRecv, Collective, Other };

template <MPITypeLevel Level, typename Ret, typename... Args>
Ret mpiCallProfiler
(
    Ret (*pmpi_func)(Args...),
    Args... args
) {
    auto start = std::chrono::high_resolution_clock::now();
    Ret ret = pmpi_func(std::forward<Args>(args)...);
    auto end = std::chrono::high_resolution_clock::now();
    long long elapsed = std::chrono::duration_cast<Duration>(end - start).count();

    std::lock_guard<std::mutex> lock(mpiCommsStatsMutex);
    if constexpr (Level == MPITypeLevel::P2PSend) {
        mpiCommsStats.p2pSendTime += elapsed;
        mpiCommsStats.nP2PSends++;
    } else if constexpr (Level == MPITypeLevel::P2PRecv) {
        mpiCommsStats.p2pRecvTime += elapsed;
        mpiCommsStats.nP2PRecvs++;
    } else if constexpr (Level == MPITypeLevel::Collective) {
        mpiCommsStats.collectiveTime += elapsed;
        mpiCommsStats.nCollectives++;
    } else if constexpr (Level == MPITypeLevel::Other) {
        mpiCommsStats.otherTime += elapsed;
        mpiCommsStats.nOthers++;
    }

    return ret;
}

}

#define WRAP_MPI_FUNCTION(Level, RetType, FuncName, Sig, ...) \
extern "C" RetType MPI_##FuncName Sig { \
    return Foam::mpiCallProfiler<Level>(PMPI_##FuncName, __VA_ARGS__); \
}

#endif

// ************************************************************************* //
