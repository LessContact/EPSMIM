#ifndef AFFINITY_H
#define AFFINITY_H

#include <cstdlib> // For std::exit

#ifdef _WIN32
#include <windows.h>
#else
#define _GNU_SOURCE
#include <sched.h>
#include <sys/types.h>
#endif

// Define the namespace
namespace affinity {

    // Function declaration for setting CPU affinity
    void setAffinity(int core);

}  // namespace Affinity

#endif // AFFINITY_H