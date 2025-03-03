#include "affinity.h"
#include <iostream>
#include <cerrno>
#include <cstring>

namespace affinity {

    void setAffinity(int core) {
#ifdef _WIN32
        // Windows-specific implementation
        DWORD_PTR affinityMask = (1 << core);
        HANDLE hProcess = GetCurrentProcess();
        if (!SetProcessAffinityMask(hProcess, affinityMask)) {
            // Retrieve the error code
            DWORD error = GetLastError();
            LPVOID errorMessage;
            FormatMessage(
                FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
                NULL,
                error,
                MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
                (LPTSTR)&errorMessage,
                0,
                NULL
            );
            // Print the error message
            std::cerr << "SetProcessAffinityMask failed with error " << error << ": " << (char*)errorMessage << std::endl;
            // Free the buffer allocated by FormatMessage
            LocalFree(errorMessage);
            std::exit(EXIT_FAILURE);
        }
#else
        // POSIX-specific implementation
        cpu_set_t mask;
        CPU_ZERO(&mask);
        CPU_SET(core, &mask);
        if (sched_setaffinity(0, sizeof(mask), &mask) == -1) {
            std::cerr << "sched_setaffinity failed: " << strerror(errno) << std::endl;
            std::exit(EXIT_FAILURE);
        }
#endif
    }

}  // namespace Affinity