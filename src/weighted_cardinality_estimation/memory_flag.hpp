#pragma once
#include <cstdint>

namespace MemoryFlag {
    constexpr uint64_t NOTHING                    = 0;
    constexpr uint64_t REGISTERS                  = 1 << 0;
    constexpr uint64_t ALL_WRITE_NO_REGISTERS     = 1 << 1;
    constexpr uint64_t FISHER_YATES_PERM_INIT     = 1 << 2;
    constexpr uint64_t FISHER_YATES_NON_PERM_INIT = 1 << 3;
    constexpr uint64_t SEEDS                      = 1 << 4;
    constexpr uint64_t TOTAL                      = 1 << 5;

    constexpr uint64_t ALL_WRITE    = REGISTERS | ALL_WRITE_NO_REGISTERS;
    constexpr uint64_t FISHER_YATES = FISHER_YATES_PERM_INIT | FISHER_YATES_NON_PERM_INIT;

    // Resolve flags: if TOTAL bit is set, other bits become EXCLUSIONS.
    inline uint64_t resolve(uint64_t flags) {
        constexpr uint64_t ALL_COMPONENTS = REGISTERS | ALL_WRITE_NO_REGISTERS
                                          | FISHER_YATES_PERM_INIT | FISHER_YATES_NON_PERM_INIT
                                          | SEEDS;
        if (flags & TOTAL) {
            return ALL_COMPONENTS & ~(flags & ~TOTAL);
        }
        return flags;
    }
}
