#pragma once
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <cassert>
#include <ctime>

// Replacement for FastPhylo log_utils.hpp: errors throw exceptions instead of exit(1).

#ifndef USE_PRINT
#define USE_PRINT 0
#endif

#ifndef PRINT
#define PRINT(EXP) \
    if(USE_PRINT){ std::cout << #EXP << " = " << (EXP) << std::endl; }
#endif

#ifndef PRINT_V
#define PRINT_V(EXP) \
    if(USE_PRINT){ std::cout << __FILE__ << ":" << __LINE__ << "  (" << #EXP << ") = " << (EXP) << std::endl; }
#endif

#ifndef PRINT_TIME
#define PRINT_TIME(EXP) \
    do { clock_t _t = clock(); EXP; \
    if(USE_PRINT) std::cout << __FILE__ << ":" << __LINE__ << "  (" << #EXP << ")  took  " \
        << double(clock()-_t)/(CLOCKS_PER_SEC/1000) << " ms" << std::endl; } while(0)
#endif

#ifndef PRINT_EXP
#define PRINT_EXP(EXP) \
    if(USE_PRINT){ std::cout << "executing: " << #EXP << std::endl; EXP; }
#endif

#ifndef LINE
#define LINE() \
    if(USE_PRINT){ std::cout << __FILE__ << ":" << __LINE__ << std::endl; }
#endif

#ifndef SEPARATOR
#define SEPARATOR() \
    if(USE_PRINT){ std::cout << __FILE__ << ":" << __LINE__ << "------------------------------" << std::endl; }
#endif

#ifndef ASSERT
#ifndef NDEBUG
#define ASSERT(EXP, P1, P2) assert(EXP)
#else
#define ASSERT(EXP, P1, P2)
#endif
#endif

#ifndef ASSERT_EQ
#ifndef NDEBUG
#define ASSERT_EQ(EXP1, EXP2) assert((EXP1)==(EXP2))
#else
#define ASSERT_EQ(EXP1, EXP2)
#endif
#endif

#ifndef ASSERT_OP
#ifndef NDEBUG
#define ASSERT_OP(EXP1, OP, EXP2) assert((EXP1)OP(EXP2))
#else
#define ASSERT_OP(EXP1, OP, EXP2)
#endif
#endif

#ifndef MEM_CHECK
#define MEM_CHECK(PTR) \
    do { if((PTR)==nullptr){ throw std::bad_alloc(); } } while(0)
#endif

#ifndef PROG_ERROR
#define PROG_ERROR(EXP) \
    do { std::ostringstream _oss; \
         _oss << "Internal error [" << __FILE__ << ":" << __LINE__ << "]: " << EXP; \
         throw std::runtime_error(_oss.str()); } while(0)
#endif

#ifndef USER_ERROR
#define USER_ERROR(EXP) \
    do { std::ostringstream _oss; \
         _oss << EXP; \
         throw std::invalid_argument(_oss.str()); } while(0)
#endif

#ifndef USER_WARNING
#define USER_WARNING(EXP) \
    std::cerr << "Warning: " << EXP << std::endl;
#endif

#ifndef USE_DEBUG
#define USE_DEBUG 0
#endif
#ifndef DEBUG
#define DEBUG(EXP) \
    if(USE_DEBUG){ std::cout << __FILE__ << ":" << __LINE__ << "  (" << #EXP << ") = " << EXP << std::endl; }
#endif
