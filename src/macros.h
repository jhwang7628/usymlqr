#ifndef MACROS_H
#define MACROS_H
#define PRINT_MAT(__stream, __A) __stream << #__A << "=\n"; \
                                 __stream <<  __A << std::endl; 
#define PRINT(__stream, __A) __stream << #__A << "="; \
                             __stream <<  __A << std::endl; 
#define SMALL_NUM 1E-10
#define DEFAULT_TOL 1E-12
#define MAX_VERBOSE_LEVEL 2

#endif
