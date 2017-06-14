#ifndef MACROS_H
#define MACROS_H
#define PRINT_MAT(__A) std::cout << #__A << "=\n"; \
                       std::cout <<  __A << std::endl; 
#define PRINT(__A) std::cout << #__A << "="; \
                   std::cout <<  __A << std::endl; 
#endif
