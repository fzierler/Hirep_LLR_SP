#ifndef __DEBUG_H_
#define __DEBUG_H_
#include <stddef.h>

// To be able to run GDB on MPI jobs, as desribed in
// https://www.open-mpi.org/faq/?category=debugging#serial-debuggers
int gethostname(char *name, size_t size);
void gdb_hook();

#endif // __DEBUG_H_
