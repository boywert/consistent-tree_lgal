#ifndef CHECK_SYSCALLS_H
#define CHECK_SYSCALLS_H

#include <stdio.h>
FILE *check_fopen(char *filename, char *mode);
void *check_realloc(void *ptr, size_t size, char *reason);

#endif /* CHECK_SYSCALLS_H */
