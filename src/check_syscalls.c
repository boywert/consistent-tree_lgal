#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <inttypes.h>

FILE *check_fopen(char *filename, char *mode) {
  FILE *res = fopen(filename, mode);
  if (res == NULL) {
    if (mode[0] == 'w')
      fprintf(stderr, "Failed to open file %s for writing!\n", filename);
    else if (mode[0] == 'a')
      fprintf(stderr, "Failed to open file %s for appending!\n", filename);
    else
      fprintf(stderr, "Failed to open file %s for reading!\n", filename);
    exit(1);
  }
  return res;
}

void *check_realloc(void *ptr, size_t size, char *reason) {
  void *res = realloc(ptr, size);
  if ((res == NULL) && (size > 0)) {
    fprintf(stderr, "[Error] Failed to allocate %"PRId64" bytes of memory (%s)!\n", (int64_t)size, reason);
    assert(0);
  }
  return res;
}
