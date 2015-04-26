#ifndef _CACHED_IO_H_
#define _CACHED_IO_H_
#include <inttypes.h>

struct linecache {
  char *line;
  struct linecache *next;
};

struct cached_io {
  char *filename;
  int64_t do_locks;
  int64_t cachesize, rpos, wpos, append;
  int64_t rcache_size;
  int64_t wcache_pos, rcache_pos;
  char *header;
  char *wcache, *rcache;
};

struct cached_io *cfopen(char *filename, int64_t cache_size);
void cfclose(struct cached_io *cio);
void cfputs(struct cached_io *cio, char *line);
int cfgets(struct cached_io *cio, char *buffer, int64_t maxlen);

#endif /* _CACHED_IO_H_ */
