#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "check_syscalls.h"
//#include <sys/file.h>
#include <sys/stat.h>
#include <assert.h>
#include "cached_io.h"

char *check_strdup(char *string) {
  char *newstring = check_realloc(NULL, (strlen(string)+1)*sizeof(char),
				"Allocating string copy space.\n");
  strcpy(newstring, string);
  return newstring;
}

struct cached_io *cfopen(char *filename, int64_t cache_size) {
  struct cached_io *cio = check_realloc(NULL, sizeof(struct cached_io),
					"Allocating new cached IO.");
  memset(cio, 0, sizeof(struct cached_io));
  cio->cachesize = cache_size;
  cio->filename = check_strdup(filename);
  return cio;
}

void _flush_wcache(struct cached_io *cio, char *extra, int64_t extra_size) {
  FILE *io;
  struct stat out_stats;
  if (cio->wcache_pos || extra_size) {
    if (cio->append) io = check_fopen(cio->filename, "a");
    else {
      io = fopen(cio->filename, "r+");
      if (!io) io = check_fopen(cio->filename, "w");
      fseeko(io, cio->wpos, SEEK_SET);
    }
    //if (cio->do_locks) flock(fileno(io), LOCK_EX);
    
    if (cio->append) fseek(io, 0, SEEK_END);
    if (cio->header) {
      fstat(fileno(io), &out_stats);
      if (!out_stats.st_size)
	fwrite(cio->header, 1, strlen(cio->header), io);
    }
    if (cio->wcache_pos) fwrite(cio->wcache, 1, cio->wcache_pos, io);
    if (extra) fwrite(extra, 1, extra_size, io);
    fflush(io);
    cio->wpos += cio->wcache_pos + extra_size;
    //if (cio->do_locks) flock(fileno(io), LOCK_UN);
    fclose(io);
    cio->wcache_pos = 0;
  }
}

void _fill_rcache(struct cached_io *cio) {
  FILE *io;
  io = check_fopen(cio->filename, "r");
  fseek(io, cio->rpos, SEEK_SET);
  cio->rcache_size = fread(cio->rcache, 1, cio->cachesize, io);
  cio->rcache[cio->rcache_size] = 0;
  cio->rpos = ftello(io);
  cio->rcache_pos = 0;
  fclose(io);
}

void cfclose(struct cached_io *cio) {
  _flush_wcache(cio, NULL, 0);
  if (cio->rcache) free(cio->rcache);
  if (cio->wcache) free(cio->wcache);
  if (cio->header) free(cio->header);
  free(cio->filename);
  free(cio);
}

void cfputs(struct cached_io *cio, char *line) {
  int64_t linelen = strlen(line);
  if (!cio->wcache) cio->wcache = check_realloc(NULL, cio->cachesize,
						"Allocating write cache.");
  if (linelen + cio->wcache_pos > cio->cachesize) {
    _flush_wcache(cio, line, linelen);
    return;
  }

  memcpy(cio->wcache + cio->wcache_pos, line, linelen);
  cio->wcache_pos += linelen;
}

int cfgets(struct cached_io *cio, char *buffer, int64_t maxlen) {
  int64_t bytes_copied = 0;
  int64_t to_copy = 0;
  char *end;
  if (!cio->rcache) {
    cio->rcache = check_realloc(NULL, cio->cachesize+1,
				"Allocating read cache.");
    _fill_rcache(cio);
  }

  if (cio->rcache_pos == cio->rcache_size) {
    if (cio->rcache_size < cio->cachesize) return 0;
    _fill_rcache(cio);
    if (cio->rcache_size == 0) return 0;
  }

  while (1) {
    end = strchr(cio->rcache + cio->rcache_pos, '\n');
    if (end) {
      to_copy = maxlen - bytes_copied - 1;
      if (to_copy > ((end - (cio->rcache + cio->rcache_pos)) + 1))
	to_copy = ((end - (cio->rcache + cio->rcache_pos)) + 1);
      memcpy(buffer + bytes_copied, cio->rcache + cio->rcache_pos, to_copy);
      buffer[bytes_copied + to_copy] = 0;
      cio->rcache_pos = (end - cio->rcache) + 1;
      return 1;
    }

    if (cio->rcache_size == 0) return 1;

    to_copy = maxlen - bytes_copied - 1;
    if (to_copy > (cio->rcache_size - cio->rcache_pos))
      to_copy = cio->rcache_size - cio->rcache_pos;
    memcpy(buffer+bytes_copied, cio->rcache + cio->rcache_pos, to_copy);
    bytes_copied += to_copy;
    buffer[bytes_copied] = 0;
    _fill_rcache(cio);
  }
}

