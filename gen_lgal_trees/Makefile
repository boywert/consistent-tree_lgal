CFLAGS=-m64 -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -D_BSD_SOURCE -D_POSIX_SOURCE -D_POSIX_C_SOURCE=200809L -D_SVID_SOURCE -D_DARWIN_C_SOURCE -Wall  
SUPPORT= masses.c distance.c universe_time.c stringparse.c halo_io.c check_syscalls.c
EXTRA_SUPPORT=$(SUPPORT) read_config.c grav_config.c
OFLAGS = -lm -O3 -std=c99
DEBUGFLAGS = -lm -g -O -std=c99
PROFFLAGS = -lm -g -pg -O2 -std=c99
#CC = gcc-mp-4.8

all:
	@make reg EXTRA_FLAGS="$(OFLAGS)"

debug:
	@make reg EXTRA_FLAGS="$(DEBUGFLAGS)"

reg:
	$(CC) $(CFLAGS) gen_lgal_trees.c read_tree.c $(EXTRA_SUPPORT) -o ../gen_lgaltrees  $(EXTRA_FLAGS)

clean:
	rm -f *~
	rm -f ../gen_lgaltrees
