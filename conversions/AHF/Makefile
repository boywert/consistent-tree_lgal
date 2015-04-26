CFLAGS=-m64 -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -D_BSD_SOURCE -D_POSIX_SOURCE -D_POSIX_C_SOURCE=200809L -D_SVID_SOURCE -D_DARWIN_C_SOURCE -Wall
OFLAGS = -lm -O3 -std=c99
CC = gcc

CTREES_SRC_PATH = ../../src
SUPPORT_DESC = src/calc_desc.c $(CTREES_SRC_PATH)/check_syscalls.c $(CTREES_SRC_PATH)/inthash.c

all: get_descendants

get_descendants: get_descendants.c $(SUPPORT_DESC)
	$(CC) $(CFLAGS) $^ -o $@ -I$(CTREES_SRC_PATH) $(OFLAGS) 

clean:
	rm -f get_descendants
