/* AHF Conversion from Yao-Yuan Mao */
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <stdint.h>
#include <string.h>
#include "src/calc_desc.h"

int64_t * realloc_int64(int64_t* ptr, size_t count, char * varname){
    int64_t * res = (int64_t *) realloc(ptr, sizeof(int64_t)*count);
    if(!count) return NULL;
    if(!res){
        fprintf(stderr, "[Error] Failed to allocate %s!\n", varname);
        exit(4);
    }
    return res;
}

FILE * fopen_check(char *filename, char *mode) {
    FILE * res = fopen(filename, mode);
    if (!res) {
        fprintf(stderr, "[Error] Failed to open file %s\n", filename);
        exit(2);
    }
    return res;
}

int main(int argc, char **argv) {

    FILE * f_list, * f_snap, * f_out;
    char buffer[1024];
    struct snap s1, s2;
    int64_t i, j, p_count;
    int snap_count = -1;
    int64_t * desc = NULL;

    if (argc < 2) {
        fprintf(stderr, "Usage: %s <snapshot list>\n", argv[0]); 
        exit(1);
    }

    f_list = fopen_check(argv[1], "r");

    memset(&s1, 0, sizeof(struct snap));
    memset(&s2, 0, sizeof(struct snap));
    
    while (fgets(buffer, 1024, f_list)) {

        buffer[strlen(buffer)-1] = '\0';
        f_snap = fopen_check(buffer, "r");
        fgets(buffer, 1024, f_snap);
        sscanf(buffer, "%"PRId64"", &(s2.num_h));
        s2.hid = realloc_int64(NULL, s2.num_h, "s2.hid");
        s2.num_p_halo = realloc_int64(NULL, s2.num_h, "s2.num_p_halo");
        s2.pid = NULL;

        p_count = 0;
        for (i=0; i<s2.num_h; ++i){
            fgets(buffer, 1024, f_snap);
            sscanf(buffer, "%"PRId64" %"PRId64"", &(s2.num_p_halo[i]), 
                    &(s2.hid[i]));
            p_count += s2.num_p_halo[i];
            s2.pid = realloc_int64(s2.pid, p_count, "s2.pid");
            for (j=p_count-s2.num_p_halo[i]; j<p_count; ++j){
                fgets(buffer, 1024, f_snap);
                sscanf(buffer, "%"PRId64"", &(s2.pid[j]));
            }
        }
        s2.num_p_total = p_count;
        fclose(f_snap);

        if (snap_count > -1){
            desc = realloc_int64(desc, s1.num_h, "desc");
            calculate_descendants(desc, &s1, &s2);
            free(s1.pid);
            free(s1.num_p_halo);

            sprintf(buffer, "desc_%d.bin", snap_count);
            f_out = fopen(buffer, "wb");
            fwrite(s1.hid, sizeof(int64_t), s1.num_h, f_out);
            fwrite(desc,   sizeof(int64_t), s1.num_h, f_out);
            fflush(f_out);
            fclose(f_out);

            free(s1.hid);
        }
        s1 = s2;
        ++snap_count;
    }

    free(desc);
    return 0;
}
