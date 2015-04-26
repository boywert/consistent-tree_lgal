#include <stdlib.h>
#include <stdio.h>
#include "calc_desc.h"
#include "check_syscalls.h"
#include "inthash.h"

int compare_int64(const void *a, const void *b) {
    int64_t c = *((int64_t *)a);
    int64_t d = *((int64_t *)b);
    if (c < d) return -1;
    if (c > d) return 1;
    return 0;
}

void calculate_descendants(int64_t * desc, struct snap *s1, struct snap *s2) {
    /* s2 is the later snap */

    int64_t i, j, k, p_count, max_num_p, cand, cand_count, flag;
    int64_t * new_hid_idx = NULL;
    struct inthash * hash = NULL;

    /* initialize descendants */
    if (!s1->num_p_total || !s1->num_h) return;
    for (i=0; i<s1->num_h; ++i) desc[i] = -1;

    /* generate hash table */
    if (!s2->num_p_total || !s2->num_h) return;
    hash = new_inthash();
    ih_prealloc(hash, s2->num_p_total);
    p_count = 0;
    for (i=0; i<s2->num_h; ++i) {
        for (j=p_count; j<p_count+s2->num_p_halo[i]; ++j){
            k = ih_getint64(hash, s2->pid[j]);
            if (k == IH_INVALID || s2->num_p_halo[k] > s2->num_p_halo[i])
                ih_setint64(hash, s2->pid[j], i);
        }
        p_count += s2->num_p_halo[i];
    }

    /* match particles in s1 to s2 */
    p_count = 0;
    max_num_p = 0;
    for (i=0; i<s1->num_h; ++i) {
        if (s1->num_p_halo[i] > max_num_p) {
            max_num_p = s1->num_p_halo[i];
            new_hid_idx = check_realloc(new_hid_idx, sizeof(int64_t)*max_num_p,
			        "calc_desc.c: Allocate new_hid_idx.");
        }

        k=0;
        for (j=p_count; j<p_count+s1->num_p_halo[i]; ++j){
            new_hid_idx[k] = ih_getint64(hash, s1->pid[j]);
            if (new_hid_idx[k] < IH_INVALID) ++k;
        }
        p_count += s1->num_p_halo[i];
        if (!k) continue;

        /* find the new halo with most particles */
        qsort(new_hid_idx, k, sizeof(int64_t), compare_int64);
        cand = new_hid_idx[0];
        cand_count = 0;
        flag = 0;
        for (j=1; j<k; ++j){
            if (new_hid_idx[j] != new_hid_idx[flag]){
                if (j-flag > cand_count){
                    cand = new_hid_idx[flag];
                    cand_count = j-flag;
                }
                flag = j;
            }
        }
        desc[i] = s2->hid[cand];
    }

    free_inthash(hash);
    new_hid_idx = check_realloc(new_hid_idx, 0, 
            "calc_desc.c: free new_hid_idx.");
}

