#ifndef _CALC_DESC_H_
#define _CALC_DESC_H_

#include <stdint.h>

struct snap {
    int64_t num_h, num_p_total;
    int64_t *pid, *hid, *num_p_halo;
};

void calculate_descendants(int64_t *, struct snap *, struct snap *);

#endif /* _CALC_DESC_H */
