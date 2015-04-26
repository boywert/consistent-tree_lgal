#ifndef _GRAVITATIONAL_CONSISTENCY_VARS_H_
#define _GRAVITATIONAL_CONSISTENCY_VARS_H_

#include <stdint.h>

#define string(a,b)  char *  a = b;
#define real(a,b)    double  a = b;
#define real3(a,b)   double  a[3] = {0,0,0};
#define integer(a,b) int64_t a = b;

#include "config.template.h"

#endif /* _GRAVITATIONAL_CONSISTENCY_VARS_H_ */
