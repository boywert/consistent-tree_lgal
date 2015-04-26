#ifndef _MASSES_H_
#define _MASSES_H_
#include <math.h>

#ifndef M_PI
#define M_PI        3.14159265358979323846264338327950288   /* From math.h */
#endif /* M_PI */

extern double Om, Ol;

double kravtsov_f(double x);
double delta_vir(double a);
double calculate_mass_factor(double mvir, double rvir, double rs);
void gen_ff_cache(void);
double concentration(double mvir, double scale);

/* Returns the filling fraction as a function of r_h / r_s */
float ff_cached(float x);
void convert_mvir_to_delta(float mvir, float rvir, float rs, float scale,
			   float delta, char type,
			   float *mdelta, float *rdelta);
#endif /* _MASSES_H_ */
