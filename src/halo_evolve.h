#ifndef HALO_EVOLVE_LIB_H
#define HALO_EVOLVE_LIB_H

#define VEL_RESOLUTION 0.5 /* In km/s */
#define NUM_TIMESTEPS 10

//Newton's gravitational constant 
/* Google query: "G in Mpc^2 / 1e6 yr / (mass of the sun) * km/s" */
#define Gc 4.39877036e-15 /* In Mpc^2 / Myr / Msun * km/s */


void evolve_halos(float a1, float a2, struct halo_stash *h);

inline float halo_grav_range(float mvir, float dt);

/* Advances simulation from a2 to a3; old a1 needed for velocity updates. */
/* special_step = 1 (first step), 0 (normal step), -1 (last step) */
void do_timestep(double a1, double a2, double a3, int special_step);

void print_evolved_halos(void);


#endif /* HALO_EVOLVE_LIB_H */
