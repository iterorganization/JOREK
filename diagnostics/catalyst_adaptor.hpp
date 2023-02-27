#ifndef CATALYST_ADAPTOR_H_
#define CATALYST_ADAPTOR_H_

#if USE_CATALYST

extern "C" {
void catalyst_adaptor();

void catalyst_adaptor_initialise(char *a_catalyst_scripts);

void catalyst_adaptor_execute(int *a_step_index, double *a_time);

void catalyst_adaptor_finalise();

// ********** Fortran functions below *********

void catalyst_get_params(int *a_nsub, int *a_n_elements, int *a_n_scalars);

void catalyst_get_scalar_name(char *a_scalar_name, int *a_iscalar);

void catalyst_interp_grid(int *a_nnos, int *a_nel, float *a_coords_R,
                          float *a_coords_Z, int *a_cell_points);

void catalyst_interp_scalar(float *a_scalars, int *a_iscalar);
}

#endif /* USE_CATALYST */
#endif