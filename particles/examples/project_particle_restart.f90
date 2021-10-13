program project_particle_restart
use particle_tracer
use mod_particle_io
implicit none

type(projection)                 :: proj
type(event)                      :: fieldreader

! Start up MPI, jorek
call sim%initialize(num_groups=1)

! Read a sim
call read_simulation_hdf5(sim, 'part_restart.h5')

! Set up the field reader
fieldreader = event(read_jorek_fields_interp_linear(basename='jorek', i=-1))
call with(sim, fieldreader)

proj = new_projection(sim%fields%node_list, sim%fields%element_list, filter=1d-5, filter_hyper=1d-10, &
                      f=[proj_f(proj_one, group=1), &
                         proj_f(proj_q, group=1)], &
                      to_h5=.true.,basename='proj')

call with(sim, proj)
end program project_particle_restart
