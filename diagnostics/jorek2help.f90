!> Command line function, print help or hardcoded parameter values to screen
!!
!! NOTE: To remove this function a -DNO_HELP flag should be added to the
!! compiler flags
!!
!! When one or more commandline argument is given this function will loop
!! through those options and will respond to each one until one is recognized.
!! The following argument are supported:
!!
!!   -h or --help: print some basic info about JOREK_2
!!   -p or --param: print hardcoded parameters
!!   -v or --version: print the JOREK version info
!!   -a or --about: print information about JOREK
!!
!! unknown argument give a short message that the argument is invalid. When
!! none of the given arguments is recognized a sugestion is given for valid
!! arguments.
!!
!! writer: Mark Verbeek
!! data  : 11-02-2016
subroutine jorek2help(n_cpu, nbthreads)
  use constants
  use mod_parameters
  use mod_log_params

#include "version.h"
  
  implicit none

  integer, intent(in) :: n_cpu, nbthreads
  
#ifndef NO_HELP

  integer           :: narg, cptArg !< for commandline arguments
  character(len=20) :: ArgName      !< Argument name

  narg = command_argument_count() ! Get number of commandline agruments
    
  if ( narg > 0 )then
    do cptArg = 1, narg
      call get_command_argument(cptArg,ArgName)
      select case(adjustl(ArgName))
        case("--about","-a")
          call print_about()
        case("--help","-h")
          call print_help()
        case("--version","-v")
          call print_version()
        case("--param","-p")
          call log_parameters(0, .true.)
        case default
          write(*,'(A, A, A)') 'Option ', adjustl(ArgName), 'unknown'
          call print_help()
          stop
      end select
    end do
    stop
  else
    call print_version()
  end if

#endif

return

contains

  subroutine print_about()
    200 format(79('*'))
    write(*,200)
    write(*,'(A)') '*                            JOREK ' // trim(JOREK_VERSION) // '                                 *'
    write(*,'(A)') '*                            https://www.jorek.eu/                            *'
    write(*,200)
    write(*,'(A)') '* Solves the (reduced) MHD equations in 3D toroidal geometry                  *'
    write(*,'(A)') '*                                                                             *'
    write(*,'(A)') '* - solvers implemented:                                                      *'
    write(*,'(A)') '*   - MUMPS                                                                   *'
    write(*,'(A)') '*   - PastiX                                                                  *'
    write(*,'(A)') '*   - STRUMPACK                                                               *'
    write(*,'(A)') '*   - GMRES or BICGSTAB (+direct solver in the preconditioner)                *'
    write(*,'(A)') '*                                                                             *'
    write(*,'(A)') '* - required libraries :                                                      *'
    write(*,'(A)') '*   - MPI                                                                     *'
    write(*,'(A)') '*   - HDF5                                                                    *'
    write(*,'(A)') '*   - MUMPS and/or PaStiX and/or STRUMPACK                                    *'
    write(*,'(A)') '*   - SCOTCH (metis)                                                          *'
    write(*,'(A)') '*   - FFTW                                                                    *'
    write(*,'(A)') '*   - SCALAPACK (BLACS)                                                       *'
    write(*,'(A)') '*   - LAPACK, BLAS                                                            *'
    write(*,'(A)') '*   - PPPLIB                                                                  *'
    write(*,'(A)') '*                                                                             *'
    write(*,'(A)') '* Original author : Guido Huysmans (Euratom / CEA Association)                *'
    write(*,'(A)') '*          authors: various, see https://www.jorek.eu                         *'
    write(*,'(A)') '* start date: 18-7-2008                                                       *'
    write(*,'(A)') '*                                                                             *'
    write(*,200)
  end subroutine print_about

  subroutine print_help()
    200 format(79('*'))
    write(*,200)
    write(*,'(A)') '* List of command line argument options                                       *'
    write(*,'(A)') '*   -a or --about for info about JOREK in general                             *'
    write(*,'(A)') '*   -p or --param for JOREK hardcoded parameters                              *'
    write(*,'(A)') '*   -v or --version for JOREK compile version                                 *'
    write(*,200)
  end subroutine print_help

  subroutine print_version()
    111 format(2x,a,': ',a)
    write(*,*) '*************************************************'
    write(*,'(A)') '*                 JOREK ' // trim(JOREK_VERSION) // '               *'
    write(*,*) '*                https://www.jorek.eu/          *'
    write(*,*) '*************************************************'
    write(*,*) ' MPI processes       : ', n_cpu
    write(*,*) ' OpenMP threads      : ', nbthreads
    write(*,*) ' GIT revision        : ', trim(adjustl(RCS_VERSION))
    write(*,*) ' GIT revision label  : ', trim(adjustl(RCS_LABEL))
    write(*,*) ' GIT revision time   : ', trim(adjustl(RCS_TIME))
    write(*,111) 'compile_user        ', trim(adjustl(compile_user))
    write(*,111) 'compile_machine     ', trim(adjustl(compile_machine))
    write(*,111) 'compile_dir         ', trim(adjustl(compile_dir))
    write(*,111) 'compile_command     ', trim(adjustl(compile_command))
    write(*,111) 'compile_flags       ', trim(adjustl(compile_flags))
    write(*,111) 'compile_includes    ', trim(adjustl(compile_includes))
    write(*,111) 'compile_defines     ', trim(adjustl(compile_defines))
    write(*,111) 'compile_libs        ', trim(adjustl(compile_libs))
    write(*,111) 'compile_modules     ', trim(adjustl(compile_modules))
  end subroutine print_version

end subroutine jorek2help
