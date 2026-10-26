import atexit

import imas
import numpy as np
from imas.util import get_data_dictionary_version
import argparse
import getpass
from scipy import interpolate
from scipy.interpolate import griddata

# --- Routine to find index with nearest float value in an array
def find_nearest(arr, target_value):
    idx = abs(arr - target_value).argmin()
    return arr.flat[idx],idx

# --- Routine to calculate Coulmb's log ei
def coulomb_log(Te, ne):  # Te in eV, ne in SI
    if np.ndim(Te)>0:
        mask1 = Te < 10
        mask2 = Te > 10
        Coulomb_Log = np.copy(Te)
        Coulomb_Log[mask1] = 23.0000 - np.log((ne*1e-6)**0.5*Te[mask1]**(-1.5))
        Coulomb_Log[mask2] = 24.1513 - np.log((ne*1e-6)**0.5*Te[mask2]**(-1.0))
    else:
        if (Te<10):
            Coulomb_Log = 23.0000 - np.log((ne*1e-6)**0.5*Te**(-1.5))
        else:
            Coulomb_Log = 24.1513 - np.log((ne*1e-6)**0.5*Te**(-1.0))
    return Coulomb_Log

# --- Routine to calculate Spitzer's resistivity
def eta_spitzer(Te, ne, Zeff):  # Te in eV, ne in SI, eta in SI
    coef_Zeff = Zeff*(1.+1.198*Zeff+0.222*Zeff**2)/(1.+2.966*Zeff+0.753*Zeff**2) / ((1.+1.198+0.222)/(1.+2.966+0.753))
    eta = 1.65e-9 * coulomb_log(Te, ne) * (Te*0.001)**(-1.5) * coef_Zeff # From Wesson
    return eta
  
# --- Routine to calculate Spitzer-Haerm conductivity
def kappa_par_Spitzer(Te):  # Te in eV, ne in SI, kappa_par in SI
    return 3.6e29 * (Te*1e-3)**2.5

# --- Routine to construct the initial boundary of JOREK's grid
def build_JOREK_boundary(tokamak_name):

  R_scale = 1

  if (tokamak_name == 'ITER'):
    #--------------------close fit to ITER wall
    ellip  = 2.0
    tria_u = 0.55
    tria_l = 0.65
    quad_u = -0.1
    quad_l = 0.15
    n_tht  = 257
    r0     = 6.2  * R_scale
    z0     = 0.1  * R_scale
    a0     = 2.25 * R_scale
  
    #-------------------- contour outside ITER wall
    ellip  = 2.1
    tria_u = 0.58
    tria_l = 0.65
    quad_u = -0.12
    quad_l = -0.
    n_tht  = 257
    r0     = 6.2   * R_scale
    z0     = -0.05 * R_scale
    a0     = 2.34  * R_scale
    a0     = 2.30  * R_scale
  
  elif (tokamak_name == 'JET'):
    
    #-------------------- contour outside JET wall
    # blue contour in https://www.jorek.eu/wiki/doku.php?id=eqdsk2jorek.f90
    ellip  = 1.85
    tria_u = 0.4
    tria_l = 0.4
    quad_u = -0.2
    quad_l = -0.2
    n_tht  = 257
    r0     = 2.9  * R_scale
    z0     = 0.1  * R_scale
    a0     = 1.08 * R_scale
  
    #-------------------- contour to avoid too long divertor legs
    # red contour in https://www.jorek.eu/wiki/doku.php?id=eqdsk2jorek.f90
    ellip  = 1.7
    tria_u = 0.4
    tria_l = 0.4
    quad_u = -0.4
    quad_l = -0.2
    n_tht  = 257
    r0     = 2.85 * R_scale
    z0     = 0.15 * R_scale
    a0     = 1.1  * R_scale

  elif (tokamak_name == 'inxflow'):
    
    #-------------------- contour outside JET wall
    # blue contour in https://www.jorek.eu/wiki/doku.php?id=eqdsk2jorek.f90
    ellip  = 1.7
    tria_u = 0.0
    tria_l = 0.0
    quad_u = 0
    quad_l = -0.4
    n_tht  = 257
    r0     = 3.0  * R_scale
    z0     = 0.0  * R_scale
    a0     = 0.98 * R_scale
  
  elif (tokamak_name == 'DIII-D'):
  
    #-------------------- contour outside DIII-D wall
    ellip  = 1.85
    tria_u = 0.4
    tria_l = 0.4
    quad_u = -0.2
    quad_l = -0.2
    n_tht  = 257
    r0     = 1.7 * R_scale
    z0     = 0.  * R_scale
    a0     = 0.7 * R_scale
    
    #-------------------- Atomic physics JOREK/NIMROD/M3D-C1 benchmark case (paper by B. Lyons)
    ellip  = 1.35/0.7
    tria_u = 0.3
    tria_l = 0.3
    quad_u = 0.
    quad_l = 0.
    n_tht  = 257
    r0     = 1.7 * R_scale
    z0     = 0.  * R_scale
    a0     = 0.7 * R_scale
  
  else:
    print( 'Tokamak name not available or wrongly specified, stopping' )
 
  r_bnd = np.zeros( n_tht )
  z_bnd = np.zeros( n_tht )

  n_mid = int(n_tht/2)
 
  for i in range(1, n_mid + 1):
    angle    = 2 * np.pi * float(i-1)/float(n_tht-1)
    r_bnd[i-1] = r0 + a0 * np.cos(angle + tria_u*np.sin(angle) + quad_u*np.sin(2*angle))
    z_bnd[i-1] = z0 + a0 * ellip * np.sin(angle)
  for i in range(n_mid + 1,n_tht+1):
    angle    = 2 * np.pi * float(i-1)/float(n_tht-1)
    r_bnd[i-1] = r0 + a0 * np.cos(angle + tria_l*np.sin(angle) + quad_l*np.sin(2*angle))
    z_bnd[i-1] = z0 + a0 * ellip * np.sin(angle)

  return r_bnd, z_bnd, r0, z0


def parse_dd_version(version):
  """Return the major and minor components of an IMAS DD version."""
  try:
    components = tuple(int(value) for value in version.split('.'))
  except ValueError as error:
    raise argparse.ArgumentTypeError(
        "DD version must use the form MAJOR[.MINOR[.PATCH]]"
    ) from error

  if not components:
    raise argparse.ArgumentTypeError("DD version cannot be empty")

  return components + (0,) * (2 - len(components))


print(" ")
print(" Example of usage: ")
print("    python imas2jorek.py --uri 'imas:hdf5?path=/path/to/imas/data' -t 54.5")
print("    python imas2jorek.py -u public -d ITER -p 105033 -r 1 -t 54.5")
print(" ")
print(" To see options and default values do: ")
print("    python imas2jorek.py -h ")
print(" ")

# Import shot
parser = argparse.ArgumentParser(description="Create a JOREK input file from an equilibrium IDS in a given IMAS database",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--uri", type=str, default=None,
                    help="IMAS data-entry URI. When set, it takes precedence over backend/database/pulse/run/user")
parser.add_argument("-p", "--pulse", type=int, default=1, help="Pulse number")
parser.add_argument("-r", "--run", type=int, default=7, help="Run number")
parser.add_argument("-u", "--user", type=str, default=getpass.getuser(),
                    help="Location of ~$USER/public/imasdb")
parser.add_argument("-d", "--database", type=str, default="test_db", help="Database name under public/imasdb/")
parser.add_argument("-dd", "--dd-version", "--data_dictionary", dest="dd_version", type=str, default=None,
                    help="Optional Data Dictionary version override; by default use the DBEntry version")
parser.add_argument("-o", "--occurrence", type=int, default=0, help="Occurrence number")
parser.add_argument("-tk", "--tokamak", type=str, default="ITER", help="Name of the tokamak (to construct R,Z boundary)")
parser.add_argument("-f", "--backend", type=int, default=imas.ids_defs.HDF5_BACKEND,
                    help="Database format: 12=MDSPLUS, 13=HDF5")
parser.add_argument("-t", "--time", type=float, default=-1, help="The requested time in seconds")
args = parser.parse_args()

entry_options = {} if args.dd_version is None else {"dd_version": args.dd_version}
if args.uri:
  entry = imas.DBEntry(args.uri, "r", **entry_options)
  entry_description = args.uri
else:
  entry = imas.DBEntry(args.backend, args.database, args.pulse, args.run, args.user,
                        **entry_options)
  entry.open()
  entry_description = "legacy entry: user=%s, database=%s, pulse=%s, run=%s" % (
      args.user, args.database, args.pulse, args.run)

atexit.register(lambda: entry.close())
print("Reading IMAS entry: %s" % entry_description)

# Detect the DD version stored in the entry before doing normal reads. Reading
# this metadata with autoconvert=False works even when the stored and default
# DD versions have different major versions.
if args.dd_version is None:
  stored_version_ids = entry.get("equilibrium", occurrence=args.occurrence,
                                 lazy=True, autoconvert=False)
  dd_version_string = get_data_dictionary_version(stored_version_ids)

  # Reopen with the on-disk version so the generated IDS structure matches the
  # data, including entries stored with DD 3.x while DD 4.x is installed.
  if dd_version_string != entry.dd_version:
    entry.close()
    if args.uri:
      entry = imas.DBEntry(args.uri, "r", dd_version=dd_version_string)
    else:
      entry = imas.DBEntry(args.backend, args.database, args.pulse, args.run, args.user,
                           dd_version=dd_version_string)
      entry.open()
else:
  # --dd-version is an explicit target-version override.
  dd_version_string = entry.dd_version

dd_version = parse_dd_version(dd_version_string)
dd_major = dd_version[0]
print("Using Data Dictionary version: %s" % dd_version_string)

# cocos factors
if (dd_major <=3):
  cocos_psi  =  1.0/(2*np.pi)       # Transform to COCOS convention 11 --> 8
else:
  cocos_psi  =  -1.0/(2*np.pi)      # Transform to COCOS convention 17 --> 8

cocos_curr = -1.0 
cocos_Bphi = -1.0
it=0

time = np.asarray(entry.get("equilibrium", occurrence=args.occurrence, lazy=True).time)
equilibrium = entry.get_slice("equilibrium", args.time, imas.ids_defs.CLOSEST_INTERP,
                              occurrence=args.occurrence)
try:
  pf_active = entry.get_slice("pf_active", args.time, imas.ids_defs.CLOSEST_INTERP,
                              occurrence=args.occurrence)
except Exception as error:
  print("Warning: pf_active occurrence %d is unavailable; no coil currents will be imported (%s)" %
        (args.occurrence, error))
  pf_active = None
# Use plasma_profiles if available
if dd_major > 3:
  use_core_prof  = False
  try:
    profiles = entry.get_slice("plasma_profiles", args.time, imas.ids_defs.CLOSEST_INTERP,
                               occurrence=args.occurrence)
    if (len(profiles.profiles_1d) == 0 or
        len(profiles.profiles_1d[it].electrons.temperature) == 0):
      use_core_prof = True
    else:
      print('Use plasma_profiles')
  except Exception as error:
    print('Could not read usable plasma_profiles (%s)' % error)
    use_core_prof = True
else:
  use_core_prof = True
if (use_core_prof):
  print('Use core_profiles')
  profiles = entry.get_slice("core_profiles", args.time, imas.ids_defs.CLOSEST_INTERP,
                             occurrence=args.occurrence)
    
# Find out array index of the requested time
tc   = equilibrium.time[0]

print(' **************** Found time and index *****************')
print(' Time  = %.5fs with t=%.5f - %.5fs '% (tc, time[0], time[-1]))
print(" ")

########## Create boundary of the JOREK domain ###############
boundary_line = build_JOREK_boundary(args.tokamak)
R_bnd = boundary_line[0] 
Z_bnd = boundary_line[1] 
R_geo = boundary_line[2] 
Z_geo = boundary_line[3] 
###############################################################

# Constants
mu0      = 12.566370614359e-7
AMU      = 1.660539040e-27
e_ch     = 1.6021766e-19

# Read 0D parameters
if dd_major < 4:
  a_min      = equilibrium.time_slice[it].boundary_separatrix.minor_radius
else:
  a_min      = equilibrium.time_slice[it].boundary.minor_radius

eps          = a_min / R_geo
B_geo        = equilibrium.vacuum_toroidal_field.r0 * equilibrium.vacuum_toroidal_field.b0[0] / R_geo * cocos_Bphi
xip          = equilibrium.time_slice[it].global_quantities.ip           * cocos_curr
if dd_version >= (4, 1):
  psi_axis     = equilibrium.time_slice[it].global_quantities.psi_magnetic_axis * cocos_psi
  psi_boundary = equilibrium.time_slice[it].boundary.psi                        * cocos_psi
else:
  psi_axis     = equilibrium.time_slice[it].global_quantities.psi_axis     * cocos_psi
  psi_boundary = equilibrium.time_slice[it].global_quantities.psi_boundary * cocos_psi
beta_p       = equilibrium.time_slice[it].global_quantities.beta_pol                
beta_tor     = equilibrium.time_slice[it].global_quantities.beta_tor
if dd_major < 4:
  beta_normal  = equilibrium.time_slice[it].global_quantities.beta_normal
else:
  beta_normal  = equilibrium.time_slice[it].global_quantities.beta_tor_norm
li3          = equilibrium.time_slice[it].global_quantities.li_3                    
volume       = equilibrium.time_slice[it].global_quantities.volume                  
area         = equilibrium.time_slice[it].global_quantities.area                    
q_axis       = equilibrium.time_slice[it].global_quantities.q_axis                  
q_95         = equilibrium.time_slice[it].global_quantities.q_95                    
Wth          = equilibrium.time_slice[it].global_quantities.energy_mhd              
R_axis       = equilibrium.time_slice[it].global_quantities.magnetic_axis.r         
Z_axis       = equilibrium.time_slice[it].global_quantities.magnetic_axis.z         
#R_xpoint     = equilibrium.time_slice[it].boundary_separatrix.x_point.r        
#Z_xpoint     = equilibrium.time_slice[it].boundary_separatrix.x_point.z        

print(' **************** Equilibrium quantities *****************')
print( ' Ip                = %s MA'%(xip* 1e-6) )
print( ' B_geo             = %s T'%(B_geo) )
print( ' a_min             = %s m'%(a_min) )
print( ' area              = %s m^2'%(area) )
print( ' volume            = %s m^3'%(volume) )
print( ' psi_axis, psi_bnd, difference  = %s , %s , %s Wb/rad'%(psi_axis, psi_boundary, psi_boundary-psi_axis) )
print( ' beta_p, beta_tor, beta_norm    = %s , %s , %s '%(beta_p, beta_tor, beta_normal) )
print( ' li3, q_axis, q_95              = %s , %s , %s '%(li3, q_axis, q_95) )
print( ' R_axis, Z_axis                 = %s , %s '%(R_axis, Z_axis) )
#print( ' R_xpoint, Z_xpoint             = %s , %s '%(R_xpoint, Z_xpoint) )
print( ' ')


# Read 1D profiles  
psi_1d      = equilibrium.time_slice[it].profiles_1d.psi                 * cocos_psi
qprof       = equilibrium.time_slice[it].profiles_1d.q     #             * cocos_psi
pressure_1d = equilibrium.time_slice[it].profiles_1d.pressure
ffprime_1d  = equilibrium.time_slice[it].profiles_1d.f_df_dpsi     *(-1.0)      / cocos_psi
ions        = profiles.profiles_1d[it].ion
Te_ids      = profiles.profiles_1d[it].electrons.temperature
psi_1d_core = profiles.profiles_1d[it].grid.psi * cocos_psi

# Get ion density profile
n_tot    = np.zeros( len(ions[0].density))
rho_tot  = np.zeros( len(ions[0].density))

for ion in ions:
    if (len(ion.density)<1):
       continue
    n_tot   += ion.density
    rho_tot += ion.density*ion.element[0].a * AMU

psi_norm    = (psi_1d      - psi_axis) / (psi_boundary - psi_axis)
psi_norm2   = (psi_1d_core - psi_axis) / (psi_boundary - psi_axis)

psi_1d = psi_1d[psi_norm<=1.]
psi_norm = psi_norm[psi_norm<=1.]
# Export 1d profiles to compare with JOREK
np.savetxt('profiles_ids.txt', np.transpose([psi_norm2, Te_ids, n_tot, qprof]   ) )

# Read 2D profiles  
R_2d   = equilibrium.time_slice[it].profiles_2d[0].r
Z_2d   = equilibrium.time_slice[it].profiles_2d[0].z
psi_2d = equilibrium.time_slice[it].profiles_2d[0].psi  * cocos_psi

# Read PF coil currents when the optional pf_active IDS is available
coils  = [] if pf_active is None else pf_active.coil

# Get poloidal flux at the JOREK boundary
Ra = R_2d.flatten()
Za = Z_2d.flatten()
pa = psi_2d.flatten()

psi_bnd = griddata((Ra, Za), pa, (R_bnd, Z_bnd))  

n0              = n_tot[0]
central_density = n0*1e-20
central_mass    = rho_tot[0] / (n0 * AMU)
rho0            = central_mass * AMU * central_density * 1e20 
rho_unit        = rho_tot / rho_tot[0]

# Extend profiles into the SOL
n_sol          = 80
n_core         = len(psi_1d)
n_prof         = n_sol + n_core
psin_sol       = 2.0
psin_sep       = psi_norm[-1] + (psi_norm[-1] - psi_norm[-2])

psi_norm_ext    = np.zeros( n_prof )
ffprime_1d_ext  = np.zeros( n_prof )
pressure_1d_ext = np.zeros( n_prof )
rho_1d_ext      = np.zeros( n_prof )

# Fill core values in
psi_norm_ext   [0:n_core] = psi_norm[0:n_core]
ffprime_1d_ext [0:n_core] = ffprime_1d[0:n_core]
pressure_1d_ext[0:n_core] = pressure_1d[0:n_core]
rho_1d_ext     [0:n_core] = rho_unit [0:n_core]

# Extend psi
psi_norm_ext[n_core:n_prof] = np.linspace(psin_sep,psin_sol,num=n_sol) 

# Fill SOL with last values
ffprime_1d_ext  [n_core:n_prof] = np.full( n_sol,  ffprime_1d[n_core-1] ) 
pressure_1d_ext [n_core:n_prof] = np.full( n_sol, pressure_1d[n_core-1] ) 
rho_1d_ext      [n_core:n_prof] = np.full( n_sol,    rho_unit[n_core-1] ) 

# Ramp-down values in the SOL
ramp_down   = 0.5*(1 - np.tanh((psi_norm_ext - 1.01)/0.01) )

rho_bnd     = rho_unit[-1]
rho_1d_ext  = (rho_1d_ext-rho_bnd) * ramp_down + rho_bnd
#rho_1d_ext  = np.full( n_tot,    1.0 ) 

temp_1d_ext = pressure_1d_ext*mu0 /rho_1d_ext  # p/rho_mass (as specified in the JOREK normalization wiki)
T_bnd       = 2 * (e_ch * mu0 * n0 )           # 2 eVs
temp_1d_ext = (temp_1d_ext-T_bnd) * ramp_down + T_bnd

# Export profiles
np.savetxt('jorek_density',         np.transpose( [psi_norm_ext, rho_1d_ext]     ) )
np.savetxt('jorek_temperature',     np.transpose( [psi_norm_ext, temp_1d_ext]    ) )
np.savetxt('jorek_ffprime',         np.transpose( [psi_norm_ext, ffprime_1d_ext] ) )

# Write namelist files
namelist = open('jorek_namelist', 'w')
  
namelist.write( "***********************************************\n")
namelist.write( "* namelist from imas2jorek.py                 *\n")
if args.uri:
  namelist.write("* URI %s\n" % args.uri)
else:
  namelist.write("* pulse %06i        run %02i          *\n"%(args.pulse,args.run))
  namelist.write("* database %s user %s*\n"%(args.database,args.user))
namelist.write("* time %.6f                           \n*" %tc)
namelist.write( "***********************************************\n")

namelist.write(" &in1"+"\n")

namelist.write("  tstep = 1."+"\n")
namelist.write("  nstep = 0"+"\n")
namelist.write("  nout = 20"+"\n")
namelist.write("\n")
namelist.write("  freeboundary          = .f." +"\n")
namelist.write("  resistive_wall        = .f." +"\n")
namelist.write("  freeboundary_equil    = .f." +"\n")
namelist.write("  wall_resistivity_fact = 1.d0"+"\n")
namelist.write("\n")
namelist.write("  psi_axis_init = %18.9e \n"%(psi_axis))
namelist.write("  amix          = 0.d0"+"\n")
namelist.write("  amix_freeb    = 0.d0"+"\n")
namelist.write("  use_mumps_eq  = .t."+"\n")
namelist.write("\n")
for i in range(0, len(coils)):
    namelist.write("  pf_coils(%3d)%%current = %18.9e  ! %s \n"%(i+1, coils[i].current.data[it]*cocos_curr, coils[i].name))
    
namelist.write("\n")
namelist.write("! --- Grid parameters -----------------------\n")
namelist.write("  n_R      = 0"+"\n")
namelist.write("  n_Z      = 0"+"\n")
namelist.write("  n_radial = 41"+"\n")
namelist.write("  n_pol    = 64"+"\n")
namelist.write("  n_flux   = 0"+"\n")
namelist.write("  n_tht    = 64"+"\n")
namelist.write("  n_open   = 15"+"\n")
namelist.write("  n_leg    = 15"+"\n")
namelist.write("  n_private = 9"+"\n")
namelist.write("  dPSI_open    = 0.04"+"\n")
namelist.write("  dPSI_private = 0.02"+"\n")
namelist.write(" "+"\n")
namelist.write("! --- Initial conditions --------------------\n")
namelist.write("  rho_file     = 'jorek_density'"+"\n")
namelist.write("  T_file       = 'jorek_temperature'"+"\n")
namelist.write("  ffprime_file = 'jorek_ffprime'"+"\n")
namelist.write(" "+"\n")

namelist.write(" "+"\n")
namelist.write("! --- Physical parameters -------------------\n")

eta_JU  = eta_spitzer(Te_ids[0], n_tot[0], 1) * np.sqrt(rho0/mu0)
namelist.write("  keep_current_prof     = .false. ! Fix current profile\n")
namelist.write("  eta                   = %18.9e    !!! Spitzer value at Zeff=1 \n"%(eta_JU))
namelist.write("  eta_ohmic             = %18.9e    !!! Spitzer value at Zeff=1 \n"%(eta_JU))
namelist.write("  eta_T_dependent       = .t.  \n" )
namelist.write("  eta_coul_log_dep      = .t.  \n" )
Te_norm_fact = e_ch * mu0 * n0 
namelist.write("  T_max_eta             = %18.9e  !!! You may want to modify this, taken at T_tot_axis \n"%( Te_ids[0]*Te_norm_fact*2) )
namelist.write("  T_max_eta_ohm         = 1d99    !!! No upper value for ohmic heating \n" )
namelist.write("  T_min                 = %18.9e  !!! T_tot_min = 1eV, double check!! \n"%( 1*Te_norm_fact) )
namelist.write("  T_min_neg             = %18.9e  !!! T_tot_min = 1eV, double check!! \n"%( 1*Te_norm_fact) )
namelist.write("  Te_eV_min             = 1.d0    !!! Used for radiation and impurity functions \n" )
namelist.write("  corr_neg_temp_coef(1) = 1.0     !!! These are used to make sure the minimum corrected temperature is T_min_neg \n" )
namelist.write("  corr_neg_temp_coef(2) = 0.1 \n" )
namelist.write("  rho_min               = %18.9e  !!! double check!! \n"%( rho_1d_ext[-1]*0.5 ) )
namelist.write("  rho_min_neg           = %18.9e  !!! double check!! \n"%( rho_1d_ext[-1]*0.5 ) )

diffusion_coef_default = 1.0   # 1 m2/s
namelist.write("\n")
namelist.write("  F0     = %18.9e \n"%(R_geo * B_geo))
namelist.write("  xpoint = .t."+"\n")
namelist.write("  gamma_stangeby  = 8.d0"+"\n")
namelist.write("  bc_natural_open = .t."+"\n")
namelist.write("\n")
namelist.write("  D_par     = %18.9e    !!! 1m2/s - Double check! \n"%(diffusion_coef_default * np.sqrt(rho0*mu0)))
namelist.write("  D_perp(1) = %18.9e    !!! 1m2/s - Double check! \n"%(diffusion_coef_default * np.sqrt(rho0*mu0)))
namelist.write("  D_perp(2) = 0.85d0"+"\n")
namelist.write("  D_perp(3) = 0.d0"+"\n")
namelist.write("  D_perp(4) = 0.01d0"+"\n")
namelist.write("  D_perp(5) = 1.d4"+"\n")
namelist.write("\n")
fact_norm_kappa = np.sqrt(mu0/rho0) * (5.0/3.0-1.0) * central_mass * AMU 
Kappa_par_JU    = kappa_par_Spitzer(Te_ids[0]) * fact_norm_kappa
namelist.write("  ZKpar_T_dependent = .t. "+"\n")
namelist.write("  ZK_par            = %18.9e    !!! Spitzer value in the core \n"%(Kappa_par_JU))
namelist.write("  ZK_par_max        = %18.9e    !!! double check!! you may want to lower if Te_axis > 1keV \n"%(Kappa_par_JU))
namelist.write("  T_min_ZKpar       = %18.9e    !!! set to 1 eV, you may want to rise it to resolve sharp gradients \n"%(Te_norm_fact*1))
namelist.write("\n")
namelist.write("  ZK_perp(1) = %18.9e    !!! 1m2/s - Double check! \n"%(diffusion_coef_default * fact_norm_kappa * n0 ))
namelist.write("  ZK_perp(2) = 0.85d0"+"\n")
namelist.write("  ZK_perp(3) = 0.d0"+"\n")
namelist.write("  ZK_perp(4) = 0.01d0"+"\n")
namelist.write("  ZK_perp(5) = 1.d4"+"\n")
namelist.write("\n")
visco_JU = diffusion_coef_default * np.sqrt(mu0/rho0) * rho0
namelist.write("  visco_T_dependent = .f. "+"\n")
namelist.write("  visco             = %18.9e    !!! 1 m2/s - Double check! \n"%(visco_JU))
namelist.write("  visco_heating     = %18.9e    !!! 1 m2/s - Double check! \n"%(visco_JU))
namelist.write("  visco_par         = %18.9e    !!! 10m2/s - Double check! \n"%(visco_JU*10))
namelist.write("  visco_par_heating = %18.9e    !!! 10m2/s - Double check! \n"%(visco_JU*10))
namelist.write(""+"\n")
namelist.write("  heatsource     = 0.d0"+"\n")
namelist.write("  particlesource = 0.d0"+"\n")
namelist.write(""+"\n")
namelist.write("  central_mass    = %18.9e \n"%(central_mass))
namelist.write("  central_density = %18.9e \n"%(central_density))
namelist.write(""+"\n")
namelist.write("  R_geo = %18.9e \n"%(R_geo))
namelist.write("  Z_geo = %18.9e \n"%(Z_geo))
namelist.write(""+"\n")
namelist.write("! --- Initial psi ------------------------\n")
namelist.write("  mf = 0 \n")
namelist.write("  n_boundary = "+str(len(R_bnd))+"\n")

for i in range(0, len(R_bnd)):
    namelist.write("  R_boundary(%3d)= %18.9e, Z_boundary(%3d)= %18.9e, psi_boundary(%3d)= %18.9e \n" %(i+1, R_bnd[i], i+1, Z_bnd[i], i+1, psi_bnd[i]) )

namelist.write( "/"+"\n")
namelist.close()
