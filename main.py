import numpy as np
import make_halos as init #NFW_dist(r, rs), cusp_dist(r,rs), NFW_vel(r, rs, M, Rvir), cusp_vel(r,rs,M,Rvir), get_halos(fMH, fDENS, fVELS, rs, delta, npart, M_type, Rv_typ, bigrat, axis)
import methast as mh
import argparse

# parse arguments
parser = argparse.ArgumentParser()

parser.add_argument('profile', type=str,
                    help='string of the profile type: "NFW" or "cusp"')
parser.add_argument('npart', type=int,
                    help='target total number of particles; note: this might change by +/-1')
parser.add_argument('M', type=float,
                    help='typical galactic dark matter halo mass')
parser.add_argument('R', type=float,
                    help='typical galactic dark matter halo virial radius')
parser.add_argument('rs', type=float,
                    help='scale radius for dark matter halo profiles')
parser.add_argument('delta', type=float,
                    help='delta value to use in metropolis hastings')
parser.add_argument('bigrat', type=float,
                    help='ratio of npart in larger halo (ex. 0.5 for equal size halos); note: must be less than 1')
parser.add_argument('sep', type=float,
                    help='separation between halos')
parser.add_argument('axis', type=bool,
                    help='true for sep on axis in plane of rotation, false for sep on axis perpendicular to rotation')

args = parser.parse_args()

#initialize the problem
if args.profile == 'NFW':
    fDENS = init.NFW_dist
    fVELS = init.NFW_vel
elif args.profile == 'cusp':
    fDENS = init.cusp_dist
    fVELS = init.cusp_vel

# get inital conditios
halos_init = init.get_halos(mh.methast, fDENS, fVELS, args.rs, args.delta, args.npart, args.M, args.R, args.bigrat, args.sep, args.axis)
print(halos_init)

# run integrator.py and get all that good stuff
# format output in the way we need it for analysis
# graphs?? plt.plot(physics)
