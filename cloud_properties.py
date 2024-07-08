#import sys
#import os
import numpy as np
#import matplotlib.pyplot as plt
#from astropy.io import fits
#from astropy.wcs import WCS
import argparse
#from tqdm import tqdm
import time

#######################################
# Create command line argument parser
#######################################

def create_parser():

	# Handle user input with argparse
    parser = argparse.ArgumentParser(
        description="Compute thermal velocity.")

    parser.add_argument('-T', '--temperature',
        dest='T',
        default=1.0e4,
        type=float,
        help='Gas temperature in Kelvin.')

    parser.add_argument('-m', '--mass',
        dest='m',
        default=1.67355675e-24,
        type=float,
        help='Particle mass in grams (default H).')

    parser.add_argument('-mu', '--mu',
        dest='mu',
        default=1.0,
        type=float,
        help='Molecular weight (default 1).')

    parser.add_argument('-g', '--gamma',
        dest='gamma',
        default=1.667,
        type=float,
        help='Adiabatic index (default 5/3).')

    parser.add_argument('-R',
        dest='R',
        default=1.0,
        type=float,
        help='Radius in pc (default R=1).')

    parser.add_argument('-n',
        dest='n',
        default=1.0,
        type=float,
        help='Number density [cm^-3] (default n=1).')

    parser.add_argument('-v', '--verbose',
        dest='verbose',
        action='store_true',
        help='Print helpful information to the screen? (default: False)',
        default=False)

    return parser

#######################################
# main() function
#######################################
def main():

    #begin timer
    time_global_start = time.time()

    #create the command line argument parser
    parser = create_parser()

    #store the command line arguments
    args   = parser.parse_args()

    #boltzmann in cgs
    kb = 1.3807e-16 #cm^2 g s^-2 K

    #define thermal velocity as sqrt total velocity
    vth = (3*kb*args.T/(args.mu*args.m))**0.5
    vth /=1.0e5 #convert cm/s to km/s

    #define sound speed cs^2 = \gamma P/rho
    cs = (args.gamma*kb*args.T/(args.mu*args.m))**0.5
    cs /=1.0e5 #convert cm/s to km/s

    #gravitational constant
    G = 6.67430e-8 #cgs -- dyn cm^2 * g

    #year in seconds
    year_in_seconds = 3.15576e7 #year in seconds
    myr_in_seconds = year_in_seconds *1e6

    #pc in km
    pc_in_km = 30856775812800

    #mass density
    rho = args.m*args.mu*args.n

    #free fall time in Myr
    tff = 1./np.sqrt(G*rho)
    tff /=myr_in_seconds

    #sound crossing time
    ts = (args.R * pc_in_km)/cs 
    ts /= myr_in_seconds

    #column density
    NH = 2*args.n*args.R*(pc_in_km*1e5)

    #cloud volume in cm
    Vcloud = 4.0*np.pi/3. * (args.R*pc_in_km*1e5)**3

    #solar mass in grams
    msun = 1.9884e33

    #cloud mass in solar masses
    mcloud = Vcloud*args.n*args.m*args.mu / msun

    print(f'Temperature          = {args.T:5.4e} K.')
    print(f'Particle mass        = {args.m:5.4e} g.')
    print(f'Molecular weight     = {args.mu:5.4e}.')
    print(f'Thermal velocity vth = {vth:5.4e} [km/s].')
    print(f'Adiabatic index gam  = {args.gamma:5.4f}.')
    print(f'Ideal sound speed cs = {cs:5.4e} [km/s].')
    print(f'Mass density rho     = {rho:5.4e} [g/cm^3].')
    print(f'Free fall time tff   = {tff:5.4e} [Myr].')
    print(f'Radius               = {args.R:5.4e} [pc].')
    print(f'Sound crossing time  = {ts:5.4e} [Myr].')
    print(f'Column density       = {np.log10(NH):5.4e} [log cm^-2].')
    print(f'Cloud mass           = {mcloud:5.4e} [Solar masses].')

    #end timer
    time_global_end = time.time()
    if(args.verbose):
    	print(f"Time to execute program: {time_global_end-time_global_start}s.")

#######################################
# Run the program
#######################################
if __name__=="__main__":
	main()