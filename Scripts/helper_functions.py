import os
import tfs
import numpy as np
from math import log10, floor

import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

########################################################################
# Matplotlib global plotting parameters
########################################################################
plt.rcParams['figure.figsize'] = [8.0, 5.0]
plt.rcParams['figure.dpi'] = 200
plt.rcParams['savefig.dpi'] = 200

plt.rcParams['axes.titlesize'] = 14
plt.rcParams['axes.labelsize'] = 14

plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10

plt.rcParams['font.size'] = 10
plt.rcParams['legend.fontsize'] = 8

plt.rcParams['lines.linewidth'] = 1.5
plt.rcParams['lines.markersize'] = 5

########################################################################
# Used to calculate initial beam parameters for a given proton energy
########################################################################
class MADX_Proton_Beam_Parameters:                
	mass = 938.272E6 # in eV
	energy = -1. # in eV   
	beta = -1.
	gamma = -1.
	total_energy = -1.
	momentum = -1.
	
	def __init__(self, energy):
		self.energy = energy
		self.total_energy = self.get_total_energy()
		self.gamma = self.get_gamma()
		self.beta = self.get_beta()
		self.momentum = self.get_momentum()
	
	def get_total_energy(self): return (self.energy + self.mass)
	def get_gamma(self): return (self.total_energy / self.mass)
	def get_beta(self): return(np.sqrt( 1 - (1/self.gamma)**2 ))
	def get_momentum(self): return(self.gamma * self.mass * self.beta)
	
	def print_beam(self):
		print('M_proton = ', round_sig(self.mass/1E6) , 'MeV')
		print('Energy = ', round_sig(self.energy/1E9) , 'GeV')
		print('Total Energy = ', round_sig(self.total_energy/1.E9), 'GeV')
		print('Gamma = ', round_sig(self.gamma))
		print('Beta = ', round_sig(self.beta))
		print('Momentum = ', round_sig(self.momentum/1E9, 8), 'GeV/c')


########################################################################
# Check if file exists
########################################################################
def check_if_file_exists(name):
    ret_val = False
    if os.path.isfile(name):
        print (name, ' already exists')
        ret_val = True
    return ret_val
    
########################################################################
# Make directory
########################################################################  
def make_directory(path):
    if os.path.isdir(path):
        print ("Directory %s already exists" % path)  
    else:
        try:
            os.mkdir(path)
        except OSError:
            print ("Creation of the directory %s failed" % path)
        else:
            print ("Successfully created the directory %s" % path)  

########################################################################
# Relativistic functions
########################################################################             
        
def LorentzGamma(E_tot, E_rest=938.27208816E6):
    return (E_tot / E_rest)
    
def LorentzGamma_from_beta(beta):
    return (1./np.sqrt(1.-beta**2))    

def LorentzBeta(gamma):
    return np.sqrt( 1. - (1./gamma**2) )

def RelativisticMomentum(gamma, E_rest=938.27208816E6):
    return (gamma * E_rest * LorentzBeta(gamma))

    
def E_from_gamma(gamma, E_rest=938.27208816E6):
    return (gamma*E_rest)

########################################################################
# Delta P over P from dE or vice versa
# dp_ov_p = dE_ov_E/beta^2
########################################################################
def dpp_from_dE(dE, E, beta):
    return (dE / (E * beta**2))
    
def dE_from_dpp(dpp, E, beta):
    return (dpp * E * beta**2)

def z_to_time(z, beta): 
    c = 299792458
    return z / (c * beta)
    
########################################################################
# Round number to n significant figures
########################################################################    
def round_sig(x, sig=3):
    return round(x, sig-int(floor(log10(abs(x))))-1)      

########################################################################
# Read PTC Twiss and return dictionary of columns/values
########################################################################
def Read_PTC_Twiss_Return_Dict(filename, verbose=True):
    # Dictionary for output
    d = dict()
    d['HEADER_FILENAME'] = filename
    keywords = ''
    
    # First we open and count header lines
    fin0=open(filename,'r').readlines()
    headerlines = 0
    for l in fin0:
        # Store each header line
        headerlines = headerlines + 1
        # Stop if we find the line starting '* NAME'
        if '* NAME' in l:
            keywords = l
            break
        # Store the headers as d['HEADER_<name>'] = <value>
        else:
            if '"' in l:
                d[str('HEADER_'+l.split()[1])]=[str(l.split('"')[1])]
            else:
                d[str('HEADER_'+l.split()[1])]=[float(l.split()[-1])]                 
    headerlines = headerlines + 1    
    
    if verbose: print ('\nRead_PTC_Twiss_Return_Dict found Keywords: \n',keywords)
    
    # Make a list of column keywords to return (as an aid to iterating)
    dict_keys = []
    for key in keywords.split():
        dict_keys.append(key)
    dict_keys.remove('*')
    
    if verbose: print ('\nRead_PTC_Twiss_Return_Dict Dict Keys: \n',dict_keys)
    
    # Initialise empty dictionary entries for column keywords 
    for key in dict_keys:
        d[key]=[]
        
    if verbose: print ('\nRead_PTC_Twiss_Return_Dict header only dictionary \n', d)
    
    # Strip header
    fin1=open(filename,'r').readlines()[headerlines:]   
    
    # Populate the dictionary line by line
    for l in fin1:
        i = -1        
        for value in l.split():
            i = i+1
            if 'KEYWORD' or 'NAME' in dict_keys[i]:
                d[dict_keys[i]].append(str(value))
            else:
                d[dict_keys[i]].append(float(value))    
                
    # Return list of column keywords 'dict_keys', and dictionary 'd'
    return dict_keys, d
   
# Adapted from Alex Huschauer's (CERN ABP) webtools found at:
# https://gitlab.cern.ch/acc-models/acc-models-ps/-/tree/2021/_scripts/web    
def plot_lattice_elements(ax, twiss_in, suppress_x=True):

    # Set plot limits
    ax.set_ylim(-1.5,1.5)
    
    # Suppress ticks on y-axis
    #ax.set_yticks([])
    #ax.set_yticklabels([])
    ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
        
    if suppress_x:
        ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        #ax.set_xticks([])
        #ax.set_xticklabels([])

    # Extract Twiss Header
    twiss = tfs.read(twiss_in)
    twissHeader = dict(twiss.headers)
    print('plot_lattice_elements for sequence ', twissHeader['SEQUENCE'])
    
    # Positions and lengths of elements
    pos = twiss.S.values - twiss.L.values/2
    lengths = twiss.L.values
    total_length = (pos[-1]+lengths[-1])
    print('Full length of accelerator lattice = ', total_length, 'm')
    
    # modify lengths in order to plot zero-length elements
    lengths[np.where(lengths == 0)[0]] += 0.001
    
    # Plot line through centre
    ax.plot([0, total_length], [0., 0.], color='grey', linestyle='-', linewidth=0.5)
    
    # Markers - black 0.1m centred lines    
    idx = np.array([idx for idx, elem in enumerate(twiss.KEYWORD.values) if 'MARKER' in elem])
    for i in idx:
        ax.add_patch(Rectangle((pos[i], -0.5), width = 0.1, height = 1., angle=0.0, ec='k', fc='k', lw=0.0))
    
    # BENDS - blue centred rectangles   
    idx = np.array([idx for idx, elem in enumerate(twiss.KEYWORD.values) if 'BEND' in elem])
    for i in idx:
        ax.add_patch(Rectangle((pos[i], -0.5), width = lengths[i], height = 1., angle=0.0, ec='k', fc='b', lw=0.0))
    
    # Kickers - cyan centred rectangles 
    idx = np.array([idx for idx, elem in enumerate(twiss.KEYWORD.values) if 'HKICKER' in elem])
    for i in idx:
        ax.add_patch(Rectangle((pos[i], -0.5), width = 0.1, height = 1., angle=0.0, ec='k', fc='c', lw=0.0))
    idx = np.array([idx for idx, elem in enumerate(twiss.KEYWORD.values) if 'VKICKER' in elem])       
    for i in idx:
        ax.add_patch(Rectangle((pos[i], -0.5), width = 0.1, height = 1., angle=0.0, ec='k', fc='c', lw=0.0))
              
    # QUADRUPOLES - red offset rectangles indicating Focussing or Defocussing
    idx = np.array([idx for idx, elem in enumerate(twiss.KEYWORD.values) if 'QUADRUPOLE' in elem])
    name = np.array(twiss.NAME.values)[idx]
    if (twissHeader['SEQUENCE'] == 'SYNCHROTRON'):        
        idx_1 = idx[np.array([i for i, n in enumerate(name) if 'QD' in n])]
        idx_2 = idx[np.array([i for i, n in enumerate(name) if 'QF' in n])]
        
    elif (twissHeader['SEQUENCE'] == 'RING'):    
        idx_1 = idx[np.array([i for i, n in enumerate(name) if (n.startswith('D') or 'QTD' in n)])]
        idx_2 = idx[np.array([i for i, n in enumerate(name) if (n.startswith('F') or 'QTF' in n)])]
        
    else: 
        idx_1 = idx[np.array([i for i, n in enumerate(name) if 'D' in n ])]
        idx_2 = idx[np.array([i for i, n in enumerate(name) if 'F' in n ])]        
        
    offset = [-0.5, 0.5]
    for i in idx_1:
        ax.add_patch(Rectangle((pos[i], (-0.5 + offset[0])), width = lengths[i], height = 1., angle=0.0, ec='k', fc='r', lw=0.0))
    for i in idx_2:
        ax.add_patch(Rectangle((pos[i], (-0.5 + offset[1])), width = lengths[i], height = 1., angle=0.0, ec='k', fc='r', lw=0.0))
   
