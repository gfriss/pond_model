#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import math
import matplotlib.pyplot as plt
import sys


# temperature
T_celsius = 65
T = 65 + 273.15

HCN_mass_influx_ox = []
HCN_mass_influx_red = []
t = []
with open('datafile') as f:
    # assigning the HCN rainout rate and time array
    ...

H2CO_mass_influx_red = 6.01494508e-12
H2CO_mass_influx_ox = 1.15285175e-09#4.11794768e-09
# UV flux
F = 0.4 # W/m^2


# Calculate delta_t
delta_t = t[2] - t[1] # have to check if there's data on this.....


plt.clf()

#Variable Declarations
w_i = 60.7e-9
m_dot_I = 6e8
f_s = 0.32
r = 40.
rho = 2185.
r_p = 1.
A_p = math.pi*r_p**2
r_g = 500.
tau_d_1cm = 4.9e-3
tau_d_5cm = 0.12
tau_d_10cm = 0.48
R_plus = 6371000
gamma = 31557600
rho_w = 1000.

tau_s = 1.
P = {'IDN': 4.5, 'COL': 6., 'CAM': 3.5}
delta = {'IDN': 0.2, 'COL':0.5, 'CAM':0.5}
sp = {'IDN': 0.85, 'COL': 0.3, 'CAM':0.3}
min_water = 0.001 # 1mm
S = 0.95 #Seepage rate 0.36525
Phi = 1e-4
lambda_uv = 225e-9
h = 6.626e-34
c = 2.9979e8
N_A = 6.022e23 #Avogadro's number

#Organic Haze Values
#Trainer et al. (2006) 1e13 to 3e15 g/yr
#m_dot_haze_5 = 522.e-6/(6/365.25)  # 522 mg / 6 days to kg/yr our experiments
#m_dot_haze_05 = 72.e-6/(6/365.25)
#m_dot_haze_5_min = 1.e11*7.
m_dot_haze_05_min = 4.9e11
#m_dot_haze_5_max = 1.e11*7.
m_dot_haze_05_max = 8.9e14

#initialising species (going to be dict in dict)
species = ['HCN', 'adenine', 'guanine', 'uracil', 'cytosine', 'thymine', '2amino', 'ribose', 'formaldehyde', 'xanthine', 'hypoxanthine']

#and mu values in kg/mol
# 'HCN', 'adenine', 'guanine', 'uracil', 'cytosine', 'thymine', '2amino', 'ribose', 'formaldehyde', 'xanthine', 'hypoxanthine'
MU = [0.0270253, 0.13513, 0.15113, 0.1120868, 0.1111, 0.1261133, 0.084077, 0.15013, 0.030031, 0.15211, 0.1361115]

#and densities to calc d 
rho_species = [687, 1470, 2200, 1320, 1550, 1230, 800, 1200, 815, 1600, 2000]
def calc_d(mu_val, rho_val):
    ''' Calculates molecular density(?)'''
    return 2*(3*mu_val/(4*math.pi*N_A*rho_val))**(1./3)

# haze thingies, missing is 0.
c_05_min = np.array([0., 4.6e-6, 0., 0., 0., 0., 0., 0., 0., 0., 0.])
c_05_max = np.array([0., 6.3e-6, 17.5e-6, 3.46e-6, 4.e-6, 1.7e-6, 0., 0., 0., 65e-6, 21.4e-6])
def calc_m_dot_haze(min_or_max, c_spec):
    ''' Calculates the haze mass changing rate'''
    return min_or_max * c_spec

# photo dissociation rates (kg/yr/m^2)
def photodis_rate(muval):
    ''' Calculates the photodissociation rate for a given molecule'''
    return ((Phi*F*lambda_uv*gamma*muval)/(h*c*N_A))

# k values, eq part for non-nucleobase is missing so k is set to 0 for such species
k = [0., 10**(-5902/T + 8.15), 10**(-6330/T + 9.40), 10**(-7649/T + 11.76), 10**(-5620/T + 8.69), 10**(-7709/T + 11.24), 0., 0., 0., 10**(-6230/T + 9.42), 10**(-5270/T + 7.95)]


def calc_m_dot_red(iteration, spec, mu_val, mu_val_hcn):
    ''' Calculates the mass increase due to the incoming mass flux of rain out of HCN and formaldehyde
        in the reducing model'''
    if spec == 'HCN':
        return HCN_mass_influx_red[iteration] * 4 * math.pi * R_plus**2
    elif spec == 'formaldehyde':
        return H2CO_mass_influx_red[iteration] * 4 *math.pi * R_plus**2
    else:
        return HCN_mass_influx_red[iteration] * 4 * math.pi * R_plus**2 * mu_val/mu_val_hcn
    
def calc_m_dot_ox(iteration, spec, mu_val, mu_val_hcn):
    ''' Calculates the mass increase due to the incoming mass flux of rain out of HCN and formaldehyde
        in the oxidising model'''
    if spec == 'HCN':
        return HCN_mass_influx_ox[iteration] * 4 * math.pi * R_plus**2
    elif spec == 'formaldehyde':
        return H2CO_mass_influx_ox[iteration] * 4 *math.pi * R_plus**2
    else:
        return HCN_mass_influx_ox[iteration] * 4 * math.pi * R_plus**2 * mu_val/mu_val_hcn

# initialising dict in dict to store the all the info
constants_and_rates = {}
for spec,mu_spec,rho_spec,c_min,c_max,k_spec in zip(species,MU,rho_species,c_05_min,c_05_max,k):
    constants_and_rates[spec] = {'mu': mu_spec,
                       'd': calc_d(mu_spec, rho_spec),
                       'rho': rho_spec,
                       'm_dot_min': calc_m_dot_haze(m_dot_haze_05_min, c_min),
                       'm_dot_max': calc_m_dot_haze(m_dot_haze_05_max, c_max),
                       'M_uv_dot': photodis_rate(mu_spec),
                       'm_dot_red': calc_m_dot_red(spec, mu_spec, MU[0]),
                       'm_dot_ox': calc_m_dot_ox(spec, mu_spec, MU[0]),
                       'k': k_spec}
    for i in range(0, nt-1):
        ...
    

#Fraction of surviving organics during entry
f_s = {'IDP': 0.06, 'Met': 0.32}

E = S-0.12 + 0.06*T_celsius
tmax = 8 #years
level = 16

nt = (2**level) + 1 #Choosing nt to have twice as many grid points as nx
    

#Constant seepage mass per year
m_seepage_rate = math.pi*rho_w*r_p**2*S

m_i0 = (4./3)*w_i*f_s['Met']*r**3*rho*A_p/r_g**2

pause_Met = {'IDN': 0, 'COL': 0}

#initialising, original mixing a lot oxidising for some reason and last two nucleotides...now included
values = {'IDN':{}, 'COL':{}}
names = ['L_IDP', 'L_Met', 'm_IDP', 'm_Met', 'm_IDP_A', 'm_Met_A', 'C_IDP', 'C_Met'] +\
['m_' + spec + '_red' for spec in species] + ['m_' + spec + '_ox' for spec in species] + ['C_' + spec + '_red' for spec in species] + ['C_' + spec + '_ox' for spec in species] +\
['m_' + spec + '_haze05_min' for spec in species] + ['m_' + spec + '_haze05_max' for spec in species] + ['C_' + spec + '_haze05_min' for spec in species] + ['C_' + spec + '_haze05_max' for spec in species]
for prec_mod in values:
    for name in names:
        values[prec_mod][name] = np.zeros(shape=nt)
    values[prec_mod]['L_IDP'][0] = r_p - min_water
    values[prec_mod]['L_Met'][0] = r_p - min_water
    values[prec_mod]['m_IDP'][0] = math.pi*rho_w*r_p**2*(min_water) #og was r_p-(r_p-min_water)....
    values[prec_mod]['m_Met'][0] = math.pi*rho_w*r_p**2*(min_water)
    

def precipitation(model, iteration):
    ''' Calculates the precipitation rate that is sinusoidal to represent the seasonal cycle for the current iteration.
        P (mean precipitation rate) and sp (seasonal phase shift) are model dependent.'''
    return (delta_t*P[model])*(1 + delta[model]*np.sin(2*math.pi*(t[iteration] - sp[model])/tau_s))

def water_decrease(impactor, model, iteration):
    ''' Calculates the rate of water decrese by adding up evaporation and seepage then subtraction the precipitation.
        This is also done for current iteration steps and all precipitation models.'''
    return E*delta_t + values[model]['L_'+impactor][iteration] - precipitation(model, iteration)

def nucleobase_outflow(iteration, pause, kval, model):
    ''' Calculates the nucleobase outflow rate from meteorites for current iteration step. This only happens when the pond is wet, hence the introduction of the variable pause.
        In case of IDP this is set to zero as outflow from small IDPs are basically instantenaous compared to the simulation time.'''
    if model == 'IDP':
        return 0.
    elif model == 'Met':
        return delta_t * m_i0 * np.e**(-t[iteration-pause]*(gamma*kval + (1./tau_d_1cm))) / tau_d_1cm
    
def water_mass(model, impactor, iteration):
    ''' Calculates the current water mass.'''
    return math.pi*rho_w*r_p**2*(r_p-values[model]['L_'+impactor][iteration])

def mass_increase(mdot, model = 'IDP', wi = 1., fs = 1.):
    ''' Calculates the mass increase using the timestep (delta_t) and unit pond surface area (A_p/(4pi*R_plus**2))
        This term is not present when working with Meteorites. For generalisation porpuses, though, it appaers there but is set to 0.'''
    if model == 'Met':
        return 0
    elif model == 'IDP':
        return (delta_t * wi * mdot * fs * A_p) / (4*np.pi * R_plus**2)
    
def photo_destruction(muvdot, mass, RHO, den):
    ''' Calculates the photo destruction rate depending on the amount of nucleobes:
        if there's enough nucleobase to cover the base of the pond then we multiply by A_p;
        if there isn't, we multiply by the cross-sectional area'''
    if mass / (RHO*den) < A_p:
        return -delta_t * muvdot * mass / (RHO*den)
    else:
        return -delta_t * muvdot * A_p

def decomposition(K_val, mass):
    ''' Calculates the decomposition/hydrolisis rate using the rate constant (k) adn the constant gamma'''
    return -delta_t * gamma * K_val * mass

def seepage(mass, mass_Met_or_IDP):
    ''' Calculates the seepage rate for each specias using a constant overall seepage rate'''
    return -delta_t*mass*m_seepage_rate/mass_Met_or_IDP

def get_mdot(NAME, spec_type):
    ''' Returns the m_dot value for a given species depending on the environment hardcoded in the paper (reducing vs oxidising)
        + treats hazes the same way'''
    if spec_type == 'red':
        return constants_and_rates[NAME]['m_dot_red']
    elif spec_type == 'ox':
        return constants_and_rates[NAME]['m_dot_ox']
    elif spec_type == 'min':
        return constants_and_rates[NAME]['m_dot_min'] # haze min
    elif spec_type == 'max':
        return constants_and_rates[NAME]['m_dot_max'] # haze max


# Solve ODE numerically
# Biomolecule evolution from meteorites, IDPs, and aqueous production from atmospheric precursors
for n in range(0,nt-1):
    
    for prec_mod in values: # cycle all precipitation models (IDN and COL)
        # first cycle through the Meteorite and IDP equations
        for impact in ['Met', 'IDP']:
            values[prec_mod]['L_'+impact][n+1] = water_decrease(impact, prec_mod, n)
            if (values[prec_mod]['L_'+impact][n+1] < 0):
                values[prec_mod]['L_'+impact][n+1] = 0
            if (values[prec_mod]['L_'+impact][n+1] >= (r_p - min_water)):
                values[prec_mod]['L_'+impact][n+1] = r_p - min_water
                values[prec_mod]['m_'+impact+'_A'][n+1] = values[prec_mod]['m_'+impact+'_A'][n] + mass_increase(m_dot_I, model = impact, wi = w_i, fs = f_s[impact]) + photo_destruction(constants_and_rates['adenine']['M_uv_dot'], values[prec_mod]['m_'+impact+'_A'][n], constants_and_rates['adenine']['rho'], constants_and_rates['adenine']['d'])
                if impact == 'Met':
                    pause_Met[prec_mod] += 1
            else:
                values[prec_mod]['m_'+impact+'_A'][n+1] = values[prec_mod]['m_'+impact+'_A'][n] + mass_increase(m_dot_I, model = impact, wi = w_i, fs = f_s[impact]) + nucleobase_outflow(n, pause_Met[prec_mod], constants_and_rates['adenine']['k'], model = impact) + decomposition(constants_and_rates['adenine']['k'], values[prec_mod]['m_'+impact+'_A'][n]) + seepage(values[prec_mod]['m_'+impact+'_A'][n], values[prec_mod]['m_'+impact][n])
            if values[prec_mod]['m_'+impact+'_A'][n+1] < 0:
                values[prec_mod]['m_'+impact+'_A'][n+1] = 0
            values[prec_mod]['m_'+impact][n+1] = water_mass(prec_mod, impact, n+1)
            # update concentration
            values[prec_mod]['C_'+impact][n+1] = values[prec_mod]['m_'+impact+'_A'][n+1] / values[prec_mod]['m_'+impact][n+1]    
        # then the aqueous equations
        for name in names:
            # we only work with masses in the pond, not from Met or IDP here (= aqueous production)
            if name[0] == 'm' and any([x in name for x in ['Met', 'IDP']]) == False:
                #print(name)
                spec_name = name.split('_')[1] # the format is m_spec_etc
                spec_mdot = get_mdot(spec_name, name.split('_')[-1]) # type is the last part of the srting
                
                if (values[prec_mod]['L_IDP'][n+1] >= (r_p - min_water)):
                    values[prec_mod]['L_IDP'][n+1] = r_p - min_water
                    values[prec_mod][name][n+1] = values[prec_mod][name][n] + mass_increase(spec_mdot) + photo_destruction(constants_and_rates[spec_name]['M_uv_dot'], values[prec_mod][name][n], constants_and_rates[spec_name]['rho'], constants_and_rates[spec_name]['d'])
                else:
                    values[prec_mod][name][n+1] = values[prec_mod][name][n] + mass_increase(spec_mdot) + decomposition(constants_and_rates[spec_name]['k'], values[prec_mod][name][n]) + seepage(values[prec_mod][name][n], values[prec_mod]['m_IDP'][n])
                if values[prec_mod][name][n+1] < 0:
                    values[prec_mod][name][n+1] = 0
                # update concentrations (biomolecule mass / water mass)
                conc_name = 'C' + name[1:]
                values[prec_mod][conc_name][n+1] = values[prec_mod][name][n+1] / values[prec_mod]['m_IDP'][n+1]
            


#Conversion from molar to mass mixing ratios
def molar2mass(x):
    return x * 1.e3 * constants_and_rates['adenine']['mu']

def mass2molar(x):
    return x / 1.e3 / constants_and_rates['adenine']['mu']

#Experimental yields
Adenine_lower = 0.005
Adenine_upper = 0.18
Guanine_lower = 6.7e-5
Guanine_upper = 0.2
Cytosine = 0.036
Uracil_lower = 1.7e-5
Uracil_upper = 0.018
Thymine = 0.012
Two_Amino_oxazole = 0.0011
Ribose = 3.6e-4
Formaldehyde = 0.036

#f, ax1 = plt.subplots(1, 1, figsize=(15,10))
f, (ax1, ax2) = plt.subplots(1, 2, figsize=(28,10))
        
p1 = ax1.fill_between(t, values['IDN']['C_adenine_red']*Adenine_lower*1e6/constants_and_rates['adenine']['mu'], values['IDN']['C_adenine_red']*Adenine_upper*1e6/constants_and_rates['adenine']['mu'], linestyle='-', color='#34b33a', lw=3.5, label=r'$\it{in}$ $\it{situ}$ Production - Early Hadean (reducing)', alpha=.60)
         
p2, = ax1.plot(t, values['IDN']['C_Met']*1e6/constants_and_rates['adenine']['mu'], linestyle='-', color='#9d620c', lw=3.5, label='Meteorites - Intermediate Env.')
         
p3, = ax1.plot(t, values['IDN']['C_IDP']*1e6/constants_and_rates['adenine']['mu'], linestyle='-', color='#9f53e6', lw=3.5, label='IDPs - Intermediate Env.')  
        
p4 = ax1.fill_between(t,values['IDN']['C_adenine_ox']*Adenine_lower*1e6/constants_and_rates['adenine']['mu'], values['IDN']['C_adenine_ox']*Adenine_upper*1e6/constants_and_rates['adenine']['mu'], linestyle='--', color='#095909', lw=3.5, label=r'$\it{in}$ $\it{situ}$ Production - Late Hadean (oxidizing)', alpha=.60)
                        
p5, = ax1.plot(t, values['COL']['C_Met']*1e6/constants_and_rates['adenine']['mu'], linestyle='--', color='#9d620c', lw=3.5, label='Meteorites - Wet Env.')      
               
#p6, = ax1.plot(t, C_IDN_65_Adenine_haze5*1e6/mu_Adenine, linestyle='-', color='#a30313', lw=3.5, label='Organic Hazes (5% methane)') 
  
p7 = ax1.fill_between(t, values['IDN']['C_adenine_haze05_min']*1e6/constants_and_rates['adenine']['mu'], values['IDN']['C_adenine_haze05_max']*1e6/constants_and_rates['adenine']['mu'], linestyle='-', hatch = '/',lw=3.5, color="#999999", label='Organic Hazes (0.5% methane)',alpha=0.4)                       
     
secax = ax1.secondary_yaxis('right', functions=(molar2mass, mass2molar))

secax.set_ylabel("Adenine Mass Fraction (ppb)",fontsize=18)

ax1.set_ylim(1e-11,50.)

ax1.legend([p1,p4,p2,p5,p3,p7], [r'aq. production - Early Hadean (reducing)', r'aq. production - Late Hadean (oxidizing)','Meteorites', 'Meteorites - no UV', 'IDPs',r'Organic Hazes (Early Hadean, 0.5% CH$_4$)'], loc=1,ncol=2, fontsize=16)
ax1.set_yscale('log')

  
ax2.plot(t, values['IDN']['C_uracil_haze05_max']*1e6/constants_and_rates['uracil']['mu'], linestyle='-', color='#AA3377', lw=3.5, label='Uracil')     
ax2.plot(t, values['IDN']['C_hypoxanthine_haze05_max']*1e6/constants_and_rates['hypoxanthine']['mu'], linestyle='--', color='#66CCEE', lw=3.5, label='Hypoxanthine')  
ax2.plot(t, values['IDN']['C_xanthine_haze05_max']*1e6/constants_and_rates['xanthine']['mu'], linestyle='-', color='#EE6677', lw=3.5, label='Xanthine*')   
ax2.plot(t, values['IDN']['C_thymine_haze05_max']*1e6/constants_and_rates['thymine']['mu'], linestyle='--', color='#BBBBBB', lw=3.5, label='Thymine')     
ax2.plot(t, values['IDN']['C_guanine_haze05_max']*1e6/constants_and_rates['guanine']['mu'], linestyle='-', color='#4477AA', lw=3.5, label='Guanine*')      
ax2.plot(t, values['IDN']['C_cytosine_haze05_max']*1e6/constants_and_rates['cytosine']['mu'], linestyle='--', color='#CCBB44', lw=3.5, label='Cytosine')    


ax2.set_yscale('log')
ax2.set_ylim(1.e-2,20e0)

ax2.text(-0.2,1.3e1,r"Maxima from organic hazes",fontweight='bold',fontsize=20)
ax2.text(-0.2,1.e1,r"Early Hadean, 0.5% CH$_4$",fontweight='bold',fontsize=20)

ax2.legend(loc=1,ncol=2,fontsize=16)

#for tick in ax1.yaxis.get_major_ticks()[::2]:
 #   tick.set_visible(False)
ax1.set_xlabel('Time (yr)', fontsize=18)
ax1.set_ylabel('Adenine Molar Concentration ($\mu$M)', fontsize=18)
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(18) 
for tick in ax1.yaxis.get_major_ticks():
    tick.label1.set_fontsize(18)
    
ax2.set_xlabel('Time (yr)', fontsize=18)
ax2.set_ylabel('Nucleobase Molar Concentration ($\mu$M)', fontsize=18)
for tick in ax2.xaxis.get_major_ticks():
    tick.label1.set_fontsize(18) 
for tick in ax2.yaxis.get_major_ticks():
    tick.label1.set_fontsize(18)

secax.tick_params(axis='both',labelsize=18)

ax1.text(-1.25, 1.7e1,'A)', fontsize=30, weight='bold', color='black')
ax2.text(-1.2, 1.52e1,'B)', fontsize=30, weight='bold', color='black')

def plot_save(fig, figtitle, folder, fonts = 35):
    ''' Saves the previously generated plot with a given title and within a given folder'''
    fig.suptitle(figtitle, fontsize = fonts)
    plt.tight_layout()
    if int(sys.argv[-1]) < 10:   
        plt.savefig(folder+'/Nucleobase_Comparison_0'+sys.argv[-1]+'.png',dpi=300)
    else:
        plt.savefig(folder+'/Nucleobase_Comparison_'+sys.argv[-1]+'.png',dpi=300)

if sys.argv[1] == 'temp':
    plot_save(f, u'T = {} \u00B0C'.format(sys.argv[2]), 'temperatures')
elif sys.argv[1] == 'hcn':
    plot_save(f, r'$mdot-HCN,in,red={:.3E}$   $mdot-HCN,in,ox={:.3E}$'.format(float(sys.argv[2]),float(sys.argv[3])), 'hcn_influx')
elif sys.argv[1] == 'h2cn':
    plot_save(f, r'$mdot-H2CN,in,red={:.3E}$   $mdot-H2CN,in,ox={:.3E}$'.format(float(sys.argv[2]),float(sys.argv[3])), 'h2cn_influx')
elif sys.argv[1] == 'uv':
    plot_save(f, r'$Fuv = {:.3f} W/m^2$'.format(float(sys.argv[2])), 'uv')
