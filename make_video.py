''' A python script to run the pond model with a varying paramter and then make a short video out of it.'''

import subprocess
import glob
from PIL import Image
import numpy as np
import sys

HCN_mass_influx_ox = 3.42599458e-12
HCN_mass_influx_red = 1.25776186e-07
H2CO_mass_influx_red = 6.01494508e-12
H2CO_mass_influx_ox = 1.15285175e-09#4.11794768e-09

if sys.argv[1] == 'temp':
    temperatures = np.linspace(20, 99, 15)
    subprocess.run(['rm', '-rf', 'temperatures'])
    subprocess.run(['mkdir', 'temperatures'])
    for i,temp in enumerate(temperatures):
        subprocess.run(['python', 'Wet_Dry_Cycling_Pond_Model_interactive.py', 'temp', str(temp), str(i)], check = True, text = True)
elif sys.argv[1] == 'hcn':
    hcn_influx_red = np.linspace(HCN_mass_influx_red/10,HCN_mass_influx_red*10, 15) # varying the HCN influx by an order of magnitude
    hcn_influx_ox = np.linspace(HCN_mass_influx_ox/10,HCN_mass_influx_ox*10, 15) # for both the reducing and oxidising atmosphere
    subprocess.run(['rm', '-rf', 'hcn_influx'])
    subprocess.run(['mkdir', 'hcn_influx'])
    for i in range(len(hcn_influx_red)):
        subprocess.run(['python', 'Wet_Dry_Cycling_Pond_Model_interactive.py', 'hcn', str(hcn_influx_red[i]), str(hcn_influx_ox[i]), str(i)], check = True, text = True)
elif sys.argv[1] == 'h2co':
    h2cn_influx_red = np.linspace(H2CO_mass_influx_red/10,H2CO_mass_influx_red*10, 15) # varying the formaldehyde influx by an order of magnitude
    h2cn_influx_ox = np.linspace(H2CO_mass_influx_ox/10,H2CO_mass_influx_ox*10, 15) # for both the reducing and oxidising atmosphere
    subprocess.run(['rm', '-rf', 'h2co_influx'])
    subprocess.run(['mkdir', 'h2co_influx'])
    for i in range(len(h2cn_influx_red)):
        subprocess.run(['python', 'Wet_Dry_Cycling_Pond_Model_interactive.py', 'h2co', str(h2cn_influx_red[i]), str(h2cn_influx_ox[i]), str(i)], check = True, text = True)
elif sys.argv[1] == 'uv':
    flux = np.linspace(0, 2, 15) # og was 0.4
    subprocess.run(['rm', '-rf', 'uv'])
    subprocess.run(['mkdir', 'uv'])
    for i,f in enumerate(flux):
        subprocess.run(['python', 'Wet_Dry_Cycling_Pond_Model_interactive.py', 'uv', str(f), str(i)], check = True, text = True)


def make_gif(frame_folder, gif_name):
    frames = [Image.open(image) for image in sorted(glob.glob(f"{frame_folder}/*.png"))]
    frame_one = frames[0]
    frame_one.save(gif_name, format="GIF", append_images=frames,
               save_all=True, duration=1000, loop=0)
    
if sys.argv[1] == 'temp':
    make_gif('temperatures', "temperature_dependence.gif")
elif sys.argv[1] == 'hcn':
    make_gif('hcn_influx', 'hcn_influx_dependencies.gif')
elif sys.argv[1] == 'h2cn':
    make_gif('h2cn_influx', 'h2cn_influx_dependencies.gif')
elif sys.argv[1] == 'uv':
    make_gif('uv', 'uv_dependencies.gif')