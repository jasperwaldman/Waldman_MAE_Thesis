import os
import subprocess
import shutil
import numpy as np
from concurrent.futures import ThreadPoolExecutor
import time
import itertools
import pandas as pd

FLOW_TEMPLATE = """
##########################################################
#####       Input file for computations               ###
##########################################################

+
++ Name ++
+++
simulation_foil

+
++ No. bodies ++ # number of bodies
+++
1

+
++ Body 1: Mesh++ # mesh file
+++
FLAP_MESH_FILE

+
++ Body 1: Mesh: Translations++ # translation to apply (default 0 0 0)
+++
0 0 0

+
++ Body 1: Mesh: Rotations++ # rotation to apply (default 0 0 0)
+++
HEEL_VALUE PITCH_VALUE LEEWAY_VALUE

+
++ Body 1: Mesh: Center of rotations ++ # center of rotation (default 0 0 0)
+++
-0.376 0.0 -1.060

+
++ Body 1: Lifting ++ # should be 1 for foils
+++
1

+
++ Free stream velocities ++ # velocity coordinates in x,y,z-directions (m/s)
+++
-6 0.0 0.0
    
+
++ Fluid density ++ # kg/m3
+++
1027

+
++ Number of time steps ++ 
+++
200

+
++ Time step value ++ # in seconds
+++
0.004

+
++ Wake: limitation method ?++ # 0=no, 1= limit the number of panels in the wake
+++
1

+
++ Wake: limitation: value  ++ # 0=no, 1= limit the number of panels in the wake
+++
100

+
++ Free Surface ? ++ # 0= no FS (default), >0= FS considered (optional) should be 5 for Neumann-Kelvin condition !
+++
0

+
++ Free Surface: initial position ++ # z coordinate of the initial FS (optional) (default=0.0)
+++
FREE_SURFACE_VALUE

+
++ Free Surface: boundaries ++ # size of the FS domain min(x), max(x), min(y), max(y) 
+++
-10 2 -4 +4

+
++ Reference point ++
+++ 
-0.376 0 -1.060

+
++ Free Surface: no. panels ++ # no. panels in x-y directions 
+++
60 40
+
++ Probes interval ++ # 0=no probes (default), or provide the number of time step between probes (optional)
+++    
10

+
++ Steady ? ++ # remove time dependance # should be used for steady computations
+++
1

+
++ Check convergence ? ++  # computation stop if convergence is obtained within n iterations
+++
70
+
++ Convergence tolerance ? ++
+++
0.001

+
++ No. cores ++ # number of cores to use, defaults is 2
+++    
8

+
++ Results format ++ # 0=tecplot (default), 1=vtk for paraview
+++
1
"""

def modify_flow_dat(heel, pitch, leeway, flap, free_surface, case_dir):
    """Modifies flow.dat for a given set of conditions and stores it in the case directory."""
    data = FLOW_TEMPLATE.replace('HEEL_VALUE', str(heel)) # - heel in puffin is positive heel in FS. Flipped 
    data = data.replace('PITCH_VALUE', str(pitch))
    data = data.replace('LEEWAY_VALUE', str(leeway)) # POSITIVE LEEWAY IN VPP IS TO PORT. 
    data = data.replace('FREE_SURFACE_VALUE', str(free_surface))
    mesh_filename = f'mainfoil_flap_{int(flap)}.dat'
    data = data.replace('FLAP_MESH_FILE', mesh_filename)
    
    os.makedirs(case_dir, exist_ok=True)
    with open(os.path.join(case_dir, 'flow.dat'), 'w') as f:
        f.write(data)
    if os.path.exists(mesh_filename):
        shutil.copy(mesh_filename, case_dir)

def run_puffin(case_dir):
    """Runs PUFFin in the case directory and handles errors."""
    try:
        print(f"Running {case_dir}")
        result = subprocess.run(['puffin.exe', rf"flow= .\{case_dir}\flow.dat"], capture_output=False)
        return result.returncode == 0
    except subprocess.TimeoutExpired:
        return False

def extract_last_line(filename):
    """Extracts the last line of a file and returns only the required columns."""
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
            last_line = lines[-1].strip().split()
            
            if 'forces_body' in filename:
                print("Extracting Forces")
                return list(map(float, last_line[1:4]))  # Extract Fx, Fy, Fz
            elif 'moments_body' in filename:
                print("Extracting Moments")
                return list(map(float, last_line[1:4]))  # Extract Mx, My, Mz
    except Exception:
        return None


def process_case(case_index, heel, pitch, leeway, flap, free_surface, failed_cases, ref_area, ref_chord, dyn_pressure):
    case_dir = f'case_{case_index+1}'
    modify_flow_dat(heel, pitch, leeway, flap, free_surface, case_dir)
    
    if run_puffin(case_dir):
        forces = extract_last_line(os.path.join(case_dir, 'forces_body_0001_PUFFIn.dat'))
        moments = extract_last_line(os.path.join(case_dir, 'moments_body_0001_PUFFIn.dat'))
        if forces and moments:
            CFx, CFy, CFz = [F / (dyn_pressure * ref_area) for F in forces]
            CMx, CMy, CMz = [M / (dyn_pressure * ref_area * ref_chord) for M in moments]
            return (heel, pitch, leeway, flap, free_surface,  CFx, CFy, CFz, CMx, CMy, CMz)
    print(f"Case failed: Heel={heel}, Pitch={pitch}, Leeway={leeway}, Flap={flap}, Free Surface={free_surface}")
    failed_cases.append((heel, pitch, leeway, flap, free_surface))
    return None



def main():
    start_time = time.time()
    #Define ranges for operating conditions
    heel_angles = [0] # Convention: >0 Windward, <0 Leeward
    pitch_angles = [1,0.6,0.2,-0.4,-0.8, -1.2, -1.6] # Convention: >0 Bow Down, <0 Bow Up
    leeway_angles = [0] # Convention: >0 Bow to Windward, <0 Bow to leeward
    flap_angles = [0] # Positive Down
    free_surface_heights = [-0.79] # Equivalent to Sink. Max = 0 
    
    # User-provided reference values
    ref_area = 0.0815  # Reference Area
    ref_chord = 0.08986  # Reference Chord
    rho = 1027 # kg/m^3
    u = 6 # m/s
    dyn_pressure = 0.5*rho*u**2  # Example dynamic pressure


    # Generate all combinations of operating conditions
    operating_conditions = list(itertools.product(heel_angles, pitch_angles, leeway_angles, flap_angles, free_surface_heights))
    #df = pd.read_csv("256_sobol_samples.csv")
    #operating_conditions = df.values.tolist()
    print(operating_conditions)
    print(f"Generating Hydro Data for {len(operating_conditions)} Cases!")
   

    
    results = []
    failed_cases=[]
    with ThreadPoolExecutor(max_workers=4
                            ) as executor:
        futures = [executor.submit(process_case, i, *conditions, failed_cases, ref_area, ref_chord, dyn_pressure) for i, conditions in enumerate(operating_conditions)]
        for future in futures:
            result = future.result()
            if result:
                results.append(result)
    
    for idx, coeff_name in enumerate(['CFx', 'CFy', 'CFz', 'CMx', 'CMy', 'CMz'], start=5):
        with open(f'{coeff_name}_results.txt', 'w') as f:
            f.write('Heel Pitch Leeway Flap FreeSurface {}\n'.format(coeff_name))
            for res in results:
                f.write(f'{res[0]} {res[1]} {res[2]} {res[3]} {res[4]} {res[idx]}\n')

    np.savetxt('All_results.txt', results, fmt='%s', header='Heel Pitch Leeway Flap FreeSurface CFx CFy CFz CMx CMy CMz')

    with open('cases.txt', 'w') as f:
        f.write('Case Heel Pitch Leeway Flap FreeSurface\n')
        for i, (heel, pitch, leeway, flap, free_surface) in enumerate(operating_conditions):
            f.write(f'{i+1} {heel} {pitch} {leeway} {flap} {free_surface}\n')

    if failed_cases:
        print("\nFailed Cases:")
        for case in failed_cases:
            print(f"Heel={case[0]}, Pitch={case[1]}, Leeway={case[2]}, Flap={case[3]}, Free Surface={case[4]}")
    
    print("Completed all other cases.")
    end_time = time.time()
    print(f"Execution time: {end_time - start_time:.2f} seconds")

if __name__ == "__main__":
    main()
