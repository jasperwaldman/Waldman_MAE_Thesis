### CREATES A VSP GEOMETRY FROM SYN107 RESULTS. 
### AIRFOIL PREPROCESSING REQUIRED IN XFOIL, NEEDS INTEGRATION, IS NOT DONE, BANDAID FIX FOR THESIS APRIL 16 2025
### AUTHOR: JASPER WALDMAN
### DATE: APRIL 16 2025


import numpy as np
from pathlib import Path
import openvsp as vsp
import re
import matplotlib.pyplot as plt
def parse_foil_file(filepath):
    """
    Parses a .foil file with multiple sections, extracts and normalizes airfoil data.

    Returns a list of tuples with (chord, deltaX, deltaY) for each section.
    Writes a normalized coordinate file per section.
    """
    filepath = Path(filepath)  # Ensure it's a Path object
    with open(filepath, 'r') as file:
        lines = file.readlines()

    results = []
    section_count = 0
    i = 0
    while i < len(lines):
        if "(upper surface)" in lines[i]:
            section_count += 1
            index = section_count - 1

            # Extract xl, yl, chord from the appropriate line before the surface points
            info_line = lines[i - 5].split()
            xl, yl, chord = float(info_line[1]), float(info_line[2]), float(info_line[3])

            # Extract number of points from the nl/nu line
            nl_line = lines[i - 3].split()
            n_upper = int(float(nl_line[2]))
            n_lower = int(float(nl_line[1]))

            # Read upper surface points
            upper_lines = lines[i + 1:i + 1 + n_upper]
            upper_surface = np.array([[float(val) for val in line.split()[:2]] for line in upper_lines])

            # Read lower surface points (skipping label line)
            lower_start = i + 2 + n_upper
            lower_lines = lines[lower_start:lower_start + n_lower]
            lower_surface = np.array([[float(val) for val in line.split()[:2]] for line in lower_lines])

            # Process the points
            upper_sorted = upper_surface[np.argsort(-upper_surface[:, 0])]
            upper_sorted = upper_sorted[:-1]  # remove (0,0)
            coords = np.vstack((upper_sorted, lower_surface))
            coords_norm = coords / chord

            # Save to file
            out_path = filepath.parent / f"normalized_section_{index}.dat"
            #np.savetxt(out_path, coords_norm, fmt="%.6f", header="x_norm y_norm", comments='')

            # Normalize xl and yl
            deltaX = xl/chord
            deltaY = yl/chord

            results.append((chord, deltaX, deltaY))

            # Move index past this section
            i = lower_start + n_lower
        else:
            i += 1

    return results

def makefoil(results, flap):
    mainfoil = vsp.AddGeom("WING")
    vsp.SetDriverGroup(mainfoil, 1, vsp.SPAN_WSECT_DRIVER, vsp.ROOTC_WSECT_DRIVER, vsp.TIPC_WSECT_DRIVER)
    vsp.SetParmValUpdate(mainfoil, "Sweep", "XSec_1", 0)
    vsp.SetParmValUpdate(mainfoil, "TotalArea", "WingGeom", 0.085)
    vsp.SetParmValUpdate(mainfoil, "Span", "XSec_1", 0.01625)
    vsp.SetParmValUpdate(mainfoil, "Root_Chord", "XSec_1", results[0][0])
    vsp.SetParmValUpdate(mainfoil, "Tip_Chord", "XSec_1", results[1][0])
    vsp.SetParmValUpdate(mainfoil, "DeltaX", "XSecCurve_0", results[0][1])
    vsp.SetParmValUpdate(mainfoil, "DeltaY", "XSecCurve_0", results[0][2])
    vsp.SetParmValUpdate(mainfoil, "DeltaX", "XSecCurve_1", results[1][1])
    vsp.SetParmValUpdate(mainfoil, "DeltaY", "XSecCurve_1", results[1][2])

    
    # Set up first section (first two airfoil sections)
    xsecsurf = vsp.GetXSecSurf(mainfoil, 0)
    vsp.ChangeXSecShape(xsecsurf, 0, vsp.XS_FILE_AIRFOIL)
    vsp.ChangeXSecShape(xsecsurf, 1, vsp.XS_FILE_AIRFOIL)
    xsec0 = vsp.GetXSec(xsecsurf, 0)
    xsec1 = vsp.GetXSec(xsecsurf, 1)
    vsp.ReadFileAirfoil(xsec0, f"smoothed_section_0_flap_{flap}.dat")
    vsp.ReadFileAirfoil(xsec1, f"smoothed_section_1_flap_{flap}.dat")



    #Insert XSecs after first section and set airfoil files
    for i in range(2, 33):
        vsp.InsertXSec(mainfoil, i-1, vsp.XS_FILE_AIRFOIL)
        xsecid = vsp.GetXSec(xsecsurf, i)
        vsp.ReadFileAirfoil(xsecid, f"smoothed_section_{i}_flap_{flap}.dat")
        vsp.SetParmValUpdate(mainfoil, "Span", f"XSec_{i}", 0.01625)
        vsp.SetParmValUpdate(mainfoil, "SectTess_U", f"XSec_{i}", 2)
        vsp.SetParmValUpdate(mainfoil, "Tip_Chord", f"XSec_{i}", results[i][0])
        vsp.SetParmValUpdate(mainfoil, "Root_Chord", f"XSec_{i}", results[i-1][0])
        vsp.SetParmValUpdate(mainfoil,  "DeltaX", f"XSecCurve_{i}", results[i][1])
        vsp.SetParmValUpdate(mainfoil,  "DeltaY", f"XSecCurve_{i}", results[i][2])
        #vsp.SetParmValUpdate(mainfoil,  "Theta", f"XSecCurve_{i}", 180)


        
    vsp.WriteVSPFile(f"107_flap_{flap}.vsp3", vsp.SET_ALL)

flap = -10
results = parse_foil_file("jasper.foil")
print(results[0], results[1])
makefoil(results, flap)
# Run the function on the file
#all_sections_data = parse_foil_file("./jasper.foil")


