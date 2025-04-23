# Handles parametric geometry generation via OpenVSP API and mesh conversion for PUFFIn BEM Code
# Author: Jasper Waldman 

import openvsp as vsp
import numpy as np
def VSP2PUFFin(input_file, output_file, name):
    # Open the input file
    try:
        with open(input_file, 'r') as file:
            raw_text = file.readlines()
    except Exception as e:
        raise RuntimeError(f"Could not open input file: {input_file}. Error: {e}")
    
    points = []
    cross_sections = []
    current_component = 0
    current_section = 0

    for line in raw_text:
        line = line.strip()
        
        # Skip empty lines or lines starting with comments
        if not line or line.startswith('HERMITE') or line.startswith('NUMBER') or line.startswith('GROUP') or line.startswith('TYPE'):
            continue

        # Parse "WingGeom" component and section info
        if 'WingGeom' in line:
            current_component += 1
            current_section = 0
            continue

        elif 'CROSS SECTIONS' in line:
            # Read number of cross-sections
            num_u = int(line.split('=')[-1].strip())
            continue
        elif 'PTS/CROSS SECTION' in line:
            # Read number of points per cross-section
            num_w = int(line.split('=')[-1].strip())
            continue
        
        # Parse the coordinates in the subsequent lines
        if line:
            coords = np.fromstring(line, sep=' ')
            if len(coords) == 3:
                points.append(coords)
            if len(points) == num_w:
                current_section += 1
                if len(cross_sections) < current_section:
                    cross_sections.append([])
                cross_sections[current_section - 1].append(np.array(points))  # Add current cross-section to the list
                points = []  # Reset for the next section

    num_nodes = (2 * num_u - 1) * (num_w)

    # Now, write the parsed data to the output file
    try:
        with open(output_file, 'w') as f_out:
            f_out.write(f"{name}\n")
            f_out.write(f"{2 * num_u - 2} {num_w - 1} {num_nodes}\n")
            f_out.write('nodes\n')

            node_count = 1
            # Write node data for the first component
            for i in range(num_u):
                station_points = cross_sections[-(i + 1)][0]  # Extract station points from CrossSections
                for j in range(num_w):
                    f_out.write(f"{node_count}, {station_points[j, 0]:.16f} {station_points[j, 1]:.16f} {station_points[j, 2]:.16f}\n")
                    node_count += 1

            # Write node data for the second component
            for i in range(1, num_u):
                station_points = cross_sections[i][1]  # Extract station points from CrossSections
                for j in range(num_w):
                    f_out.write(f"{node_count}, {station_points[j, 0]:.16f} {station_points[j, 1]:.16f} {station_points[j, 2]:.16f}\n")
                    node_count += 1

            # Connectivity
            f_out.write('connectivities\n')
            num_panels = (2 * num_u - 2) * (num_w - 1)
            panel_count = 1
            i = 0
            while panel_count <= num_panels:
                if (i % num_w) != (num_w - 1):  # Skip the last node
                    f_out.write(f"{panel_count}, {i + 1}, {i + num_w + 1}, {i + num_w + 2}, {i + 2}\n")
                    panel_count += 1
                i += 1

            f_out.write('0 0 0 0 0')  # Ending code

    except Exception as e:
        raise RuntimeError(f"Could not open output file: {output_file}. Error: {e}")

def CreateMainFoil(sectionarea, sectionspan, taper, dihedral, sweep, sweep_location, twist, airfoilfile, flap):
    mainfoil = vsp.AddGeom("WING")
    vsp.SetDriverGroup(mainfoil, 1, vsp.SPAN_WSECT_DRIVER, vsp.AREA_WSECT_DRIVER, vsp.TAPER_WSECT_DRIVER)
    # Wing Planform
    vsp.SetParmValUpdate(mainfoil, "X_Rel_Location", "XForm", -0.376)
    vsp.SetParmValUpdate(mainfoil, "Z_Rel_Location", "XForm", -1.060)
    vsp.SetParmValUpdate(mainfoil, "Z_Rel_Rotation", "XForm", 180)
    vsp.SetParmValUpdate(mainfoil, "Span", "XSec_1", sectionspan)
    vsp.SetParmValUpdate(mainfoil, "Area", "XSec_1", sectionarea)
    vsp.SetParmValUpdate(mainfoil, "Taper", "XSec_1", taper)
    vsp.SetParmValUpdate(mainfoil, "Dihedral", "XSec_1", dihedral)
    vsp.SetParmValUpdate(mainfoil, "Sweep", "XSec_1", sweep)
    vsp.SetParmValUpdate(mainfoil, "Sec_Sweep_Location", "XSec_1", 0.25)
    vsp.SetParmValUpdate(mainfoil, "Sweep_Location", "XSec_1", sweep_location)
    vsp.SetParmValUpdate(mainfoil, "Sec_Sweep_Location", "XSec_1", 0)
    vsp.SetParmValUpdate(mainfoil, "Twist", "XSec_1", twist)
    
    # Root and End Caps
    vsp.SetParmValUpdate(mainfoil, "CapUMinOption", "EndCap", 0) 
    vsp.SetParmValUpdate(mainfoil, "CapUMaxOption", "EndCap", 0)

    # Tesselation
    vsp.SetParmValUpdate(mainfoil, "OutCluster", "XSec_1", 0.1)
    vsp.SetParmValUpdate(mainfoil, "SectTess_U", "XSec_1", 29)
    vsp.SetParmValUpdate(mainfoil, "Tess_W", "Shape", 121)

    # Airfoil Root and Tip
    xsec_surf = vsp.GetXSecSurf(mainfoil, 0)
    vsp.ChangeXSecShape(xsec_surf, 0, vsp.XS_FILE_AIRFOIL)
    xsec = vsp.GetXSec(xsec_surf, 0)
    vsp.ReadFileAirfoil(xsec, airfoilfile) 
    vsp.ChangeXSecShape(xsec_surf, 1, vsp.XS_FILE_AIRFOIL)
    xsec = vsp.GetXSec(xsec_surf, 1)
    vsp.ReadFileAirfoil(xsec, airfoilfile)

    vsp.Update()

    # Blending
    vsp.SetParmValUpdate(mainfoil, "OutLEMode", "XSec_0", 1)
    vsp.SetParmValUpdate(mainfoil, "OutLESweep", "XSec_0", 0)
    vsp.SetParmValUpdate(mainfoil, "OutTEMode", "XSec_0", 1)
    vsp.SetParmValUpdate(mainfoil, "OutTESweep", "XSec_0", 0)

    # Export
    vsp.ExportFile(f"mainfoil_flap_{flap}.hrm", vsp.SET_ALL, vsp.EXPORT_XSEC)
    vsp.WriteVSPFile(f"mainfoil_flap_{flap}.vsp3", vsp.SET_ALL) #don't need this for batch use. only need the mesh file. 
    vsp.VSPRenew()

if __name__ == "__main__":
    totalarea = 0.0815
    totalspan = 1.04
    # RootChord = 0.120
    # TipChord = 0.05
    sectionspan = totalspan/2
    sectionarea = totalarea/2
    Dihedral = 0
    Sweep = 0
    SweepLocation = 1
    Twist = 0
    taper = 0.5
    flap = -10
    #CreateMainFoil(sectionarea, sectionspan, taper, Dihedral, Sweep, SweepLocation, Twist, f'LightAir_V8_ClosedTE_flap_{flap}.dat', flap)
    VSP2PUFFin(f'opt_chord_flap_{flap}.hrm', f'opt_chord_flap_{flap}.dat', 'mainfoil')