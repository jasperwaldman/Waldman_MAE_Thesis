### Calculate Parasite Drag for range of Re for a VSP Geometry

import numpy as np
import openvsp as vsp
import pandas as pd


def CalcParasiteDrag(fname, reynoldsnumber, mach, lampercent):
    # Import Geometry
    vsp.ReadVSPFile(fname)
    geomid = vsp.FindGeoms()
    SRef = vsp.GetParmVal(geomid[0], "TotalArea", "WingGeom")
    print("Loading Geometry from File:", fname)
    print(SRef)
    # Set Parasite Drag Settings
    PD_analysis = "ParasiteDrag"
    vsp.SetAnalysisInputDefaults(PD_analysis) # Set the defaults for the analysis type
    vsp.SetIntAnalysisInput(PD_analysis, "LengthUnit", (2,))  # meters
    vsp.SetDoubleAnalysisInput(PD_analysis, "Sref", (SRef,))
    vsp.SetIntAnalysisInput(PD_analysis, "FreestreamPropChoice", (5,))
    vsp.SetDoubleAnalysisInput(PD_analysis, "SpecificHeatRatio", (1.4,))  # For air, gamma = Cp/Cv = 1.4
    vsp.SetDoubleAnalysisInput(PD_analysis, "Re_L", (reynoldsnumber, ))
    vsp.SetDoubleAnalysisInput(PD_analysis, "Mach", (mach, ))  # Mach number for the flow
    vsp.SetParmVal(geomid[0], "PercLam", "ParasiteDragProps", lampercent)
    #vsp.SetDoubleAnalysisInput(PD_analysis, "PercLam", (lampercent, ))
    vsp.PrintAnalysisInputs(PD_analysis)

    vsp.Update()

    # Run analysis and return results
    print("Running Parasite Drag Analysis...")
    ridpd = vsp.ExecAnalysis(PD_analysis)

    dat = vsp.GetDoubleResults( ridpd, "Total_CD_Total", 0 )
    print(dat)
    vsp.ClearVSPModel()
    return dat[0]

def ExtractLamPercent(re, flap):
    df = pd.read_csv("re_flap_xtr_summary.csv")
    # Define the input values (adjust these as needed)
    # Find the matching row
    match = df[(df["reynolds"] == re) & (df["flap_angle"] == flap)]

    if not match.empty:
        xtr_value = match.iloc[0]["xtr_top"]
        print(f"xtr_top for Re = {re}, flap = {flap}°: {xtr_value}")
        return xtr_value*100
    else:
        print(f"No match found for Re = {re}, flap = {flap}")
        return 0

if __name__ == "__main__":
    
    reynoldsnumbers = range(200000, 2000000, 100000)
    flap_angles = range(0, 1)
    results = []
    ref_length = 0.083
    mach = 0.01
    
    fname = f"107.vsp3"
    for flap in flap_angles:
        for re in reynoldsnumbers:
                #lampercent = ExtractLamPercent(re, flap)
                re_l = re/ref_length
                cd0 = CalcParasiteDrag(fname, re_l, mach, 0)  # Calculate parasite drag for the given file, Re, and Mach number
                results.append([re, cd0]) # Append the result to the list
    
    print(results)
    np.savetxt("mainfoilresults.txt", results, header="Re, Cd0", 
              fmt=["%.2e", "%.6f"], delimiter=", ")
  