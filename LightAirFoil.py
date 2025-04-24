# ======================================================================
#         Import modules
# ======================================================================
# rst imports (beg)
import os
import numpy as np
from mpi4py import MPI
from baseclasses import AeroProblem
from pygeo import DVConstraints, DVGeometryCST
from pyoptsparse import Optimization, OPT
from multipoint import multiPointSparse
from cmplxfoil import CMPLXFOIL, AnimateAirfoilOpt
import argparse 

# Use Python's built-in Argument parser to get commandline options
parser = argparse.ArgumentParser()
parser.add_argument("--output", type=str, default="output")
parser.add_argument("--opt", type=str, default="SLSQP", choices=["SLSQP", "SNOPT"])
parser.add_argument("--gridFile", type=str, default="n0012.cgns")
#parser.add_argument("--optOptions", type=ast.literal_eval, default={}, help="additional optimizer options to be added")
args = parser.parse_args()

# rst imports (end)

# ======================================================================
#         Specify parameters for optimization
# ======================================================================
# rst params (beg)
mycl = [0.5124, 0.445, 0.1965]  # lift coefficient constraint
mycpmin = [-0.75, -0.62, -0.44]
alpha = [2.18939, 1.57766, -0.6363] # if mycl == 0.0 else 1.0  # initial angle of attack (zero if the target cl is zero)
mach = [0.00404, 0.00428, 0.00687]   
Re = [394016, 417578, 669595]  # Reynolds number
T = 288.15  # 1976 US Standard Atmosphere temperature @ sea level (K)
nFlowCases = len(mach)
# rst params (end)

# ======================================================================
#         Create multipoint communication object
# ======================================================================
# rst procs (beg)
nGroup = 1
nProcPerGroup = MPI.COMM_WORLD.size
MP = multiPointSparse(MPI.COMM_WORLD)
MP.addProcessorSet("cruise", nMembers=nGroup, memberSizes=nProcPerGroup)
comm, setComm, setFlags, groupFlags, ptID = MP.createCommunicators()

# rst procs (end)

# ======================================================================
#         Create output directory
# ======================================================================
# rst dir (beg)
curDir = os.path.abspath(os.path.dirname(__file__))
outputDir = os.path.join(curDir, "output")

if not os.path.exists(outputDir):
    os.mkdir(outputDir)
# rst dir (end)

# ======================================================================
#         CFD solver set-up
# ======================================================================
# rst solver (beg)
aeroOptions = {
    "writeSolution": True,
    "writeSliceFile": True,
    "writeCoordinates": True,
    "plotAirfoil": True,
    "outputDirectory": outputDir,
    "nCrit": 3.87925,
}

# Create solver
CFDSolver = CMPLXFOIL(os.path.join(curDir, "naca0012.dat"), options=aeroOptions)
# rst solver (end)

# ======================================================================
#         Set up flow conditions with AeroProblem
# ======================================================================
# rst ap (beg)
aeroProblems = []
for i in range(nFlowCases):
    ap = AeroProblem(
        name="fc%d" % i,
        alpha=alpha[i], 
        #chordRef=0.12, 
        T=288.15, 
        mach=mach[i],
        reynolds = Re[i],
        reynoldsLength = 1,
        #rho = 1.225,
        #R=100, 
        #muSuthDim=1.22e-3, 
        #TSuthDim=288.15,
        evalFuncs=["cl", "cd", "cpmin"],
    )
    ap.addDV("alpha", value=alpha[i], lower=-5, upper=7, scale=1)
    aeroProblems.append(ap)
# Add angle of attack variable

# rst ap (end)

# ======================================================================
#         Geometric Design Variable Set-up
# ======================================================================
# rst geom (beg)
nCoeff = 5 # number of CST coefficients on each surface
DVGeo = DVGeometryCST(os.path.join(curDir, "naca0012.dat"), numCST=nCoeff)

DVGeo.addDV("upper_shape", dvType="upper", lowerBound=-0.1, upperBound=0.5)
DVGeo.addDV("lower_shape", dvType="lower", lowerBound=-0.5, upperBound=0.1)

# Add DVGeo object to CFD solver
CFDSolver.setDVGeo(DVGeo)
# rst geom (end)

# ======================================================================
#         DVConstraint Setup
# ======================================================================
# rst cons (beg)
DVCon = DVConstraints()
DVCon.setDVGeo(DVGeo)
DVCon.setSurface(CFDSolver.getTriangulatedMeshSurface())

# Thickness, volume, and leading edge radius constraints
le = 0.0001
wingtipSpacing = 0.1
leList = [[le, 0, wingtipSpacing], [le, 0, 1.0 - wingtipSpacing]]
teList = [[1.0 - le, 0, wingtipSpacing], [1.0 - le, 0, 1.0 - wingtipSpacing]]
DVCon.addVolumeConstraint(leList, teList, 2, 100, lower=0.85, upper=3, scaled=True)
DVCon.addThicknessConstraints2D(leList, teList, 2, 100, lower=0.25, upper = 3, scaled=True)
#DVCon.addSurfaceAreaConstraint(lower=0.85)
le = 0.01
leList = [[le, 0, wingtipSpacing], [le, 0, 1.0 - wingtipSpacing]]
DVCon.addLERadiusConstraints(leList, 2, axis=[0, 1, 0], chordDir=[-1, 0, 0], lower=0.85, scaled=True)

fileName = os.path.join(outputDir, "constraints.dat")
DVCon.writeTecplot(fileName)
# rst cons (end)


# ======================================================================
#         Functions:
# ======================================================================
# rst funcs (beg)
def cruiseFuncs(x):
    print(x)
    # Set design vars
    DVGeo.setDesignVars(x)
    # Run CFD
    for i in range(nFlowCases):
        if i % nGroup == ptID:
            aeroProblems[i].setDesignVars(x)
            CFDSolver(aeroProblems[i])
    # Evaluate functions
    funcs = {}
    DVCon.evalFunctions(funcs)
    for i in range(nFlowCases):
        if i % nGroup == ptID:
            CFDSolver.evalFunctions(aeroProblems[i], funcs)
            CFDSolver.checkSolutionFailure(aeroProblems[i], funcs)
    if MPI.COMM_WORLD.rank == 0:
        print("functions:")
        for key, val in funcs.items():
            if key == "DVCon1_thickness_constraints_0":
                continue
            print(f"    {key}: {val}")
    return funcs


def cruiseFuncsSens(x, funcs):
    funcsSens = {}
    DVCon.evalFunctionsSens(funcsSens)
    for i in range(nFlowCases):
        if i % nGroup == ptID:
            CFDSolver.evalFunctionsSens(aeroProblems[i], funcsSens)
            CFDSolver.checkAdjointFailure(aeroProblems[i], funcsSens)
    print("function sensitivities:")
    evalFunc = ["fc0_cd", "fc1_cd", "fc0_cl","fc1_cl", "fail"]
    for var in evalFunc:
        print(f"    {var}: {funcsSens[var]}")
    return funcsSens


def objCon(funcs, printOK):
    # Assemble the objective and any additional constraints:
    funcs["obj"] = 0.0
    for i in range(nFlowCases):
        ap = aeroProblems[i]
        funcs["obj"] += funcs[ap["cd"]] /nFlowCases
        funcs["cl_con_" + ap.name] = funcs[ap["cl"]] - mycl[i]
        funcs["cpmin_con_" + ap.name] = funcs[ap["cpmin"]]
    if printOK:
        print("funcs in obj:", funcs)
    return funcs


# rst funcs (end)

# ======================================================================
#         Optimization Problem Set-up
# ======================================================================
# rst optprob (beg)
# Create optimization problem
optProb = Optimization("opt", MP.obj)

# Add objective
optProb.addObj("obj", scale=1e6)

# Add variables from the AeroProblem
for ap in aeroProblems:
    ap.addVariablesPyOpt(optProb)
# Add DVGeo variables
DVGeo.addVariablesPyOpt(optProb)

# Add constraints
DVCon.addConstraintsPyOpt(optProb)

# Add cl constraint
for ap in aeroProblems:
    optProb.addCon(f"cl_con_{ap.name}", lower=0.0, upper=0.0, scale=10, wrt=[f"alpha_{ap.name}", "upper_shape", "lower_shape"])

i = 0
for ap in aeroProblems:
    optProb.addCon("cpmin_con_" + ap.name, lower = mycpmin[i], upper = 0.0, scale = 1.0, wrt=[f"alpha_{ap.name}", "upper_shape", "lower_shape"])
    i = i+1


# Enforce first upper and lower CST coefficients to add to zero
# to maintain continuity at the leading edge
jac = np.zeros((1, nCoeff), dtype=float)
jac[0, 0] = 1.0
optProb.addCon(
    "first_cst_coeff_match",
    lower=0.0,
    upper=0.0,
    linear=True,
    wrt=["upper_shape", "lower_shape"],
    jac={"upper_shape": jac, "lower_shape": jac},
)

# The MP object needs the 'obj' and 'sens' function for each proc set,
# the optimization problem and what the objcon function is:
MP.setProcSetObjFunc("cruise", cruiseFuncs)
MP.setProcSetSensFunc("cruise", cruiseFuncsSens)
MP.setObjCon(objCon)
MP.setOptProb(optProb)
optProb.printSparsity()
optProb.getDVConIndex()
# rst optprob (end)

# rst opt (beg)
# Run optimization



#optOptions = {"IFILE": os.path.join(outputDir, "SLSQP.out")}
optOptions = {"IFILE": os.path.join(outputDir, "SLSQP.out"), "ACC": 1e-6}
opt = OPT("SLSQP", options=optOptions)
sol = opt(optProb, MP.sens, storeHistory=os.path.join(outputDir, "opt.hst"))
if MPI.COMM_WORLD.rank == 0:
    print(sol)
# rst opt (end)

# ======================================================================
#         Postprocessing
# ======================================================================
# rst postprocessing (beg)
# Save the final figure
CFDSolver.airfoilAxs[1].legend(["Original", "Optimized"], labelcolor="linecolor")
CFDSolver.airfoilFig.savefig(os.path.join(outputDir, "OptFoil.pdf"))

# Animate the optimization
AnimateAirfoilOpt(outputDir, "fc").animate(
    outputFileName=os.path.join(outputDir, "OptFoil"), fps=10, dpi=300, extra_args=["-vcodec", "libx264"]
)
# rst postprocessing (end)
