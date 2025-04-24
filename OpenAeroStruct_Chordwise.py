"""
Use Mphys/OpenVSP interface to perform inviscid drag minimization of
an initially rectangular wing with respect to the chord distribution,
subject to a lift and reference area constraint. Similar to the twist optimization,
the expected result from lifting line theory should produce an elliptical lift distrbution.
Check output directory for Tecplot solution files.
"""

import os
import openmdao.api as om
import numpy as np

from mphys.core import Multipoint, MPhysVariables
from mphys.scenarios import ScenarioAerodynamic
from pygeo.mphys import OM_DVGEOCOMP
from pygeo import DVGeometryVSP
from openaerostruct.geometry.utils import generate_vsp_surfaces
from openaerostruct.mphys import AeroBuilder

class Top(Multipoint):
    def setup(self):
        # VSP model
        vsp_file = os.path.join(os.path.dirname(__file__), "mainfoil_18_sect.vsp3")

        # Generate half-body mesh of rectangular wing
        surfaces = generate_vsp_surfaces(vsp_file, symmetry=True, include=["WingGeom"])

        surf_options = {
            "type": "aero",
            "S_ref_type": "wetted",  # how we compute the wing area,
            # can be 'wetted' or 'projected'
            "twist_cp": np.zeros(3),  # Define twist using 3 B-spline cp's
            "ref_axis_pos": 0.25,  # Define the reference axis position. 0 is the leading edge, 1 is the trailing edge.
            # Aerodynamic performance of the lifting surface at
            # an angle of attack of 0 (alpha=0).
            # These CL0 and CD0 values are added to the CL and CD
            # obtained from aerodynamic analysis of the surface to get
            # the total CL and CD.
            # These CL0 and CD0 values do not vary wrt alpha.
            "CL0": 0,  # CL of the surface at alpha=0
            "CD0":0,  # CD of the surface at alpha=0
            # Airfoil properties for viscous drag calculation
            "k_lam": 0,  # percentage of chord with laminar
            # flow, used for viscous drag
            "t_over_c": 0.12,  # thickness over chord ratio (NACA0015)
            "c_max_t": 0.303,  # chordwise location of maximum (NACA0015)
            # thickness
            "with_viscous": False,  # if true, compute viscous drag,
            "with_wave": False,
        }  # end of surface dictionary

        # Update each surface with default options
        for surface in surfaces:
            surface.update(surf_options)

        mach = 0.001
        aoa = 3.1
        beta = 0.0
        rho = 1027
        vel = 6
        re = 5.04924e+06 # With respect to chord

        dvs = self.add_subsystem("dvs", om.IndepVarComp(), promotes=["*"])
        dvs.add_output(MPhysVariables.Aerodynamics.FlowConditions.ANGLE_OF_ATTACK, val=aoa, units="deg")
        dvs.add_output(MPhysVariables.Aerodynamics.FlowConditions.YAW_ANGLE, val=beta, units="deg")
        dvs.add_output("rho", val=rho, units="kg/m**3")
        dvs.add_output(MPhysVariables.Aerodynamics.FlowConditions.MACH_NUMBER, mach)
        dvs.add_output("v", vel, units="m/s")
        dvs.add_output(MPhysVariables.Aerodynamics.FlowConditions.REYNOLDS_NUMBER, re, units="1/m")
        dvs.add_output("cg", val=np.zeros((3)), units="m")

        # Create mphys builder for aero solver
        aero_builder = AeroBuilder(surfaces)
        aero_builder.initialize(self.comm)

        # Create mesh component and connect with solver
        self.add_subsystem("mesh", aero_builder.get_mesh_coordinate_subsystem())

        # add the geometry component, we dont need a builder because we do it here.
        self.add_subsystem("geometry", OM_DVGEOCOMP(file=vsp_file, type="vsp"))
        # add pointset
        self.geometry.nom_add_discipline_coords("aero")

        self.mphys_add_scenario("cruise", ScenarioAerodynamic(aero_builder=aero_builder))
        
        self.connect(
            f"mesh.x_aero0_mesh",
            f"geometry.x_aero_in",
        )
        self.connect(
            f"geometry.x_aero0",
            f"cruise.x_aero",
       )

        # Connect dv ivc's to solver
        for dv in [
            MPhysVariables.Aerodynamics.FlowConditions.ANGLE_OF_ATTACK,
            MPhysVariables.Aerodynamics.FlowConditions.YAW_ANGLE,
            MPhysVariables.Aerodynamics.FlowConditions.MACH_NUMBER,
            MPhysVariables.Aerodynamics.FlowConditions.REYNOLDS_NUMBER,
            "rho",
            "v",
            "cg",
        ]:
            self.connect(dv, f"cruise.{dv}")

    def configure(self):
        # create geometric DV setup
        self.geometry.nom_addVSPVariable("WingGeom", "XSec_1", "Root_Chord", scaledStep=False)
        self.geometry.nom_addVSPVariable("WingGeom", "XSec_2", "Root_Chord", scaledStep=False)
        self.geometry.nom_addVSPVariable("WingGeom", "XSec_3", "Root_Chord", scaledStep=False)
        self.geometry.nom_addVSPVariable("WingGeom", "XSec_4", "Root_Chord", scaledStep=False)
        self.geometry.nom_addVSPVariable("WingGeom", "XSec_5", "Root_Chord", scaledStep=False)
        self.geometry.nom_addVSPVariable("WingGeom", "XSec_6", "Root_Chord", scaledStep=False)
        self.geometry.nom_addVSPVariable("WingGeom", "XSec_7", "Root_Chord", scaledStep=False)
        self.geometry.nom_addVSPVariable("WingGeom", "XSec_8", "Root_Chord", scaledStep=False)
        self.geometry.nom_addVSPVariable("WingGeom", "XSec_9", "Root_Chord", scaledStep=False)
        self.geometry.nom_addVSPVariable("WingGeom", "XSec_10", "Root_Chord", scaledStep=False)
        self.geometry.nom_addVSPVariable("WingGeom", "XSec_11", "Root_Chord", scaledStep=False)
        self.geometry.nom_addVSPVariable("WingGeom", "XSec_12", "Root_Chord", scaledStep=False)
        self.geometry.nom_addVSPVariable("WingGeom", "XSec_13", "Root_Chord", scaledStep=False)
        self.geometry.nom_addVSPVariable("WingGeom", "XSec_14", "Root_Chord", scaledStep=False)
        self.geometry.nom_addVSPVariable("WingGeom", "XSec_15", "Root_Chord", scaledStep=False)
        self.geometry.nom_addVSPVariable("WingGeom", "XSec_16", "Root_Chord", scaledStep=False)
        self.geometry.nom_addVSPVariable("WingGeom", "XSec_17", "Root_Chord", scaledStep=False)
        self.geometry.nom_addVSPVariable("WingGeom", "XSec_18", "Root_Chord", scaledStep=False)
        self.geometry.nom_addVSPVariable("WingGeom", "XSec_18", "Tip_Chord", scaledStep=False)
        


prob = om.Problem()
prob.model = Top()

# Set optimizer as model driver
prob.driver = om.ScipyOptimizeDriver()
prob.driver.options["debug_print"] = ["nl_cons", "objs", "desvars"]

# Setup problem and add design variables, constraint, and objective
prob.model.add_design_var("geometry.WingGeom:XSec_1:Root_Chord", lower=1e-3, upper=5.0)
prob.model.add_design_var("geometry.WingGeom:XSec_2:Root_Chord", lower=1e-3, upper=5.0)
prob.model.add_design_var("geometry.WingGeom:XSec_3:Root_Chord", lower=1e-3, upper=5.0)
prob.model.add_design_var("geometry.WingGeom:XSec_4:Root_Chord", lower=1e-3, upper=5.0)
prob.model.add_design_var("geometry.WingGeom:XSec_5:Root_Chord", lower=1e-3, upper=5.0)
prob.model.add_design_var("geometry.WingGeom:XSec_6:Root_Chord", lower=1e-3, upper=5.0)
prob.model.add_design_var("geometry.WingGeom:XSec_7:Root_Chord", lower=1e-3, upper=5.0)
prob.model.add_design_var("geometry.WingGeom:XSec_8:Root_Chord", lower=1e-3, upper=5.0)
prob.model.add_design_var("geometry.WingGeom:XSec_9:Root_Chord", lower=1e-3, upper=5.0)
prob.model.add_design_var("geometry.WingGeom:XSec_10:Root_Chord", lower=1e-3, upper=5.0)
prob.model.add_design_var("geometry.WingGeom:XSec_11:Root_Chord", lower=1e-3, upper=5.0)
prob.model.add_design_var("geometry.WingGeom:XSec_12:Root_Chord", lower=1e-3, upper=5.0)
prob.model.add_design_var("geometry.WingGeom:XSec_13:Root_Chord", lower=1e-3, upper=5.0)
prob.model.add_design_var("geometry.WingGeom:XSec_14:Root_Chord", lower=1e-3, upper=5.0)
prob.model.add_design_var("geometry.WingGeom:XSec_15:Root_Chord", lower=1e-3, upper=5.0)
prob.model.add_design_var("geometry.WingGeom:XSec_16:Root_Chord", lower=1e-3, upper=5.0)
prob.model.add_design_var("geometry.WingGeom:XSec_17:Root_Chord", lower=1e-3, upper=5.0)
prob.model.add_design_var("geometry.WingGeom:XSec_18:Root_Chord", lower=1e-3, upper=5.0)
prob.model.add_design_var("geometry.WingGeom:XSec_18:Tip_Chord", lower=1e-3, upper=5.0)

prob.model.add_design_var("angle_of_attack", lower=-10.0, upper=10.0)
prob.model.add_constraint("cruise.WingGeom.CL", equals=0.5124)
prob.model.add_constraint("cruise.WingGeom.S_ref", equals=0.0815)
prob.model.add_objective("cruise.WingGeom.CD", scaler=1e4)


prob.setup()

# Create a n^2 diagram for user to view model connections
# om.n2(prob)
prob.run_model()
prob.run_driver()
# Write optimized geometry to vsp file
#test = prob.model.geometry.nom_getDVGeo
#test.writePlot3d("opt.fmt")
#DVGeometryVSP.writeVSPFile("opt_chord.vsp3")

vspdesignvars = prob.model.geometry.nom_getDVGeo()
vspdesignvars.writeVSPFile("opt_chord.vsp3")



print("CL", prob["cruise.WingGeom.CL"][0])
print("CD", prob["cruise.WingGeom.CD"][0])

