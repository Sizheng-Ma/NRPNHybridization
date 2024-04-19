import numpy as np
import matplotlib.pyplot as plt
from Hybridization import get_abd, hyb_quantites, PNParameters, Hybridize, fix_BMS, abd_to_WM
import scri
import spherical_functions as sf
import matplotlib.pyplot as pl
from scipy.optimize import least_squares
from Hybridization import Optimize12D
from scipy.optimize import approx_fprime
import quaternion
from Hybridization import SquaredError
from Hybridization import Align

OptArg = 1
maxiter = 30

cce_dir = "/gpfs/sma2/memory/gw150914/cce"
abd, t0 = get_abd(cce_dir, None)

t_end = 3000 + t0
length = 2500
hyb = hyb_quantites(t_end, length)

data_dir = "/gpfs/sma2/memory/gw150914/horizon"
PN = PNParameters(data_dir, hyb, t0)

if hyb.PNIter>maxiter:
    hyb.PNIter = 10

out_name = "Output.npz"
nOrbits = None

W_NR = abd_to_WM(abd, lmin=0)

print("hi")

N_itr_maxes = {"superrest":1, "CoM_transformation": 10, "rotation" : 10, "supertranslation":10}
#trans = abd.map_to_superrest_frame(t_0=hyb.t_start+hyb.length/2, N_itr_maxes=N_itr_maxes,print_conv=True)

res=scri.bms_transformations.BMSTransformation()
res.from_file("full_bms.h5")
abd_smaller = abd.interpolate(abd.t[::10])
trans = abd_smaller.transform(supertranslation=res.supertranslation,frame_rotation=res.frame_rotation.components,boost_velocity=res.boost_velocity)

W_NR = W_NR.transform(supertranslation=res.supertranslation,frame_rotation=res.frame_rotation.components,boost_velocity=res.boost_velocity)

hyb.get_window_NR(W_NR)

PN.Parameterize_to_Physical(hyb)

minima, W_PN, chiA, chiB = Align(PN, hyb)
logR_delta = np.append(minima.x[0], minima.x[1:] + hyb.omega_mean*minima.x[0]/2)

PN.PhyParas = np.append(PN.PhyParas, logR_delta)

PN.Physical_to_Parameterize(hyb)

scale = np.array([0.05,0.02,0.1,0.1,np.pi*2,np.pi*2,np.pi,np.pi/4,np.pi/hyb.omega_i/2.0,np.pi/4,np.pi/4,np.pi/4])

lowbound12D = PN.OptParas - scale

upbound12D = PN.OptParas + scale

minima12D = least_squares(Optimize12D, PN.OptParas, bounds=(lowbound12D, upbound12D), ftol=3e-15, xtol=3e-15, gtol=3e-15, x_scale='jac', args=(PN, hyb))

minima12D.jac = approx_fprime(minima12D.x, Optimize12D, np.full_like(minima12D.x, 1.49e-8), PN, hyb)

PN.OptParas = minima12D.x

PN.Parameterize_to_Physical(hyb)

hyb.t_PNStart, hyb.t_PNEnd = -80000, 1000 - hyb.t_start

minima, W_PN, chiA, chiB = Align(PN, hyb)

t_delta = minima.x[0]
logR_delta = np.append(0.0, minima.x[1:] + hyb.omega_mean*minima.x[0]/2)

PN.PhyParas = np.append(PN.PhyParas, logR_delta)

PN.Physical_to_Parameterize(hyb)

R_delta = np.exp(quaternion.from_float_array(logR_delta))

W_PN.t = W_PN.t - t_delta

W_NR = scri.rotate_physical_system(W_NR, R_delta)
chiA = R_delta*chiA*R_delta.conjugate()
chiB = R_delta*chiB*R_delta.conjugate()



#padding_time = 250
#t_0=hyb.t_start+hyb.length/2
#
#abd_interp = abd.interpolate(
#    abd.t[np.argmin(abs(abd.t - (t_0 - (padding_time + 200)))) : np.argmin(abs(abd.t - (t_0 + (padding_time + 200)))) + 1]
#)
#
#time_translation = scri.bms_transformations.BMSTransformation(supertranslation=[sf.constant_as_ell_0_mode(t_0)])
#
#BMS_transformation = time_translation * scri.bms_transformations.BMSTransformation().reorder(
#    ["supertranslation", "frame_rotation", "boost_velocity"]
#)
#
#abd_interp_prime = abd_interp.transform(
#    supertranslation=BMS_transformation.supertranslation,
#    frame_rotation=BMS_transformation.frame_rotation.components,
#    boost_velocity=BMS_transformation.boost_velocity,
#)
#
#new_transformation, CoM_rel_errs = scri.asymptotic_bondi_data.map_to_superrest_frame.com_transformation_to_map_to_superrest_frame(abd_interp_prime, N_itr_max=2, rel_err_tol=1e-2,print_conv=True)
#
#BMS_transformation = (new_transformation * BMS_transformation).reorder(
#    ["supertranslation", "frame_rotation", "boost_velocity"]
#)
#
#BMS_transformation.to_file("com_trans.h5")
