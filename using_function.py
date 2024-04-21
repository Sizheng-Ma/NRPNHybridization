import scri

from Hybridization import get_abd, hyb_quantites, PNParameters, Hybridize, fix_BMS, abd_to_WM, MyHybridize

maxiter = 3


cce_dir = "/gpfs/sma2/memory/gw150914/cce"
_, t0 = get_abd(cce_dir, None)

t_end = 3000 + t0
length = 2500

hyb = hyb_quantites(t_end, length)

data_dir = "/gpfs/sma2/memory/gw150914/horizon"
PN = PNParameters(data_dir, hyb, t0)

while hyb.PNIter<=maxiter:
        print("PNIter=: ", hyb.PNIter)
        W_NR, W_PN, minima12D = MyHybridize(cce_dir, hyb, PN, OptimizePNParas = 1, verbose=True)
        if hyb.PNIter >= 2 and abs(hyb.cost[-1]-hyb.cost[-2])/hyb.cost[-1]<1e-2 and abs(hyb.cost[-1]-hyb.cost[-3])/hyb.cost[-1]<1e-2:
                hyb.PNIter = maxiter + 1
        else:
                hyb.PNIter += 1


#abd_original, t0 = get_abd(cce_dir, None)
#abd = abd_original.interpolate(abd_original.t[::10])
#
##within fix_BMS
#from Hybridization import PostNewtonian
#Phys = PN.PhyParas
#W_PN = PostNewtonian.PNWaveform(
#    Phys[0], np.copy(hyb.omega_i)*Phys[1], Phys[2:5], Phys[5:8], PN.frame_i, np.copy(hyb.t_start)/Phys[1],
#    t_PNStart=hyb.t_PNStart, t_PNEnd=hyb.t_PNEnd
#)
#W_PN.t = W_PN.t*Phys[1]
##W_PN.data = np.append(0.0*W_PN.data[:,0:4], np.copy(W_PN.data), axis=1)
##W_PN.ells = 0,8
#W_PN.dataType = scri.h
#W_PN_PsiM = PostNewtonian.PNWaveform(
#    Phys[0], np.copy(hyb.omega_i)*Phys[1], Phys[2:5], Phys[5:8], PN.frame_i, np.copy(hyb.t_start)/Phys[1],
#    t_PNStart=hyb.t_PNStart, t_PNEnd=hyb.t_PNEnd, datatype="Psi_M"
#)
#W_PN_PsiM.t = W_PN_PsiM.t*Phys[1]
