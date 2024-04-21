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

W_NR, W_PN, minima12D = MyHybridize(cce_dir, hyb, PN, OptimizePNParas = 1, verbose=True)
hyb.PNIter += 1
