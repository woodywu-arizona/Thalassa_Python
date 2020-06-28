import pytha
import numpy as np
import array
import time
# theoretical objects theo in decimal code 
Satnum = 841001 # theoretical objects start with 841XXX XXX for real number
#  propagation start date
MJD0 = 59000

#  initial conditions
COE0 = [7086.8184,0.004,60,90,90,0.01] # comparable frz-osc case
#COE0 = [7086.8184,0.004,3,90,90,0.01]  # noncompara

# // spacecraft physical parameters
SCMass = 13
ADrag  = 0.196
ASRP   = 0.196
CD     = 1.28
CR     = 1

# // propagation time span and step
MU = 3.986004414498200E+05
tspan = 540.0
tstep = 8.0

# // pure zonal physical model
insgrav = 1
isun    = 1
imoon   = 1
idrag   = 4
iF107   = 1
iSRP    = 2
iephem  = 1
gdeg    = 15
gord    = 15

# // integration tolerance
tol = 1.0e-12

# // no moon collision
imcoll = 0

# // equations of montion
# 1 = Cowell, 2 = EDromo(t), 3 = EDromo(c), 4 = EDromo(l),
# 5 = KS t), 6 = KS (l), 7 = Sti-Sche (t), 8 = Sti-Sche (l)
eqs = 7

# external definition
npts = 0
mxpts = round(tspan/tstep) + 1
print(mxpts)
cart_out = [[0 for x in range(mxpts)] for x in range(7)]
orbs_out = [[0 for x in range(mxpts)] for x in range(7)]
tag = 1
exitcode = 0
rmxstep = 200000

file_output_dir = '/Users/woodywu/Desktop/Research/Proper_Element/Proper_Element_Data_Association_Study/Data_Theo'
file_output_nam_orbs = 'lngtrm_TheoCompara_orbs.csv'
file_output_nam_cart = 'lngtrm_TheoCompara_cart.csv'
file_output_orbs = file_output_dir+'/'+file_output_nam_orbs
file_output_cart = file_output_dir+'/'+file_output_nam_cart

# call python thalassa module
[cart_out,orbs_out]=pytha.thalassasub(MJD0,COE0,tspan,tstep, insgrav, isun, imoon, idrag, iF107, iSRP, iephem, gdeg, gord, rmxstep, tol, imcoll, eqs, SCMass, ADrag, ASRP, CD, CR, mxpts, npts, tag, exitcode)

# append physical parameters
length = orbs_out.shape[1] # col numbers
physc_param = np.asarray([[SCMass]*length,[ADrag]*length,[ASRP]*length,[CD]*length,[CR]*length,[Satnum]*length])
cart_out_physc = np.append(cart_out,physc_param,axis=0)
orbs_out_physc = np.append(orbs_out,physc_param,axis=0)

print(orbs_out.shape)
print(orbs_out[0],orbs_out[1])
np.savetxt(file_output_orbs,orbs_out_physc.transpose(),delimiter=',')
np.savetxt(file_output_cart,cart_out_physc.transpose(),delimiter=',')

