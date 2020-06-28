import pytha
import numpy as np
import array
import time,os
from multiprocessing import Pool

# environment setup
MU = 3.986004414498200E+05
# physical model
insgrav = 1
isun    = 0
imoon   = 0
idrag   = 0
iF107   = 0
iSRP    = 0
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

# // propagation time span and step
tspan = 365*20
tstep = 2.0
# external definition
npts = 0
mxpts = int(tspan/tstep + 1)
cart_out = [[0 for x in range(mxpts)] for x in range(7)]
orbs_out = [[0 for x in range(mxpts)] for x in range(7)]
tag = 1
exitcode = 0
rmxstep = 200000

file_input_dir = '/Users/woodywu/Desktop/Research/Proper_Element/Proper_Element_Data_Association_Study/Data_Theo'
file_input_nam = 'lngtrm_TheoCompara_orbs'
file_input = file_input_dir+'/'+file_input_nam+'.csv'
data = np.genfromtxt(file_input,delimiter=',')

start_time = time.time()
file_output_dir = file_input_dir
file_output_subdir = 'Trajectory_'+file_input_nam
file_output_general_dir = file_output_dir+'/'+file_output_subdir
if not os.path.isdir(file_output_general_dir):
    os.mkdir(file_output_general_dir)

def trajectory_from_row(row):
    print(row[-1],'\t',row[0])
    MJD0,sma,ecc,inc,raan,aop,ma,SCMass,ADrag,ASRP,CD,CR,Satnum = row
    COE0 = [sma,ecc,inc,raan,aop,ma]
    [cart_out,orbs_out]=pytha.thalassasub(MJD0,COE0,tspan,tstep, insgrav, isun, imoon, idrag, iF107, iSRP, iephem, gdeg, gord, rmxstep, tol, imcoll, eqs, SCMass, ADrag, ASRP, CD, CR, mxpts, npts, tag, exitcode)
    # write out
    file_output_nam_orbs = str(int(Satnum))+'_'+str(int(MJD0))+'_orbs.csv'
    file_output_nam_cart = str(int(Satnum))+'_'+str(int(MJD0))+'_cart.csv'
    file_output_orbs = file_output_general_dir+'/'+file_output_nam_orbs
    file_output_cart = file_output_general_dir+'/'+file_output_nam_cart
    np.savetxt(file_output_orbs,orbs_out.transpose(),delimiter=',')
    np.savetxt(file_output_cart,cart_out.transpose(),delimiter=',')
    return

with Pool(10) as p:
    p.map(trajectory_from_row,data)

total_time = time.time() - start_time
print(total_time)


