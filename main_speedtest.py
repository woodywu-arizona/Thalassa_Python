import pytha
import numpy as np
import array
import time

#  propagation start date
MJD0 = 58849

#  initial conditions
COE0 = [20000,0.65,87.9,0,45,0]

# // spacecraft physical parameters
SCMass = 150
ADrag  = 1.84
ASRP   = 1.84
CD     = 1.28
CR     = 1

# // propagation time span and step
MU = 3.986004414498200E+05
tspan = 100
tstep = 1.0

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

# // equations of montion - EDromo(t)
eqs = 2

# external definition
npts = 0
mxpts = round(tspan/tstep) + 1
cart_out = [[0 for x in range(mxpts)] for x in range(7)]
orbs_out = [[0 for x in range(mxpts)] for x in range(7)]
tag = 1
exitcode = 0
rmxstep = 1000

N = 100
start_time = time.time()
for i in range(N):
    # call python thalassa module
    [cart_out,orbs_out]=pytha.thalassasub(MJD0,COE0,tspan,tstep, insgrav, isun, imoon, idrag, iF107, iSRP, iephem, gdeg, gord, rmxstep, tol, imcoll, eqs, SCMass, ADrag, ASRP, CD, CR, mxpts, npts, tag, exitcode)
total_time = time.time() - start_time
averg_time = total_time / N

print(averg_time)

#print(cart_out[0])
#print(cart_out[0,15])
##print(cart_out[1])
#print(orbs_out[3])
#print(orbs_out[1,3])
#print(type(orbs_out[0,0]))
#print(orbs_out.shape)

