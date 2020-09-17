import math      as m
import function as fun
import numpy     as np
import sys
import random
natom     = 64
mass      = 48
lbox      = 4.0 
rc        = 2.0 
nstep     = 5000
T_i       = 3.0
beta      = 1.0/T_i
# T         = 1.0
sigma_g   = 0.07
taub      = sigma_g * beta
r         = np.zeros((natom,3))
eta       = np.zeros((natom,3))
total     = 0
R         = np.zeros(3)
r         = fun.my_pos_ini_2(natom,lbox)
F_old       = np.zeros((natom, 1))
F_new       = np.zeros((natom, 1))

OUTPUT   = open("positions_"     + str(T_i) + ".txt" , "w")
OUTPUT1  = open("positions_VMD_" + str(T_i) + ".XYZ" , "w")
OUTPUT2  = open("energy_"        + str(T_i) + ".txt" , "w")
naccept     = 0
de          = 0
natom, ndim = r.shape
r_new     = np.zeros((natom, 3))

for t in range(0,nstep):
	if t % 100 == 0:
		print t 

	OUTPUT1.write('{}  \n'.format(natom))
	OUTPUT1.write('{}  \n'.format(t))

	for i in range(0,natom):
		OUTPUT.write('{0:.8f}  ' .format(r[i][0]))
		OUTPUT.write('{0:.8f}  ' .format(r[i][1]))
		OUTPUT.write('{0:.8f}\n' .format(r[i][2]))

	for i in range(0,natom):
		OUTPUT1.write("C  ")
		OUTPUT1.write('{0:.8f}  '.format(r[i][0]))
		OUTPUT1.write('{0:.8f}  '.format(r[i][1]))
		OUTPUT1.write('{0:.8f}\n'.format(r[i][2]))

	eta = np.random.normal(0, sigma_g, [natom, 3])
	# print eta
	OUTPUT2.write('{}       '.format(t))
	OUTPUT2.write('{0:.8f}\n'.format(fun.my_potential_energy_total(r, lbox, rc)))

	for i in range(0, natom):
		total = total + 1
		U_old = fun.my_potential_energy_i(i, r, lbox, rc)
		F_old = fun.my_force_i(i, r, lbox,rc)
		for k in range(0, natom):
			if k == i:
				r_new[k][0] = r[k][0] + eta[k][0] + taub * F_old[0] 
				r_new[k][1] = r[k][1] + eta[k][1] + taub * F_old[1]
				r_new[k][2] = r[k][2] + eta[k][2] + taub * F_old[2]
			else:
				r_new[k][0] = r[k][0] 
				r_new[k][1] = r[k][1] 
				r_new[k][2] = r[k][2] 
		r_new[i] = fun.my_disp_in_box(r_new[i], lbox)
		U_new    = fun.my_potential_energy_i(i, r_new, lbox, rc)
		F_new  	 = fun.my_force_i(i, r_new, lbox,rc)
		q = np.exp( fun.my_loga_symmetric(U_old, U_new, beta) + fun.my_loga_smart(F_old, F_new, eta[i], sigma_g, beta)) 
		if q > random.uniform(0, 1):
			naccept = naccept + 1
			de      = de + (U_new - U_old)
            # U_old   = U_new
			r[i][0] = r_new[i][0] 
			r[i][1] = r_new[i][1] 
			r[i][2] = r_new[i][2] 
print (float(naccept)/float(total))

