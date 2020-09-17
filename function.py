import numpy as np
import math as m
def my_distance(drij):
    return m.sqrt(drij[0] * drij[0] + drij[1] * drij[1] + drij[2] * drij[2])

def my_disp_in_box(drij, lbox):
    drij_M = drij
    lbox_inv  = 1.0 / lbox
    lbox_half = lbox * 0.5

    drij_M[0] = drij_M[0] + lbox_half
    drij_M[1] = drij_M[1] + lbox_half
    drij_M[2] = drij_M[2] + lbox_half

    if drij_M[0] < 0.0:
        n = int(drij_M[0] * lbox_inv)
        drij_M[0] = drij_M[0] - (n - 1.0) * lbox
    if drij_M[0] >= lbox:
        n = int(drij_M[0] * lbox_inv)
        drij_M[0] = drij_M[0] - n * lbox
    if drij_M[1] < 0.0:
        n = int(drij_M[1] * lbox_inv)
        drij_M[1] = drij_M[1] - (n - 1.0) * lbox
    if drij_M[1] >= lbox:
        n = int(drij_M[1] * lbox_inv)
        drij_M[1] = drij_M[1] - n * lbox
    if drij_M[2] < 0.0:
        n = int(drij_M[2] * lbox_inv)
        drij_M[2] = drij_M[2] - (n - 1.0) * lbox
    if drij_M[2] >= lbox:
        n = int(drij_M[2] * lbox_inv)
        drij_M[2] = drij_M[2] - n * lbox
    drij_M[0] = drij_M[0] - lbox_half
    drij_M[1] = drij_M[1] - lbox_half
    drij_M[2] = drij_M[2] - lbox_half
    return drij_M

def my_potential_energy_i(i, pos, lbox, rc):
    epsilon = 1.0
    sigma   = 1.0
    rij   = np.zeros(3)
    n       = len(pos)
    U       = 0.0
    div_c   = sigma / rc
    U_c     = 4.0 * epsilon * (div_c)**6 * (div_c**6 - 1)
    for j in range(0, n):
        if j != i:
            rij[0] = -pos[i][0] + pos[j][0]
            rij[1] = -pos[i][1] + pos[j][1]
            rij[2] = -pos[i][2] + pos[j][2]
            rij    = my_disp_in_box(rij,lbox)
            dij    = my_distance(rij)
            if dij < rc:
                div     = sigma / dij
                U       = 4.0 * epsilon * (div)**6 * (div**6 - 1) + U - U_c
    return U

# def my_potential_energy(pos, lbox, rc):
#     sigma   = 1.0
#     epsilon = 1.0
#     n       = len(pos)
#     rij     = np.zeros(3)
#     U       = 0.0
#     div_c   = sigma / rc
#     U_c     = 4.0 * epsilon * (div_c)**6 * (div_c**6 - 1)
#     f
# r i in range(0,n):
#         for j in range(i + 1,n):
#             rij[0] = -pos[i][0] + pos[j][0]
#             rij[1] = -pos[i][1] + pos[j][1]
#             rij[2] = -pos[i][2] + pos[j][2]
#             rij    = my_disp_in_box(rij,lbox)
#             dij    = my_distance(rij)
#             if dij == 0:
#                 print dij
#                 sys.exit("particles overlap :-(")
#             if dij < rc:
#                 div     = sigma / dij
#                 U       = 4.0 * epsilon * (div)**6 * (div**6 - 1) + U - U_c
#     return U

def my_loga_symmetric(vold, vnew, beta):
    new = np.exp(-beta * vnew)
    old = np.exp(-beta * vold)
    return np.log(new/old)

# def my_mc_sweep(pos, lbox, rc, beta, eta, acc_check):
    
#     naccept = 0
#     de = 0
#     natom, ndim = pos.shape
#     pos_new = np.zeros((natom, ndim))
#     # pos_i   = pos
#     # pos_new = pos
#     for i in range(0, natom):
#         for j in range(0, ndim):
#             pos_new[i][j] = pos[i][j]
#     U_old     = my_potential_energy(pos, lbox, rc)
#     for i in range(0, natom):
#         U_old_i   = my_potential_energy_i(i, pos, lbox, rc)
#         for j in range(0, ndim):
#             pos_new[i][j] = pos[i][j] + eta[i][j]
#         U_new_i = my_potential_energy_i(i, pos_new, lbox, rc)
#         pos_new[i] = my_disp_in_box(pos_new[i], lbox)
#         q     = np.exp(beta * (U_old_i - U_new_i))
#         # print U_old_i
#         # print q
#         if q > acc_check[i]:
#             naccept = naccept + 1
#             for k in range(0, ndim):
#                 pos[i][k] = pos_new[i][k]
#             de  = de + (U_new_i - U_old_i)
#             # U_old_i = U_new_i
#     # de = my_potential_energy(pos, lbox, rc) - U_old
#     return naccept, de

def my_mc_sweep(pos, lbox, rc, beta, eta, acc_check):
    
    naccept     = 0
    de          = 0
    natom, ndim = pos.shape
    pos_new     = np.zeros((natom, ndim))
    for i in range(0, natom):
        U_old = my_potential_energy_i(i, pos, lbox, rc)
        for k in range(0, natom):
            if k == i:
                pos_new[k][0] = pos[k][0] + eta[k][0]
                pos_new[k][1] = pos[k][1] + eta[k][1]
                pos_new[k][2] = pos[k][2] + eta[k][2]
            else:
                pos_new[k][0] = pos[k][0] 
                pos_new[k][1] = pos[k][1] 
                pos_new[k][2] = pos[k][2] 
        pos_new[i] = my_disp_in_box(pos_new[i], lbox)
        U_new      = my_potential_energy_i(i, pos_new, lbox, rc)
        q          = np.exp( beta * (U_old - U_new))
        if q > acc_check[i]:
            naccept = naccept + 1
            de      = de + (U_new - U_old)
            # U_old   = U_new
            pos[i][0] = pos_new[i][0] 
            pos[i][1] = pos_new[i][1] 
            pos[i][2] = pos_new[i][2] 
    return naccept, de

def my_force_i(i, pos, lbox,rc):
    force = np.zeros(3)
    rij   = np.zeros(3)
    n     = len(pos)
    epsilon = 1.0
    sigma   = 1.0
    epsilon_24 = 24.0 * epsilon

    for j in range(0,n):
        if j != i:
            rij[0] = -pos[i][0] + pos[j][0]
            rij[1] = -pos[i][1] + pos[j][1]
            rij[2] = -pos[i][2] + pos[j][2]
            rij    = my_disp_in_box(rij,lbox)
            dij    = my_distance(rij)
            if dij < rc:
                div = (1.0 / dij)
                com = -48.0 * sigma**12 * div**14 + 24.0 * sigma**6 * div**8
                # com = div * div**3 * (2 * div**3 - 1.0) 
                force[0] = force[0] + com * rij[0]
                force[1] = force[1] + com * rij[1]
                force[2] = force[2] + com * rij[2]
    force[0] = force[0] * epsilon
    force[1] = force[1] * epsilon
    force[2] = force[2] * epsilon
    return force
# def my_loga_smart(fold, fnew, eta0, sigma, beta):

#     tau = sigma**2
#     new = np.exp(-tau * beta * (fnew[0] * eta0[0] + fnew[1] * eta0[1] + fnew[2] * eta0[2]))
#     old = np.exp( tau * beta * (fold[0] * eta0[0] + fold[1] * eta0[1] + fold[2] * eta0[2]))
#     return np.lslog(new/old)
def my_loga_smart(fold, fnew, eta0, sigma, beta):
    tau = sigma**2
    xaold = tau*beta*fold
    xnewa = tau*beta*fnew
    eta1 = -(eta0 + xaold + xnewa)
    return (np.dot(eta0,eta0) - np.dot(eta1,eta1))/(2*tau)

# def my_mc_sweep_smart(pos, lbox, rc, beta, eta, acc_check, sig):
#     naccept     = 0
#     de          = 0
#     tau         = sig**2
#     taub        = tau * beta
#     natom, ndim = pos.shape
#     pos_new     = np.zeros((natom, ndim))
#     F_old       = np.zeros((natom, 1))
#     F_new       = np.zeros((natom, 1))
#     for i in range(0, natom):
#         U_old     = my_potential_energy_i(i, pos, lbox, rc)
#         F_old     = my_force_i(i, pos, lbox,rc)
#         for k in range(0, natom):
#             if k == i:
#                 pos_new[k][0] = pos[k][0] + eta[k][0] #+ taub * F_old[0]
#                 pos_new[k][1] = pos[k][1] + eta[k][1] #+ taub * F_old[1]
#                 pos_new[k][2] = pos[k][2] + eta[k][2] #+ taub * F_old[2]
#             else:
#                 pos_new[k][0] = pos[k][0] 
#                 pos_new[k][1] = pos[k][1] 
#                 pos_new[k][2] = pos[k][2] 
#         pos_new[i] = my_disp_in_box(pos_new[i], lbox)
#         U_new = my_potential_energy_i(i, pos_new, lbox, rc)
#         F_new  = my_force_i(i, pos_new, lbox,rc)
#         # F        = F_old + F_new
#         # R        = 2 * eta[i] + taub * (F_new - F_old)
#         # Fdotr    = F[0] * R[0] + F[1] * R[1] + F[2] * R[2] 
#         q = np.exp( my_loga_symmetric(U_old, U_new, beta))# + my_loga_smart(F_old, F_new, eta[i], sig, beta)) 
#         print q
#         if q > acc_check[i]:
#             naccept = naccept + 1
#             de      = de + (U_new - U_old)
#             # U_old   = U_new
#             pos[i][0] = pos_new[i][0] 
#             pos[i][1] = pos_new[i][1] 
#             pos[i][2] = pos_new[i][2] 
#     return naccept, de
#(29, -47.848431989117778)

def my_mc_sweep_smart(pos, lbox, rc, beta, eta, acc_check, sig):
    
    naccept     = 0
    de          = 0
    tau         = sig**2
    taub        = tau * beta

    natom, ndim = pos.shape
    pos_new     = np.zeros((natom, ndim))
    F_old       = np.zeros((natom, 1))
    F_new       = np.zeros((natom, 1))

    for i in range(0, natom):
        U_old     = my_potential_energy_i(i, pos, lbox, rc)
        F_old     = my_force_i(i, pos, lbox,rc)
        for k in range(0, natom):
            if k == i:
                pos_new[k][0] = pos[k][0] + eta[k][0] + taub * F_old[0] 
                pos_new[k][1] = pos[k][1] + eta[k][1] + taub * F_old[1]
                pos_new[k][2] = pos[k][2] + eta[k][2] + taub * F_old[2]
            else:
                pos_new[k][0] = pos[k][0] 
                pos_new[k][1] = pos[k][1] 
                pos_new[k][2] = pos[k][2] 
        pos_new[i] = my_disp_in_box(pos_new[i], lbox)
        U_new      = my_potential_energy_i(i, pos_new, lbox, rc)
        F_new  = my_force_i(i, pos_new, lbox,rc)
        q = np.exp( my_loga_symmetric(U_old, U_new, beta)+ my_loga_smart(F_old, F_new, eta[i], sig, beta)) 
        print q
        if q > acc_check[i]:
            naccept = naccept + 1
            de      = de + (U_new - U_old)
            # U_old   = U_new
            pos[i][0] = pos_new[i][0] 
            pos[i][1] = pos_new[i][1] 
            pos[i][2] = pos_new[i][2] 
    return naccept, de
def my_pos_ini_2(natom,lbox):
    r = np.zeros((natom,3))
    n = int(natom**(1.0/3.0)) + 1

    delta      = lbox / n
    delta_half = delta * 0.5
    lbox_half  = lbox  * 0.5
    index = 0
    for i in range(0,n):
        x = delta * i
        for j in range(0,n):
            y = delta * j
            for k in range(0,n):
                z = delta * k
                r[index][0] = x - lbox_half + delta_half
                r[index][1] = y - lbox_half + delta_half
                r[index][2] = z - lbox_half + delta_half
                index = index + 1
    return r
def my_potential_energy_total(pos, lbox, rc):
    sigma   = 1.0
    epsilon = 1.0
    n       = len(pos)
    rij     = np.zeros(3)
    U       = 0.0
    div_c   = sigma / rc
    U_c     = 4.0 * epsilon * (div_c)**6 * (div_c**6 - 1)
    for i in range(0,n):
        for j in range(i + 1,n):
            rij[0] = -pos[i][0] + pos[j][0]
            rij[1] = -pos[i][1] + pos[j][1]
            rij[2] = -pos[i][2] + pos[j][2]
            rij    = my_disp_in_box(rij,lbox)
            dij    = my_distance(rij)
            if dij == 0:
                print dij
                sys.exit("particles overlap :-(")
            if dij < rc:
                div     = sigma / dij
                U       = 4.0 * epsilon * (div)**6 * (div**6 - 1) + U - U_c
    return U

def my_pair_correlation(dists, natom, nbins, dr, lbox):
    gr       = np.zeros(nbins)
    V        = np.zeros(nbins)
    PI       = 4.0 * 3.14159265359 / 3.0
    counts   = np.zeros(nbins, dtype=int)
    n        = len(dists)
    V_max    = (lbox)**3
    N_Max    = natom * (natom - 1.0) / 2.0
    density_IV = V_max / N_Max

    for i in range(0,nbins):
        min_r = i * dr
        max_r = min_r + dr
        V[i]  = PI * (max_r**3 - min_r**3)
        for j in range(0,n):
            if (dists[j] < max_r and dists[j] > min_r):
                counts[i] = counts[i] + 1
        gr[i] = counts[i] * density_IV / V[i]
    return gr
def my_legal_kvecs(maxn, lbox):
    kvecs = np.zeros(((maxn + 1)**3 , 3))
    PI    = 2.0 * 3.14159265359 / lbox
    index = 0
    for i in range(0,maxn + 1):
        for j in range(0,maxn + 1):
            for k in range(0,maxn + 1):
                kvecs[index][0] = i * PI
                kvecs[index][1] = j * PI
                kvecs[index][2] = k * PI
                index = index + 1
    return np.array(kvecs)
def my_calc_sk(kvecs, pos):
    nk   = len(kvecs)
    nr   = len(pos)
    sk   = np.zeros(nk)
    rhok = np.zeros(nk, dtype=complex)
    rhok = my_calc_rhok(kvecs, pos)
    for i in range(0,nk):
        sk[i] = (rhok[i] * rhok[i].conjugate()) / nr
    return sk
def my_calc_rhok(kvecs, pos):
    nk = len(kvecs)
    nr = len(pos)
    rhok = np.zeros(nk, dtype=complex)
    for i in range(0,nk):
        S = 0 + 0j
        for j in range(0,nr):
            a = kvecs[i][0] * pos[j][0] + kvecs[i][1] * pos[j][1] + kvecs[i][2] * pos[j][2]
            S = S + np.exp(-a * 1j)
        rhok[i] = S
    return rhok

def my_mean(a): 
    S = 0.0
    n = len(a)
    for i in range(0,n):
        S = a[i] + S
    return S/n
def my_std(a): 
    mean = 0.0
    SD   = 0.0
    n = len(a)
    for i in range(0,n):
        mean = a[i] + mean
    mean = mean/n
    for i in range(0,n):
        SD = (a[i] - mean)**2 + SD
    SD = m.sqrt(SD/(n - 1.0))
    return SD
def my_var(a): 
    mean = 0.0
    SD   = 0.0
    n = len(a)
    for i in range(0,n):
        mean = a[i] + mean
    mean = mean/n
    for i in range(0,n):
        SD = (a[i] - mean)**2 + SD
    SD = (SD/(n - 1.0))
    return SD
def my_actime(a): 
    mean     = my_mean(a)
    SD       = my_std(a)
    N  = len(a)
    St = [0 for x in range(N)] 
    min_t = N
    for t in range(0,N):
        Si = 0.0
        for i in range(0,N - t):
            Si = Si + (a[i] - mean) * (a[i + t] - mean)
        St[t] = Si / ((SD**2) * (N - t))
        if St[t]<0:
            if t<min_t:
                min_t = t
    t_cutoff = min_t
    SUM      = 0.0
    for t in range(1,t_cutoff):
        SUM = St[t] + SUM
    return (1 + 2 * SUM)
def my_stderr(a):
    N   = len(a)
    SD  = my_std(a)
    K   = my_actime(a)
    return SD / m.sqrt(N/K)
def my_error_unco(a):
    mean = my_mean(a)
    N    = len(a)
    S    = 0.0
    for i in range(0,N):
        S = S + (a[i] - mean)**2
    S = m.sqrt(S / (N * (N - 1)))
    return S
def my_error_co(a):
    mean = my_mean(a)
    K    = my_actime(a)
    N    = len(a)
    S    = 0.0
    for i in range(0,N):
        S = S + (a[i] - mean)**2
    S = m.sqrt(K * S / (N * (N - 1)))
    return S
