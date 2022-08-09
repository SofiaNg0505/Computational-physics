
#%%

from sre_constants import SUCCESS
import numpy 
import bspline
import bspline.splinelab as splinelab
import matplotlib.pyplot as plt
import numpy
%matplotlib inline
import time
import math
from numba import njit

from __future__ import division
from pylab import *
from scipy.special.orthogonal import p_roots

r = 0.00000001
R_max = 10
e = 1
p = 4-1             # order of spline i.e polynomial order 3
## Spline setup and evaluationXS
order_of_poly = 6
nknots =   51 # number of knots to generate (here endpoints count only once) endpoints = 4+4
knots = numpy.linspace(r,R_max,nknots)  # create a knot vector without endpoint repeats
k     = splinelab.augknt(knots, p)  # add endpoint repeats as appropriate for spline order p. THIS IS T
B     = bspline.Bspline(k, p)       # create spline basis of order p on knots k
tau = numpy.linspace(r,R_max,501)

@njit
def trapezoidal(f_a, f_b, a, b ):
    function =  ((b-a))*(f_a+ f_b)*0.5
    return function

def assignment5(Q,V_ee,l,n,order_of_poly,R_max, e ,p, nknots,tau): 
    knots = numpy.linspace(r,R_max,nknots)  # create a knot vector without endpoint repeats
    k     = splinelab.augknt(knots, p)  # add endpoint repeats as appropriate for spline order p. THIS IS T
    B     = bspline.Bspline(k, p)       # create spline basis of order p on knots k

    def polynom_B(B,r,w,a,b,i,j):
        spline_1 = B.collmat(r*(b-a)*0.5 + (a+b)*0.5)[i]
        spline_2 = B.collmat(r*(b-a)*0.5 + (a+b)*0.5)[j]
        function= w*(spline_1*spline_2)
        return function
    def gauss(B,n,a,b,i ,j):
        [r,w] = p_roots(n+1)
        poly = 0
        for order in range(len(r)):
            poly = polynom_B(B,r[order],w[order], a,b,i,j) + poly
        G=(b-a)*0.5*poly
        return G
    B_matrix=numpy.zeros((nknots-2 + p-1, nknots-2 + p-1))
    for i in range(len(B_matrix)):
        for j in range(len(B_matrix)):
            if abs(i-j)>(2*p+1)+1:
                a  =1
            if abs(i-j)<2*p+1:
                F = 0
                for position in range(len(tau)-1):
                    F=gauss(B,order_of_poly,tau[position],tau[position+1],i+1,j+1) + F
                    B_matrix[i,j] = F


    H_matrix = numpy.zeros((nknots-2 + p-1, nknots-2 + p-1))
    B_prim= numpy.array( [( B.collmat(t,deriv_order=1) ) for t in tau], dtype=numpy.float64 )
    B = numpy.array( [( B.collmat(t ) )for t in tau], dtype=numpy.float64 )
    B_prim = numpy.delete(B_prim, -1 , axis = 1)
    B_prim = numpy.delete(B_prim, 0 , axis = 1)
    B= numpy.delete(B, -1 , axis = 1)
    B= numpy.delete(B, 0 , axis = 1)
    
    for i in range(len(H_matrix)):
        for j in range(len(H_matrix)):
           ## print(i,j)
            term_1 =0
            term_2 = 0
            term_3 = 0
            term_4 = 0
            if abs(i-j)>(2*p+1)+1:
                a  =1
            if abs(i-j)<2*p+1:

                for position in range(len(tau)-1):
       
                    term_1 =trapezoidal(0.5*(B_prim[position,i]*B_prim[position,j]),0.5*B_prim[position+1,i]*B_prim[position+1,j],tau[position],tau[position+1]) + term_1
   
                    term_2 = trapezoidal(l*(l+1)/(tau[position]**2) *  0.5*B[position,i]*B[position,j],l*(l+1)/(tau[position+1]**2) * 0.5*B[position+1,i]*B[position+1,j],tau[position],tau[position+1]) + term_2
  
                    term_3 = trapezoidal(-(Q/(tau[position])) *B[position,i] * B[position,j],-(Q/(tau[position+1])) *B[position+1,i] * B[position+1,j],tau[position],tau[position+1]) + term_3

                    term_4 = trapezoidal(B[position,i]*V_ee[position]*B[position,j], B[position+1,i]*V_ee[position+1]*B[position+1,j],tau[position],tau[position+1])+term_4
                H_matrix[i,j] = term_1+term_2+term_3+term_4

    
    from scipy.linalg import eigh
    eigenvalues, eigenev= eigh(H_matrix,B_matrix)
   
    B     = bspline.Bspline(k, p)
    B_array = numpy.array( [( B(t) ) for t in tau], dtype=numpy.float64 )

    B_array = numpy.delete(B_array, 0 , axis = 1)
    B_array = numpy.delete(B_array, -1 , axis = 1)
  
    P_nl = np.zeros(len(tau))   

    for position in range(len(tau)):
        for j in range(len(B_array[0,:])):
            P_nl[position]=((eigenev[j,n])*B_array[position,j])+ P_nl[position] 
    return P_nl,eigenvalues,H_matrix, B_matrix


#%%Test it
"""
Q=3
n = 5
l = n-1
summation = 0
P_nl = 0
electron = 0
V_ee = np.zeros((nknots))

#for iteration in range(iterations):
for shell in range(n):
    if electron ==Q:
        print("break at n:",shell+1)
        break
    print("n",shell)
        
    for subshell in range (l): 
        if electron ==Q:
            print("break at l:",subshell)
            break
        print("l",subshell)
    
        electron_max = 2*(2*subshell+1)
        P_nl = assignment5(Q,V_ee, subshell,shell,order_of_poly,R_max, e ,p, nknots,tau ) 
        N_occ = 0
        while electron<electron_max and electron<Q:
            electron = electron +1
            N_occ = N_occ+1
            
        print("N_orb",N_occ)
        Pnl = numpy.abs(P_nl)**2

        norm_P_nl = 0
        for position in range(len(tau)-1):
            norm_P_nl= trapezoidal(Pnl[position],Pnl[position+1],tau[position],tau[position+1]) + norm_P_nl
        Pnl = Pnl/norm_P_nl
        summation =e * N_occ * (Pnl)/(tau**2) + summation

rho = (1/(4*math.pi))*summation

F = 0
for position in range(len(tau)-1):
    F= trapezoidal(4*math.pi*rho[position]*tau[position]*tau[position],4*math.pi*rho[position+1]*tau[position+1]*tau[position+1],tau[position],tau[position+1]) + F
print("Control integral for Z =",Q,":",F)
"""


#%% CALCULATE POETNTIAL FOR ION
print("Before iterations")
iterations = 8  # iterate potential #iterations times 
Q =2
electronss = Q-1
V_old_ion = numpy.zeros((len(tau)))
for iteration in range(iterations): 
    tic = time.perf_counter()
    print("Iteration: ",iteration)
    n = 5   
    l = n-1
    eta = 0.4
    electron = 0
    E_total_ion = []
    V_ee_exch_ion = 0
    P_nl_ion = 0
    summation_ion = 0
    for shell in range(n):
        if electron ==electronss:
            break
        print("n",shell+1)
        electron_max_n = 2*((shell+1)**2)
        for subshell in range (l): 
            electron_max = 2*(2*subshell+1)
            if electron ==electron_max_n or electron ==electronss :
                break
            print("l",subshell) 

            P_nl_ion,E_nl_ion,H_matrix_ion,B_matrix_ion = assignment5(Q,V_old_ion, subshell,shell,order_of_poly,R_max, e ,p, nknots,tau) 
            E_s_ion= E_nl_ion[shell]
            occ_orbitals_ion = 0
            while occ_orbitals_ion<electron_max and electron<electronss :
                electron = electron + 1
                occ_orbitals_ion = occ_orbitals_ion+1
            ##ion
            Ps_ion= (P_nl_ion)**2
            norm_Ps_ion = 0
            for position in range(len(tau)-1):
                norm_Ps_ion= trapezoidal(Ps_ion[position],Ps_ion[position+1],tau[position],tau[position+1]) + norm_Ps_ion
            Ps_ion = Ps_ion/norm_Ps_ion
            integral_for_energy_ion = 0

            for position in range(len(tau)-1):
                integral_for_energy_ion = trapezoidal(Ps_ion[position]*V_old_ion[position], Ps_ion[position+1]*V_old_ion[position+1], tau[position],tau[position+1]) + integral_for_energy_ion
            E_total_ion.append( (occ_orbitals_ion)*(E_s_ion - 0.5 * integral_for_energy_ion) )
            summation_ion =e * (occ_orbitals_ion) * (Ps_ion)/(tau**2) +summation_ion
            print("E_orbital for ion : ",E_s_ion)
            toc = time.perf_counter()
            
    E_tott_ion= numpy.sum(E_total_ion)
    print("E_tot for ion",E_tott_ion)
    rho_ion=(1/(4*math.pi))*summation_ion
    k     = splinelab.augknt(knots, p)  # add endpoint repeats as appropriate for spline order p. THIS IS T
    B2     = bspline.Bspline(k, p) 
    A0 = B2.collmat(k)                 # collocation matrix for function value at sites tau. Rows are t and columns are B1,B2..
    A2 = B2.collmat(k, deriv_order=2)  # collocation matrix for second derivative at sites tau
    A1 = B2.collmat(k,deriv_order=1)
    orig_A=A2

    for i in range(2):
        A2=numpy.delete(A2,0 ,axis = 0)
    A2=numpy.delete(A2,-1 ,axis = 0)
    A2=numpy.delete(A2,-1 ,axis = 0)
    A2=numpy.delete(A2,1 ,axis = 0)
    A2=numpy.delete(A2,0 ,axis = 1)

    A2[-2,-3]=orig_A[3,2]
    A2[-2,-2]=orig_A[3,1]
    A2[-2,-1]=orig_A[3,0]
    A2[-1,-2]=A1[3,0]
    A2[-1,-1]=A1[3,1]
    rho_t_ion = numpy.zeros((nknots))
    counter = 0
    for i in range(len(tau)):
        for j in range(len(knots)):
            if abs(tau[i] - knots[j]) < 10e-4:
                counter += 1;
                rho_t_ion[j] = rho_ion[i]
    if counter == len(knots):
        print("SUCCESS!!") 
    
    rhs_ion = rho_t_ion * 4*math.pi * -knots
    rhs_ion = numpy.insert(rhs_ion, -1 , 0 , axis = 0)
    c_ion=numpy.linalg.solve(A2,rhs_ion)
    
    B_array2 = numpy.array(  [B2(t) for t in tau], dtype=numpy.float64 )
    B_array2 = numpy.delete(B_array2, 0,axis = 1)
    B_array2 = numpy.delete(B_array2, -1,axis = 1)
    column = numpy.zeros((len(rho_ion)))
    column[-1] = 1
    B_array2 =numpy.insert(B_array2,-1,column,axis = 1) 
    
    ####Update potatential for ion
    V_err_dir_ion = numpy.zeros((len(tau)))
    for position in range(len(tau)):
        for j in range(len(B_array2[0,:])):
            V_err_dir_ion[position] = ((c_ion[j])*B_array2[position,j])+ V_err_dir_ion[position] 
    V_err_dir_ion = (V_err_dir_ion)/tau
    
    V_ee_exch_ion = -3*((3*rho_ion)/(1*8*math.pi))**(1/3) 
    V_new_ion = e* (V_err_dir_ion + V_ee_exch_ion)
    V_new_ion = (1-eta)*V_new_ion + eta*V_old_ion 
    V_old_ion = V_new_ion
    print(f"One iteration in {toc - tic:0.4f} seconds")
############################################
#%% ########CALCULATE POTENTIAL FOR ATOM######
print("Before iterations")
iterations = 8  # iterate potential #iterations times 
Q =11
V_old = numpy.zeros((len(tau)))
for iteration in range(iterations): 
    tic = time.perf_counter()
    print("Iteration: ",iteration)
    n = 5   
    l = n-1
    eta = 0.4
    E_total = []
    V_ee_exch = 0
    P_nl = 0
    summation = 0
    electron = 0
    for shell in range(n):
        if electron ==Q:
            break
        print("n",shell+1)
        electron_max_n = 2*((shell+1)**2)
        for subshell in range (l): 
            electron_max = 2*(2*subshell+1)
            if electron ==electron_max_n or electron ==Q :
                break
            print("l",subshell) 
            occ_orbitals = 0
            P_nl,E_nl,H_matrix,B_matrix =  assignment5(Q,V_old, subshell,shell,order_of_poly,R_max, e ,p, nknots,tau) 
            E_s= E_nl[shell]

            while occ_orbitals<electron_max and electron<Q :
                electron = electron +1
                occ_orbitals = occ_orbitals+1
    
            Ps= (P_nl)**2
            norm_Ps = 0
            for position in range(len(tau)-1):
                norm_Ps= trapezoidal(Ps[position],Ps[position+1],tau[position],tau[position+1]) + norm_Ps
            Ps = Ps/norm_Ps
            integral_for_energy = 0
            print(occ_orbitals)
            for position in range(len(tau)-1):
                integral_for_energy = trapezoidal(Ps[position]*V_old[position], Ps[position+1]*V_old[position+1], tau[position],tau[position+1]) + integral_for_energy
            E_total.append( occ_orbitals*(E_s - 0.5 * integral_for_energy) )
            summation =e * occ_orbitals * (Ps)/(tau**2) + summation

            
            print("E_orbital for atom : ",E_s)
            toc = time.perf_counter()

    E_tott= numpy.sum(E_total)
    print("E_tot for atom",E_tott)

    
    rho= (1/(4*math.pi))*summation 
    #Create new B-splines for solving the 4:th assignment. 
    k     = splinelab.augknt(knots, p)  # add endpoint repeats as appropriate for spline order p. THIS IS T
    B2     = bspline.Bspline(k, p) 
    A0 = B2.collmat(k)                 # collocation matrix for function value at sites tau. Rows are t and columns are B1,B2..
    A2 = B2.collmat(k, deriv_order=2)  # collocation matrix for second derivative at sites tau
    A1 = B2.collmat(k,deriv_order=1)
    orig_A=A2

    for i in range(2):
        A2=numpy.delete(A2,0 ,axis = 0)
    A2=numpy.delete(A2,-1 ,axis = 0)
    A2=numpy.delete(A2,-1 ,axis = 0)
    A2=numpy.delete(A2,1 ,axis = 0)
    A2=numpy.delete(A2,0 ,axis = 1)
    A2[-2,-3]=orig_A[3,2]
    A2[-2,-2]=orig_A[3,1]
    A2[-2,-1]=orig_A[3,0]
    
    A2[-1,-2]=A1[3,0]
    A2[-1,-1]=A1[3,1]
    rho_t = numpy.zeros((nknots))
    counter = 0;
    for i in range(len(tau)):
        for j in range(len(knots)):
            if abs(tau[i] - knots[j]) < 10e-4:
                counter += 1;
                rho_t[j] = rho[i]
                #print("poop")
    if counter == len(knots):
        print("SUCCESS!!")            
    #Right hand side for both atom and ion. Add an actra rhs_last = 0 
    rhs = rho_t * 4*math.pi * -knots
    rhs = numpy.insert(rhs, -1 , 0 , axis = 0) 

    c=numpy.linalg.solve(A2,rhs)

    #B splines. Added one extra B_last = 1 
    B_array2 = numpy.array(  [B2(t) for t in tau], dtype=numpy.float64 )
    B_array2 = numpy.delete(B_array2, 0,axis = 1)
    B_array2 = numpy.delete(B_array2, -1,axis = 1)
    column = numpy.zeros((len(rho)))
    column[-1] = 1
    B_array2 =numpy.insert(B_array2,-1,column,axis = 1) 

    V_err_dir = numpy.zeros((len(tau)))
    ####Update potatential for atom
    for position in range(len(tau)):
        for j in range(len(B_array2[0,:])):
            V_err_dir[position] = ((c[j])*B_array2[position,j])+ V_err_dir[position] 
    V_err_dir = (V_err_dir)/tau

    V_ee_exch = -3*((3*rho)/(1*8*math.pi))**(1/3) 
    V_new = e* (V_err_dir + V_ee_exch)
    V_new = (1-eta)* V_new + eta*V_old
    V_old = V_new
    print(f"One iteration in {toc - tic:0.4f} seconds")

#%%
g = list()
f = open("energies.txt","a")
g.append(E_s)
g.append(E_s_ion)
g.append(E_tott)
g.append(E_tott_ion)
numpy.savetxt("energies.txt",numpy.array(g))
f.close()
#%%
plt.plot(tau,V_ee_exch)
plt.plot(tau,V_err_dir)
plt.plot(tau,V_old)
plt.grid("True")
plt.legend(["Exchange potential","Direct potential", "Total potential"])
plt.xlabel("r [$a_0$]")
plt.ylabel("Energy [$Hartree$]")
plt.savefig("ass6_potential_knot51_order4_Neon")
# %%
