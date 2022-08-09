#%%
import numpy
%matplotlib inline
import bspline
import bspline.splinelab as splinelab
import matplotlib.pyplot as plt
from numba import njit
from __future__ import division
from pylab import *
from scipy.special.orthogonal import p_roots
from scipy import integrate
R_max = 35
r = 0.00
Q = 1
l = 0

p = 5-1             # order of spline i.e polynomial order 3
## Spline setup and evaluation
nknots = 35    # number of knots to generate (here endpoints count only once) endpoints = 4+4
#tau = numpy.linspace(r,R_max,nknots)  # collocation sites (i.e. where to evaluate)
knots = numpy.linspace(r,R_max,nknots)  # create a knot vector without endpoint repeats
k     = splinelab.augknt(knots, p)  # add endpoint repeats as appropriate for spline order p. THIS IS T
B     = bspline.Bspline(k, p)       # create spline basis of order p on knots k
order_of_poly = 8

#%%


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
        F = 0
        for position in range(len(knots)-1):

            F=gauss(B,order_of_poly,knots[position],knots[position+1],i+1,j+1) + F
            B_matrix[i,j] = F
#%%

def polynom_function(B,r,w,a,b,i,j,l,Q):
    term_1 = 0.5*(B.collmat(r*(b-a)*0.5 + (a+b)*0.5,deriv_order=1)[i]*B.collmat(r*(b-a)*0.5 + (a+b)*0.5,deriv_order=1)[j])
    term_2 = l*(l+1)/((r*(b-a)*0.5 + 0.5*(a+b))**2)*(B.collmat(r*(b-a)*0.5 + (a+b)*0.5)[i]*B.collmat(r*(b-a)*0.5 + (a+b)*0.5)[j])
    term_3 = -(Q/(r*(b-a)*0.5 + 0.5*(a+b)))* ((B.collmat(r*(b-a)*0.5 + (a+b)*0.5)[i]) * (B.collmat(r*(b-a)*0.5 + (a+b)*0.5)[j]))
    return w*(term_1+term_2+term_3)

def gauss_H(B,n,a,b,i,j,l,Q):
    [r,w] = p_roots(n+1)

    poly = 0
    for order in range(len(r)):
        
        poly=polynom_function(B,r[order],w[order],a,b,i,j,l,Q) + poly
    G=(b-a)*0.5*poly
    return G

H_matrix = numpy.zeros((nknots-2 + p-1, nknots-2 + p-1))
for i in range(len(H_matrix)):
    for j in range(len(H_matrix)):
        print(i,j)
        F = 0
        for position in range(len(knots)-1):
            F= gauss_H(B,order_of_poly,knots[position],knots[position+1],i+1,j+1,l,Q) + F
            H_matrix[i,j] = F

#%%
from scipy.linalg import eigh
eigenvalues, eigenev= eigh(H_matrix,B_matrix)
eigenev
#c does not depend on R
B     = bspline.Bspline(k, p)   
B_array = numpy.array( [( B(t) ) for t in knots], dtype=numpy.float64 )
B_array = numpy.delete(B_array, 0 , axis = 1)
B_array = numpy.delete(B_array, -1 , axis = 1)

psi = np.zeros((B_array.shape))
radial = np.zeros((len(B_array)))

for i in range(len(B_array[:,0])):
    for j in range(len(B_array[0,:])):
    #eigenev[i]=eigenev[i]/numpy.linalg.norm(eigenev[i])
    #B_array[i,:]=B_array[i,:]/numpy.linalg.norm(B_array[i,:])
        psi[i,j]=(numpy.sum((eigenev[i,j])*B_array[i,:]))
        radial[i] = B_array[i,:]*B_array[i,:]
        #radial.append(numpy.sum(numpy.abs((psi[i])**2)))
plt.plot(knots,psi[:,0])
plt.plot(knots,psi[:,1])
plt.plot(knots,psi[:,2])
plt.plot(knots,psi[:,3])

plt.grid("True")
plt.ylabel("Energy [Ha]")
plt.xlabel("r [$a_0$]")
plt.legend(["$\Psi(r)_1$","$\Psi(r)_2$","$\Psi(r)_3$","$\Psi(r)_4$"])

plt.savefig("35_0_5_7_35_8")

#plt.plot(knots,radial)
# %%
