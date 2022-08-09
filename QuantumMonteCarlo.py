#%%
from cmath import exp, sqrt
import black
import numpy as np
import matplotlib.pyplot as plt
import random
from numba import njit
import math 

N = 1
dim = 1
c = 1


@njit
def wave_function(x,alpha):
    wave_trial= math.e**(-alpha*x*x)
    return wave_trial
@njit
def random_displacement(old_position):
    random_move= np.random.uniform(-0.2,0.2,1)
    new_position = old_position + random_move
    return new_position 
@njit
def energy_local(alpha,x):
    e = alpha + x*x*((1/2 - 2*alpha*alpha))
    return e 
@njit
def var(alpha,x):
    e_x_2 = (alpha + x*x*((1/2 - 2*alpha*alpha)))*(alpha + x*x*((1/2 - 2*alpha*alpha)))
    return e_x_2


#%%
position = np.random.uniform(-1,1,1)
alphas = np.linspace(0.1, 1 ,37) #analytical = 1/2
steps = 6000
step_sufficient = 2000
run = 100
variance_variance = []
average_average =[]

points_visited= []
wave_prob = []
points_visited = np.zeros((steps+1))
for alpha in alphas: 
    print(alpha)
    average_energy =np.zeros((run))
    variance = np.zeros((run))
    average_2  = np.zeros(((run)))
    for loops in range(run): # We want to compute this for the same alpha some 100 times 
        e_l = []
        variancee = []
        for step in range(steps): #do several steps and then add local energy
            wave_old =wave_function(position,alpha)
            new_position = random_displacement(position)
            wave_new= wave_function(new_position,alpha)
            probability = np.abs(wave_new)*np.abs(wave_new) /(np.abs(wave_old)*np.abs(wave_old) )
            if probability >= 1:
                position = new_position
                points_visited[step] = (position)
            else:
                tolerence = np.random.uniform(0,1)
                if tolerence<probability:
                    position= new_position  
                    points_visited[step] = (position)
                else:
                    position = position
                    points_visited[step] = (position)
            if step>=step_sufficient:
                e_l.append(energy_local(alpha,position))
                variancee.append(var(alpha,position))

            acc_rate=len(e_l)/(steps)
        average_energy[loops] = (np.sum(e_l)/len(e_l))
        average_2[loops] = ((np.sum(e_l)/len(e_l))*(np.sum(e_l)/len(e_l)))
        variance[loops]=((np.sum(variancee)/len(variancee))-(np.sum(e_l)/len(e_l))*(np.sum(e_l)/len(e_l)))
    average_average.append( np.sum(average_energy)/len(average_energy))
    variance_variance.append((np.sum(variance)/len(variance)))
#%%  analytical wave function
analytical_wave = []
alpha = 1/2
x = np.linspace(-2, 2 ,100 ) 
dx =x[1]-x[0]
wave = 0
for steps in x:
    analytical_wave.append((wave_function(steps, alpha))*(wave_function(steps, alpha)))
wave = np.sum(analytical_wave)*dx+ wave

#n=plt.plot(points_visited,analytical_wave)
above_lambda1 = []
below_lambda1 = []
for i in range(len(variance_variance)):
    below_lambda1.append(average_average[i]-variance_variance[i])
    above_lambda1.append(average_average[i]+variance_variance[i])

plt.fill_between(alphas, below_lambda1, above_lambda1,
                 color='gray', alpha=0.2)
plt.plot(alphas,average_average,"00000")
plt.grid("True")
plt.xlabel("Alpha")
plt.yticks(np.arange(-1.5,4,step = 1))
plt.ylabel("Average of <E_L>")
plt.savefig("Energytoalpha_variance_1D")
plt.show()

n, bins, patches=plt.hist(x=points_visited,density = True,stacked = True ,bins = 200,
                            alpha=0.8, rwidth=0.7,color = "#00C957")
plt.plot(x,analytical_wave/wave,color= "00000")
plt.grid(axis='y', alpha=0.75)
plt.ylabel('Probability Density')
plt.xlabel('x')
plt.legend(["Trial wave","Points visited"])
plt.savefig("Probability_hist")


#%%
plt.errorbar(alphas,average_average,yerr= variance_variance)
plt.grid("True")
plt.xlabel("Alpha")
plt.ylabel("Average of <E_L>")
plt.savefig("Energytoalphavariance_1D")
plt.show()
plt.plot(alphas, average_average)
plt.grid("True")
plt.xlabel("Alpha")
plt.ylabel("Average of <E_L>")
plt.savefig("Energytoalpha_1D")



#%% 2 PARTICLES 2D 
N = 2
dim = 2
c = 1
x = np.zeros((N))
y = np.zeros((N))


e_l = []
@njit
def wave_function(x ,y ,alpha  ,r_12):
    wave_trial= math.e**(-0.5*(x[0]**2 + y[0]**2 + x[1]**2 + y[1]**2)) * math.e**((lambdaa*r_12)/(1+ alpha*r_12))
    return wave_trial

@njit
def random_displacement(old_position):
    random_move= np.random.uniform(-0.5,0.5,(N))
    new_position = old_position + random_move
    return new_position 
@njit
def energy_local(alpha, lambdaa, r_12):
    e = 0.5 *(1-(lambdaa/(r_12*(1+alpha*r_12)**2)))*((r_12 * 2* lambdaa/(1+alpha*r_12)**2) + 4 ) + (lambdaa/r_12)*(1+ (1+ 3*alpha*r_12)/(1+alpha*r_12)**3)
    return e 
@njit
def square(alpha,lambdaa,r_12):
    e_x_2 = (0.5 *(1-(lambdaa/(r_12*(1+alpha*r_12)**2)))*((r_12 * 2* lambdaa/(1+alpha*r_12)**2) + 4 ) + (lambdaa/r_12)*(1+ (1+ 3*alpha*r_12)/(1+alpha*r_12)**3))*(0.5 *(1-(lambdaa/(r_12*(1+alpha*r_12)**2)))*((r_12 * 2* lambdaa/(1+alpha*r_12)**2) + 4 ) + (lambdaa/r_12)*(1+ (1+ 3*alpha*r_12)/(1+alpha*r_12)**3))
    return e_x_2

x = np.random.uniform(-1,1,(N))
y = np.random.uniform(-1,1,(N))


steps = 3000
step_sufficient = 1000
run = 100
e_squared = []
variance_variance_1 = []
average_average_1 =[]
alphas = np.linspace(0.1, 1 ,28) 
lambdaa = 2
for alpha in alphas:
    print(alpha)
    average_energy_1 =np.zeros((run))
    variance= np.zeros((run))
    average_2  = np.zeros(((run)))
    acc_rate=np.zeros((run))
    for loops in range(run):
        e_l = []
        variancee = []
        count = 0

        for step in range(steps):
            r_12_old = np.sqrt((x[0]- x[1])**2 + (y[0] - y[1])**2)
            wave_old = wave_function(x,y,alpha, r_12_old)
            x_new_position = random_displacement(x)
            y_new_position = random_displacement(y)
            r_12_new = np.sqrt((x_new_position[0]- x_new_position[1])**2 + (y_new_position[0] - y_new_position[1])**2)
            wave_new= wave_function(x_new_position,y_new_position,alpha, r_12_new)
            probability = np.abs(wave_new)*np.abs(wave_new) /(np.abs(wave_old)*np.abs(wave_old) )
            if probability >= 1:
                    x = x_new_position    
                    y = y_new_position   
                    count = count+ 1
                    if step>=step_sufficient:
                        e_l.append(energy_local(alpha,lambdaa,r_12_new))
                        e_squared.append(square(alpha,lambdaa,r_12_new))

            else:
                tolerence = np.random.uniform(0,1)
                if tolerence < probability:
                    x = x_new_position    
                    y = y_new_position   
                    count = count + 1  
                    if step>=step_sufficient:
                        e_l.append(energy_local(alpha,lambdaa,r_12_new))
                        e_squared.append(square(alpha,lambdaa,r_12_new))

                else:
                    x = x
                    y = y
                    if step>=step_sufficient:
                        e_l.append(energy_local(alpha,lambdaa,r_12_old))
                        e_squared.append(square(alpha,lambdaa,r_12_old))
    
            acc_rate[loops] = count/(steps)
    
        average_energy[loops] = (np.sum(e_l)/len(e_l))
        average_2[loops] = (np.sum(e_l)/len(e_l)*np.sum(e_l)/len(e_l))
        variance[loops]=((np.sum(e_squared)/len(e_squared))-np.sum(e_l)/len(e_l)*np.sum(e_l)/len(e_l))
    average_average_1.append( np.sum(average_energy)/len(average_energy))
    variance_variance_1.append((np.sum(variance)/len(variance)))
= np.loadtxt("variance_variance_11.txt")

#%%
#Golden Search

a = (average_average_0[0])
position_a = average_average_0.index((average_average_0[0]))
b = (average_average_0[-1])
position_b = average_average_0.index((average_average_0[-1]))
search = 8
golden_ratio = (math.sqrt(5)- 1)/2

golden1 = np.zeros((search))
pos1 = np.zeros((search))

golden1[0] = a
golden1[-1] = b
pos1[0] = alphas[0]
pos1[-1] = alphas[-1]

k = 2
j = 1

for i in range(1,search-1):
    position_d = round(golden_ratio*(position_b-position_a))
    x_1=position_a+position_d
    x_2 = position_b-position_d
    f_x1 = average_average_1[x_1]
    f_x2 = average_average_1[x_2]
    

    
    if f_x1<f_x2:
        golden1[search-k] = f_x2
        pos1[search-k] = alphas[x_2]

        position_a = x_2
        a = average_average_1[x_2]
        x_2 = x_1
        f_x2=f_x1
        
        k = k+1
        
    if f_x1>f_x2:
        golden1[j] = f_x1
        pos1[j] = alphas[x_1]
        
        position_b =x_1
        b = average_average_1[x_1]
        
        x_1 = x_2
        f_x1= f_x2 

        j = j+1
    if f_x1==f_x2:
        golden1[k] = f_x2
        pos1[k]=alphas[x_2]
        print(pos1[k],golden1[k])


#Run for lambdas seperately and save variables 
plt.scatter(pos[k],golden[k],color = "r", s = 50, marker="X")
plt.axhline(golden[k], color = "r",linestyle= "--",alpha=0.6)
plt.yticks(np.arange(golden[k],max(average_average_0),1))

plt.grid("True")
plt.xlabel("Alpha")
plt.ylabel("$E / \hbar\omega$")
plt.legend(["Average <E_L>","Minimum energy"])
plt.show()

#variance 
above_lambda1 = []
below_lambda1 = []

for i in range(len(variance_variance_0)):
    below_lambda1.append(average_average_0[i]-variance_variance_0[i])
    above_lambda1.append(average_average_0[i]+variance_variance_0[i])

    
plt.plot(alphas,average_average_0,color = "r")
plt.fill_between(alphas, below_lambda1, above_lambda1,
                 color='gray', alpha=0.2)

plt.grid("True")
plt.legend(["average <E_L> with $\lambda = 1 $","average <E_L> with $\lambda = 2 $"])
plt.xlabel("Alpha")
plt.ylabel("Average of <E_L>")
plt.savefig("Energytoalpha_variance_2D")
#plt.plot(pos[0:k-1],golden[0:k-1],"^",color = "gray")
#plt.plot(pos[k-1],golden[k-1],"^",color = "gray")
#plt.plot(pos[k+1:search],golden[k+1:search],"^",color = "gray")
plt.savefig("lambda_8")

# %%
