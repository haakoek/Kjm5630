import numpy as np
import matplotlib.pyplot as plt
import sys

def exact_solution(t,x0,xe,v0,w):
	
	A = x0-xe
	B = v0/w

	x_exact = A*np.cos(w*t)+B*np.sin(w*t)+xe
	v_exact = w*(B*np.cos(w*t)-A*np.sin(w*t))

	return x_exact, v_exact

def F(x,xe,k):
	return -k*(x-xe)

#physical characteristics
m  = 1.0 #mass of particle
k  = 1.0 #spring constant
xe = 0   #equlibrium position
w  = np.sqrt(k/m)



#initial conditions
x0 = 1
v0 = 1
F0 = -k*(x0-xe)

x_list = [x0]
v_list = [v0]
F_list = [F0]

#numerical parameters
T     = 6.5   			#Total simulation time
dt    = 0.1				#Time step
nst   = int(T/dt)	    #Number of time steps
t_arr = np.linspace(0,T,nst+1) #This is spooky
x_exact, v_exact = exact_solution(t_arr,x0,xe,v0,w)

#Time integration 
#Verlet
#Compute special case, x1
x1 = x0 + v0*dt + 0.5*(F0/m)*dt**2

for i in range(1,nst+1):
	
	t1 = i*dt
	F1 = F(x1,xe,k)
	x2 = 2*x1 - x0 + (F1/m)*dt**2
	v1 = 0.5*(x2-x0)/dt

	x_list.append(x1)
	v_list.append(v1)
	F_list.append(F1)

	#update positions
	x0 = x1
	x1 = x2

x_list = np.array(x_list)
v_list = np.array(v_list)

plt.figure(1)
plt.plot(t_arr,x_list,'r')
plt.plot(t_arr,x_exact,'b')
plt.xlabel('t')
plt.ylabel('x(t)')

plt.figure(2)
plt.plot(t_arr,v_list,'or')
plt.plot(t_arr,v_exact,'b')
plt.xlabel('t')
plt.ylabel('v(t)')

plt.show()
