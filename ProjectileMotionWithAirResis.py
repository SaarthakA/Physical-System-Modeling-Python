#This program is a template for you to fill in 
from numpy import *
import matplotlib.pyplot as plt

def derivatives(variables,params,t): #We are defining a function that returns derivatives
    #Variables is an array holding our variables
    x=variables[0]
    y=variables[1]
    vx=variables[2]
    vy=variables[3]
    
    #params is an array holding parameters
    m=params[0]
    g = params[1]
    v0 = params[2]
    C = params[3]
    dt = params[4]
    vT = params[5]
    
    dragx=(C/m)*v0*vx
    dragy=(C/m)*v0*vy
    
    #Now we compute the derivatives
    dxdt=vx
    dydt=vy
    dvxdt=-vx-(dragx/m)
    dvydt=-vy-g-(dragy/m)
        
    return array([dxdt, dydt, dvxdt, dvydt])
	
#Initialize some stuff
dt=0.001 #Time step
numtimes = 10000 #How many times the loop will repeat
mass = 1.0
g = 9.8
v0=10 #initial velocity
vT = 1
C = 0.1 #air res constant
x0=0 #initial x pos
y0=0.0 #initial y pos
degree = 40 #angle in degrees
theta0= pi/180*degree #converting to radians 
vx0=v0*cos(theta0)
vy0=v0*sin(theta0)
#Package our parameters
params = array([mass, g, v0, C, dt,vT])

#Establish some arrays
times=linspace(0,(numtimes-1)*dt,numtimes)
#0 is the starting element in the array
#(numtimes-1)*dt is the last element
#numtimes is the total number of elements
variables = zeros([numtimes,4])
#initialize the arrays


variables[0] = array([x0,y0,vx0,vy0])

for t in range(1,numtimes):
    df = derivatives(variables[t-1], params, times[t-1])
    variables[t] = variables[t-1]+dt*df
    if (variables[t,1]<0):
        variables[t,1] = variables[t-1,1]
    if (variables[t,0]<0):
        variables[t,0] = variables[t-1,0]
    		
plt.figure() #Creates a figure
plt.subplot(411)
plt.plot(times,variables[:,0],label='x-position')
plt.legend(loc='upper right')
plt.subplot(412)
plt.plot(times,variables[:,1],label='y-position')
plt.legend(loc='upper right')
plt.subplot(413)
plt.plot(times,variables[:,2],label='Vx')
plt.legend(loc='upper right')
plt.subplot(414)
plt.plot(times,variables[:,3],label='Vy')#Adds graphs to it.  The first variable is on the horizontal axis, the second is on the vertical axis.
plt.legend(loc='upper right')
plt.show() #Shows it