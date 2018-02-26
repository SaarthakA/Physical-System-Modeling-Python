#This program shows how to use the Euler-Cromer method for a pendulum

from numpy import *
import matplotlib.pyplot as plt
import numpy.fft as fourier

def accelerations(variables,params,t,i): #Computes derivatives of momenta or velocities
        g=params[0]
	theta = variables[i]
	
	vdot = -g*sin(theta)
	return vdot

def energy(coordinates,velocities,params,t):
    g=params[0]
    theta = coordinates[0]
    omega = velocities[0]
    
    KE = 0.5*omega**2
    U = g*(1-cos(theta))
    
    return (KE + U)
						
#Variables that we'll need
g=2*pi**2
#The 4*pi**2 is because we work in units where omega = 2pi, so the period of small oscillations is 1
omega0 = 0
degree = 49#angle in degrees
theta0= pi/180*degree #converting to radians 
dt = 0.04 #timestep
numtimes = 1*10**3 #totl number of calculations
params=array([g])
dof=1 #Number of coordinates, or "degrees of freedom" used to describe the system

#Establish some arrays
times = linspace(0, (numtimes-1)*dt, numtimes)
coordinates = zeros([numtimes,dof])
velocities = zeros([numtimes,dof])
acceleration = zeros([numtimes,dof])
energies=zeros(numtimes)


#initialize the arrays
i=0
coordinates[i] = theta0
velocities[i] = omega0
energies[i] = energy(coordinates[i],velocities[i], params,0)
acceleration[i]= accelerations(coordinates,params,times[i],i)

for i in range(1,numtimes):
    	coordinates[i] = coordinates[i-1]+dt*velocities[i-1]+0.5*acceleration[i-1]*dt**2
	acceleration[i]=accelerations(coordinates,params,times[i],i)
	velocities[i] = velocities[i-1]+(dt/2)*(acceleration[i-1]+acceleration[i])

        energies[i] = energy(coordinates[i],velocities[i], params,times[i])

spectrum = fourier.fft(coordinates[:,0]/numtimes) #We are normalizing the spectrum
spectrum = fourier.fftshift(spectrum) #Shift things so that the center frequency is zero
freq = linspace(-0.5/dt,0.5/dt,numtimes) #Create an array for the frequency
	
plt.figure()
plt.subplot(311)
plt.xlim(0,10)
plt.plot(times,coordinates[:,0],label=r'$\theta$')  #Make a plot
plt.xlabel("Time")
plt.ylabel("Angle")
plt.legend(loc = 'upper right')
plt.subplot(312)
plt.plot(times,energies,label='Velocity Verlet Energy')
plt.xlim(0,10)
plt.xlabel("Time")
plt.ylabel("Energy")
plt.legend(loc = 'upper right')
plt.subplot(313)
plt.plot(freq, abs(spectrum)**2,label='Amplitude Spectrum')
plt.xlim(-3.1,3.1)
plt.xlabel("Frequency")
plt.ylabel("Power")
plt.legend(loc = 'upper right')
plt.show()
