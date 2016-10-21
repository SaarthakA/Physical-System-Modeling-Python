#This program shows how to use the VelocityVerlet method for a pendulum

from numpy import *
import matplotlib.pyplot as plt
import numpy.fft as fourier

def accelerations(variables,params,t): #Computes derivatives of momenta or velocities
        g=params[0]
	x = variables[0]
	
	vdot = -g*(x-x**2+(x**5)/100)

	return vdot

def energy(coordinates,velocities,params,t):
    g=params[0]
    theta = coordinates[0]
    v = velocities[0]
    
    KE = 0.5*v**2
    U = g*(0.5*theta**2-(1/3)*theta**3+(theta**6)/600)                #mass is equal to 1
    
    return KE+U
						
#Variables that we'll need
g=4*pi**2
#The 4*pi**2 is because we work in units where omega = 2pi, so the period of small oscillations is 1
v0 = 0
theta0 = 0.0927             #0.04349 gives an average of .0010003, 0.0927 gives 0.004002
dt = 0.04
numtimes = 3*10**2
params=array([g])
onehalfdtsq = .5*dt**2
onehalfdt=.5*dt
dof=1 #Number of coordinates, or "degrees of freedom" used to describe the system

#Establish some arrays
times = linspace(0, (numtimes-1)*dt, numtimes)
coordinates = zeros([numtimes,dof])
velocities = zeros([numtimes,dof])
energies=zeros(numtimes)


#initialize the arrays
i=0
coordinates[i] = theta0
velocities[i] = v0
energies[i] = energy(coordinates[i],velocities[i], params,0)

aold = accelerations(coordinates[0],params,times[0])

for i in range(1,numtimes):
        coordinates[i] = coordinates[i-1]+dt*velocities[i-1]+onehalfdtsq*aold
	anew = accelerations(coordinates[i],params,times[i])
	velocities[i] = velocities[i-1]+onehalfdt*(aold+anew)
        energies[i] = energy(coordinates[i],velocities[i], params,times[i])
        aold = anew
avgx = average(coordinates,axis=0)
print (avgx)

spectrum = fourier.fft(coordinates[:,0]/numtimes) #We are normalizing the spectrum
spectrum = fourier.fftshift(spectrum) #Shift things so that the center frequency is zero
freq = linspace(-0.5/dt,0.5/dt,numtimes) #Create an array for the frequency
	
plt.figure()
plt.subplot(311)
plt.plot(times[:26],coordinates[:26],label=r'$\theta$')  #Make a plot
plt.xlabel("Time")
plt.ylabel("Angle")
plt.legend(loc = 'upper right')
plt.subplot(312)
plt.plot(times,energies,label='Energy (Velocity-Verlet)')
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