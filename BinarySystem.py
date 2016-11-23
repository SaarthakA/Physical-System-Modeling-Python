from numpy import *
import matplotlib.pyplot as plt

def acceleration(variables,params,t):
	
	xs1=variables[0]
	ys1=variables[1]
	
	xs2=variables[2]
	ys2=variables[3]
	
	MS1G = params[0] #m-star1
	MS2G = params[1] #m-star2
	
	rStars = sqrt((xs1-xs2)**2+(ys1-ys2)**2)#distance between both stars

	as1 = array([xs2-xs1,ys2-ys1])*MS2G/(rStars**3) #acceleration of star 1
        as2 = array([xs1-xs2,ys1-ys2])*MS1G/(rStars**3) #acceleration of star 2
        	
        
        return array([as1[0],as1[1],as2[0],as2[1]])

#Star 1 conditions   
xs1 = pi
ys1 = 0
vs1x0 = 0
vs1y0 = 0.2*2*pi/xs1**0.5

#Star 2 conditions
xs20 = -xs1
ys20 = -ys1
vs2x0 = -vs1x0
vs2y0 = -vs1y0

MS1G = 4*pi**2 #mass star 1
MS2G = MS1G #mass star 2
params = array([MS1G, MS2G],dtype='float')

dt = 0.001
onehalfdtsquared = 0.5*dt**2
tmax = 210
numtimes = int(tmax/dt)


times = linspace(0,tmax,numtimes)

coordinates = zeros([numtimes,4])
velocities = zeros([numtimes,4])

coordinates[0] = array([xs1,ys1,xs20,ys20])
velocities[0]  = array([vs1x0,vs1y0,vs2x0,vs2y0])

aold = acceleration(coordinates[0],params,0)

#Velocity Verlet Calculations
for t in range(1,numtimes):
    coordinates[t] = coordinates[t-1]+velocities[t-1]*dt+onehalfdtsquared*aold
    anew = acceleration(coordinates[t],params,times[t-1])
    velocities[t] = velocities[t-1]+0.5*(anew+aold)*dt
    aold = anew

r = zeros(numtimes)

for t in range(0,numtimes):
    r[t] = sqrt(coordinates[t,0]**2+coordinates[t,1]**2)


plt.figure()
plt.subplot(311)
plt.plot(coordinates[:,0],coordinates[:,1], label = 'Star 2')
plt.plot(coordinates[:,2],coordinates[:,3], label = 'Star 1')
plt.legend(loc = 'upper right')
plt.axis('equal')
plt.subplot(312)
plt.plot(times, r)
plt.xlabel('Time')
plt.ylabel('Distance from each other, r')
plt.show()