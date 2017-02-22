from numpy import *
import matplotlib.pyplot as plt

def acceleration(variables,params,t):
	
	xs1=variables[0]
	ys1=variables[1]
	
	xs2=variables[2]
	ys2=variables[3]
	
	xp1=variables[4]
	yp1=variables[5]
	
	MS1G = params[0] #m-star1
	MS2G = params[1] #m-star2
	
	rStars = sqrt((xs1-xs2)**2+(ys1-ys2)**2)#distance between both stars
        rPlanet1 = sqrt((xs1-xp1)**2+(ys1-yp1)**2)
        rPlanet2 = sqrt((xs2-xp1)**2+(ys2-yp1)**2)
        
	as1 = array([xs2-xs1,ys2-ys1])*MS2G/(rStars**3) #acceleration of jupiter wrt sun
        as2 = array([xs1-xs2,ys1-ys2])*MS1G/(rStars**3) #acceleration of jupiter wrt sun
        ap1 = array([((xs1-xp1)*(MS1G/(rPlanet1**3)))+((xs2-xp1)*(MS2G/(rPlanet2**3))),((ys1-yp1)*(MS1G/(rPlanet1**3)))+((ys2-yp1)*(MS2G/(rPlanet2**3)))])
   
        
        return array([as1[0],as1[1],as2[0],as2[1],ap1[0],ap1[1]])

#Mass variables
MS1G = 4*pi**2
MS2G = MS1G/2
MP1G = MS1G/1350

# Star 1 Conditions    
xs1 = pi
ys1 = 0
vs1x0 = 0
vs1y0 = 0.2*2*pi/xs1**0.5

# Star 2 Conditions
xs20 = -xs1
ys20 = -ys1
vs2x0 = -vs1x0 * MS1G/MS2G
vs2y0 = -vs1y0 * MS1G/MS2G

#Planet conditions
xp1 = pi+6 #6 is good orbit, 5 ok, <=4 bad
yp1 = 0
vp1x0 = 0
vp1y0 = (2*pi/xp1**0.5)-5.5 #5.5 good speed

params = array([MS1G, MS2G, MP1G],dtype='float')

dt = 0.001
onehalfdtsquared = 0.5*dt**2
tmax = 150
numtimes = int(tmax/dt)


times = linspace(0,tmax,numtimes)

coordinates = zeros([numtimes,6])
velocities = zeros([numtimes,6])

coordinates[0] = array([xs1,ys1,xs20,ys20,xp1, yp1])
velocities[0]  = array([vs1x0,vs1y0,vs2x0,vs2y0,vp1x0,vp1y0])

aold = acceleration(coordinates[0],params,0)

for t in range(1,numtimes):
    coordinates[t] = coordinates[t-1]+velocities[t-1]*dt+onehalfdtsquared*aold
    anew = acceleration(coordinates[t],params,times[t-1])
    velocities[t] = velocities[t-1]+0.5*(anew+aold)*dt
    aold = anew

r1 = zeros(numtimes)
r2 = zeros(numtimes)
for t in range(0,numtimes):
    r1[t] = sqrt(coordinates[t,0]**2+coordinates[t,1]**2)
    r2[t] = sqrt(coordinates[t,4]**2+coordinates[t,5]**2)

plt.figure()
plt.subplot(311)
plt.plot(coordinates[:,0],coordinates[:,1], label = 'Star 2')
plt.plot(coordinates[:,2],coordinates[:,3], label = 'Star 1')
plt.plot(coordinates[:,4],coordinates[:,5], label = 'Planet 1')
plt.legend(loc = 'upper right')
plt.axis('equal')
plt.subplot(312)
plt.plot(times, r1)
plt.xlabel('Time')
plt.ylabel('Distance from star 1 to star 2, r')
plt.subplot(313)
plt.plot(times, r2)
plt.xlabel('Time')
plt.ylabel('Distance of planet from CM, r')
plt.show()