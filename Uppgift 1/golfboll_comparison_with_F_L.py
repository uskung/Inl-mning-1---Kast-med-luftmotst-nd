import numpy as np
import matplotlib.pyplot as plt

plt.rc('font', size = 11, family='serif')
#plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern')

#################
# Define data for ball
m = 45.7*10**(-3) # mass [kg]
kL = 0.00031 # air resistance coefficient [kg/m]
kM = 0.000062 # Magnus coefficient [kg/rad]
################

###############
# Planet data
g = 9.82 # gravitational acceleration [m/s^2]
###############

############### 
# Initial conditions for ball throw
v0 = 75 # [m/s]

alfa0 = 12; ## in degrees
alfa0 = alfa0*np.pi/180 # converts to [rad]

w0 = 2000 ## angular velocity, in [rpm]
w0 = w0*2*np.pi/60 ## converts to [rad/s]

h = 0 # height at t=0 [m]
###############

###############
# Defines variables that are used in calculations
dt = 0.001 # step size for Euler's method [s]
vx = v0 * np.cos(alfa0)
vy = v0 * np.sin(alfa0)
alfa = alfa0 
w = w0
x = 0
y = h
xlist = [x]
ylist = [y]
###############

##### calculates if F_L = -k_l*v**2
while y >= 0: ## i.e. while ball is in the air
    # calculates new positions
    x = x + vx*dt
    y = y + vy*dt
    xlist.append(x)
    ylist.append(y)

    # calculates new angle, if-statement is required as range for -pi/2 < atan(x < pi/2
    if vx >= 0:
        alfa = np.atan(vy/vx)
    else:
        alfa = np.pi + np.atan(vy/vx)

    # updates v for this time step
    v = np.sqrt(vx**2 + vy**2)

    # calculates acceleration components for this time step
    drag = -kL/m * v**2
    magnus = kM/m * w *v
    ax = drag*np.cos(alfa) - magnus*np.sin(alfa)
    ay = drag*np.sin(alfa) + magnus*np.cos(alfa) - g

    # calculates velocity components for this time step
    vx = vx + ax*dt
    vy = vy + ay*dt

#############

#############
# calculates if F_L = -kL*v
kL = 0.0123
kM = 0.000052
vx = v0 * np.cos(alfa0)
vy = v0 * np.sin(alfa0)
alfa = alfa0
w = w0
x = 0
y = h
xlist2 = [x]
ylist2 = [y]

while y >= 0: ## i.e. while ball is in the air
    # calculates new positions
    x = x + vx*dt
    y = y + vy*dt
    xlist2.append(x)
    ylist2.append(y)

    # calculates new angle, if-statement is required as range for -pi/2 < atan(x < pi/2
    if vx >= 0:
        alfa = np.atan(vy/vx)
    else:
        alfa = np.pi + np.atan(vy/vx)

    # updates v for this time step
    v = np.sqrt(vx**2 + vy**2)

    # calculates acceleration components for this time step
    drag = -kL/m * v
    magnus = kM/m * w *v
    ax = drag*np.cos(alfa) - magnus*np.sin(alfa)
    ay = drag*np.sin(alfa) + magnus*np.cos(alfa) - g

    # calculates velocity components for this time step
    vx = vx + ax*dt
    vy = vy + ay*dt

axis = plt.gca() # plots current axes
axis.set_aspect('equal', adjustable='box')
plt.grid()
plt.plot(xlist,ylist,'.', color='blue', label='F_L = -kL*v^2')
plt.plot(xlist2,ylist2,'.', color='red', label='F_L = -kL*v')
plt.legend()
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.show()


