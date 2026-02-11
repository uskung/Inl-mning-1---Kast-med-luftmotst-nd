import numpy as np
import matplotlib.pyplot as plt

#################
# Define data for ball
m = 2.7*10**(-3) # mass of ping pong ball [kg]
kL =  18 * 0.5 * 0.2 * 1.29 * 1.26*10**(-3)   # air resistance coefficient [kg/m] 0.00031
kM = 0.18 * 0.5 * 1.23 * 1.29 * (np.pi * 0.02**2) * 0.02  # Magnus coefficient [kg/rad] 0.000062
################
###############
# Planet data
g = 9.82 # gravitational acceleration [m/s^2]
###############

############### 
# Initial conditions for ball throw
v0 = 20 # [m/s]

alfa0 = 8; ## in degrees
alfa0 = alfa0*np.pi/180 # converts to [rad]

#w0 = 0 ## angular velocity, in [rpm]
#w0 = w0*2*np.pi/60 ## converts to [rad/s]
w0 = 754 # [rad/s]. If negative => topspin, if positiv => backspin
h = 0.20 # height at t=0 [m]
###############

###############
# Defines variables that are used in calculations
dt = 0.0001 # step size for Euler's method [s]
vx = v0 * np.cos(alfa0)
vy = v0 * np.sin(alfa0)
alfa = alfa0 
w = w0
x = 0
y = h
xlist = [x]
ylist = [y]
###############

time = 0
time_constant = False

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

    if (time_constant == False) and (time >= 1): ## i.e time is 1 s
        print(f"Farten v vid t=1 s är: {v}")
        time_constant = True

    time += dt
#############

print(f"Hela slaget tog: {time:.3f} sekunder.")

plt.rc('font', size = 11, family='serif')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern')
axis = plt.gca() # plots current axes

axis.set_aspect('equal', adjustable='box')
if xlist[-1] < 2.74: ## if ball doesn't go over table, plot the length of table
    plt.xlim(0, 2.74)
else: # ie ball goes over table
    plt.plot([2.74, 2.74], [0,0.10], color='red', label='bordskant')
    plt.xlim(0,xlist[-1])

plt.plot(xlist,ylist,'.', color='blue', label='$F_L = -kL*v^2$')
#plt.legend()
plt.plot([1.37,1.37],[0,0.1525], color='black', label='nät') ## net
plt.grid()
plt.xlabel('$x$ [m]')
plt.ylabel('$y$ [m]')
plt.show()


