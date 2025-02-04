import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.integrate import solve_ivp

# -----------------------------
# 1. Define system parameters
# -----------------------------
m1 = 1.0       # Mass of first pendulum bob
m2 = 1.0       # Mass of second pendulum bob
L1 = 1.0       # Length of first pendulum rod
L2 = 1.0       # Length of second pendulum rod
g  = 9.81      # Acceleration due to gravity

# Initial angles (in radians) and angular velocities
theta1_0 = np.pi/2  # 90 degrees
theta2_0 = np.pi/2  # 90 degrees
omega1_0 = 0.0
omega2_0 = 0.0

# Time span for the simulation
t_start = 0.0
t_end   = 20.0
dt      = 0.01  # time step for saving solution

# -----------------------------
# 2. Define the ODE function
# -----------------------------
def double_pendulum_ode(t, y, m1, m2, L1, L2, g):
    """
    Returns the derivatives [theta1_dot, omega1_dot, theta2_dot, omega2_dot]
    for the double pendulum system.
    """
    theta1, omega1, theta2, omega2 = y

    # Common terms to simplify equations
    cos12 = np.cos(theta1 - theta2)
    sin12 = np.sin(theta1 - theta2)
    denom = 2*m1 + m2 - m2 * np.cos(2*theta1 - 2*theta2)

    dtheta1 = omega1
    dtheta2 = omega2

    # omega1_dot
    num1 = -g*(2*m1 + m2)*np.sin(theta1)
    num2 = -m2*g*np.sin(theta1 - 2*theta2)
    num3 = -2*sin12*m2*(omega2**2*L2 + omega1**2*L1*cos12)
    domega1 = (num1 + num2 + num3) / (L1 * denom)

    # omega2_dot
    num4 = 2*sin12*(omega1**2*L1*(m1 + m2) + g*(m1 + m2)*np.cos(theta1) + omega2**2*L2*m2*cos12)
    domega2 = num4 / (L2 * denom)

    return [dtheta1, domega1, dtheta2, domega2]


# -----------------------------
# 3. Numerical integration
# -----------------------------
# Initial condition vector
y0 = [theta1_0, omega1_0, theta2_0, omega2_0]

# Time points for dense output
t_eval = np.arange(t_start, t_end, dt)

# Integrate using solve_ivp (RK45 by default)
sol = solve_ivp(
    fun=double_pendulum_ode,
    t_span=(t_start, t_end),
    y0=y0,
    t_eval=t_eval,
    args=(m1, m2, L1, L2, g)
)

# Extract solutions
theta1_vals = sol.y[0]
omega1_vals = sol.y[1]
theta2_vals = sol.y[2]
omega2_vals = sol.y[3]
time_vals   = sol.t

# -----------------------------
# 4. Convert to Cartesian coordinates for plotting
# -----------------------------
x1 = L1 * np.sin(theta1_vals)
y1 = -L1 * np.cos(theta1_vals)

x2 = x1 + L2 * np.sin(theta2_vals)
y2 = y1 - L2 * np.cos(theta2_vals)

# -----------------------------
# 5. Create the animation
# -----------------------------
fig, ax = plt.subplots(figsize=(5, 5))
ax.set_xlim(- (L1 + L2) * 1.2, (L1 + L2) * 1.2)
ax.set_ylim(- (L1 + L2) * 1.2, (L1 + L2) * 1.2)
ax.set_aspect('equal', 'box')
ax.set_title("Double Pendulum")

line, = ax.plot([], [], 'o-', lw=2, markersize=8)  # rod + bob
trail, = ax.plot([], [], 'r-', alpha=0.5)          # trail of bob2

# For storing the trail of the second bob
max_trail = 1000  # how many points to keep in the trail
trail_x = []
trail_y = []

def init():
    line.set_data([], [])
    trail.set_data([], [])
    return line, trail

def update(frame):
    # frame is the index in the time array
    x1_curr = x1[frame]
    y1_curr = y1[frame]
    x2_curr = x2[frame]
    y2_curr = y2[frame]

    # Update rod and bob positions
    line.set_data([0, x1_curr, x2_curr],
                  [0, y1_curr, y2_curr])

    # Update trail (store bob2 positions)
    trail_x.append(x2_curr)
    trail_y.append(y2_curr)
    # Keep the trail length to max_trail
    if len(trail_x) > max_trail:
        trail_x.pop(0)
        trail_y.pop(0)

    trail.set_data(trail_x, trail_y)

    return line, trail

ani = FuncAnimation(fig, update, frames=len(time_vals),
                    init_func=init, blit=True, interval=25)

plt.show()
