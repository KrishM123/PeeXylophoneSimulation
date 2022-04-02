"""
By: Krish Modi
Date: October 30 2021

This is a particle motion simulator for Science Fair 2022
Input the initial values and the program will perform the motion experienced by the urine
This program works on the Simulation Stepper concept where the characteristics of the urine projectile will update every time step
"""
import math
import matplotlib.pyplot as plt

# The following INIT variables are the initial conditions of the urine projectile
# Change the following 4 variables to test different conditions
INIT_VELOCITY = float(input("\n\nEnter initial velocity (m/s): "))  # m/s
INIT_ANGLE = float(input("Enter initial angle (degrees): "))  # deg
INIT_HEIGHT = float(input("Enter initial height (m): "))  # m

# Smaller TIME_STEP results in higher accuracy
TIME_STEP = 0.00001  # s

# DO NOT CHANGE ANYTHING FROM AFTER THIS POINT

# The initial coordinates of the urine projectile
V_POS = INIT_HEIGHT
H_POS = 0  # m

TIME_ELAPSED = 0
all_h_pos = []  # Stores all x values to be plotted
all_v_pos = []  # Stores all y values to be plotted

# Extracts the X and Y velocity components of INIT_VELOCITY
H_VEL = math.cos(math.radians(INIT_ANGLE)) * INIT_VELOCITY
V_VEL = math.sin(math.radians(INIT_ANGLE)) * INIT_VELOCITY

# The only acceleration acting upon the projectile is gravity in the vertical direction
V_ACC = -9.81

# This loop runs until the projectile hits the ground
while True:
    # Stores all occupied positions to be plotted later
    all_h_pos.append(H_POS)
    all_v_pos.append(V_POS)
  
    # Updates the position of the projectile based on the velocity as the time
    H_POS += H_VEL * TIME_STEP
    V_POS += V_VEL * TIME_STEP

    # Updates the vertical velocity based on gravity
    V_VEL += V_ACC * TIME_STEP

    # TIME_ELAPSED will show how long the projectile took to hit the ground
    TIME_ELAPSED += TIME_STEP

    # If the projectile is on the ground, stop the loop
    if V_POS < 0:
        break
      
# The following 4 lines basically rename variables to match the variable names from the related derivations
R = H_POS
g = 9.81
h = INIT_HEIGHT
V3 = V_VEL
print(f"\nImpact velocity is: {V_VEL}. Range is: {H_POS}")

# Checks the validity of the V2 calculations from Derivation 1.docx
# print(f"The calculated initial velocity assuming initial angle is 0 is: {(round(R, 2) * g) / math.sqrt(2 * g * h)}")

# Creates the two functions required for Newton's Method from Derivation 2.docx
fu = lambda V: -(math.sqrt((math.sin(math.atan((V ** 2 - math.sqrt(V ** 4 - R ** 2 * g ** 2 + 2 * g * h * V ** 2)) / (R * g))) * V) ** 2 + (2 * g * h))) - V3
# The following function is the derivative of fu
fup = lambda V: (V * math.sqrt(-(R ** 2 * math.sqrt(V ** 4 - R ** 2 * g ** 2 + 2 * g * h * V ** 2) + (-2 * h ** 2 - R ** 2) * V ** 2 - 4 * g * h ** 3 - 3 * R ** 2 * g * h) / (2 * h ** 2 + 2 * R ** 2)) * (R ** 2 * g * math.sqrt(V ** 4 - R ** 2 * g ** 2 + 2 * g * h * V ** 2) - 2 * h * V ** 4 - 4 * g * h ** 2 * V ** 2 + 2 * R ** 2 * g ** 2 * h)) / ((2 * h * V ** 2 + 4 * g * h ** 2 + R ** 2 * g) * (V ** 4 - R ** 2 * g ** 2 + 2 * g * h * V ** 2))

fd = lambda V: -(math.sqrt((math.sin(math.atan((V ** 2 + math.sqrt(V ** 4 - R ** 2 * g ** 2 + 2 * g * h * V ** 2)) / (R * g))) * V) ** 2 + (2 * g * h))) - V3
# The following function is the derivative of fd
fdp = lambda V: -(V * math.sqrt((R ** 2 * math.sqrt(V ** 4 + (2 * g * h * V ** 2) - (R ** 2 * g ** 2))) + ((2 * h ** 2 + R ** 2) * V ** 2) + (4 * g * h ** 3) + (3 * R ** 2 * g * h)) * ((R ** 2 * g * math.sqrt(V ** 4 + (2 * g * h * V ** 2) - (R ** 2 * g ** 2))) + (2 * h * V ** 4) + (4 * g * h ** 2 * V ** 2) - (2 * R ** 2 * g ** 2 * h))) / (math.sqrt(2) * math.sqrt(h ** 2 + R ** 2) * ((2 * h * V ** 2) + (4 * g * h ** 2) + (R ** 2 * g)) * (V ** 4 + (2 * g * h * V ** 2) - (R ** 2 * g ** 2)))

# This is the first guess that will be inputted in Newton's Method
oldV = math.sqrt(((-2 * g * h) + math.sqrt(((2 * g * h) ** 2) - (4 * ((-R ** 2) * (g ** 2))))) / 2)

# This is Newton's Method applied
f = fu
p = fup
while True:
    try:
        if fd(oldV) >= 0:
            f = fd
            p = fdp
        break
    except ValueError:
        oldV = oldV * 1.000000001

while True:
    try:
        newV = oldV - (f(oldV) / p(oldV))
        if abs(oldV - newV) < 0.000000001:
            break
        oldV = newV
    except ValueError:
        oldV = oldV * 1.000000001

# Checks the validity of the V2 calculations from Derivation 2.docx
print(f"The calculated initial velocity is: {newV}")


# The following lines show the trajectory of the droplet in a visual manner
plt.plot(all_h_pos, all_v_pos)
plt.xlabel("Horizontal Position (m)")
plt.ylabel("Vertical Position (m)")
plt.title("Trajectory of Urine Droplet")
plt.show()