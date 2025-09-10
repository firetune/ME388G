import math


mass_neutron_amu = 1.008664904 # amu
mass_deuteron_amu = 2.01410178 # amu
mass_triton_amu = 3.01604927 # amu
KE_neutron_MeV = 4.0 # MeV
c_squared_MeV_amu = 931.49410242 # MeV*amu
c_MeV_amu = math.sqrt(c_squared_MeV_amu) # MeV^0.5 * amu^0.5
theta_deg = 30.0 # degrees

Q = (mass_neutron_amu + mass_deuteron_amu - mass_triton_amu) * c_squared_MeV_amu
print(f"Q value: {Q:.6f} MeV")

a = 1 / (2 * mass_triton_amu * c_squared_MeV_amu)
print(f"a: {a:.6f} 1/MeV")

b = (1 - ((math.sqrt(2 * mass_neutron_amu * KE_neutron_MeV) * math.cos(math.radians(theta_deg))) / (mass_triton_amu * c_MeV_amu)))
print(f"b: {b:.6f}")

c = KE_neutron_MeV * ((mass_neutron_amu / mass_triton_amu) - 1) - Q
print(f"c: {c:.6f} MeV")

##Solve quadratic equation
discriminant = b**2 - 4*a*c
if discriminant < 0:
    print("No real solution for the kinetic energy of the triton.")
else:
    KE_triton_1 = (-b + math.sqrt(discriminant)) / (2*a)
    KE_triton_2 = (-b - math.sqrt(discriminant)) / (2*a)
    print(f"Kinetic energy of triton solutions: {KE_triton_1:.6f} MeV, {KE_triton_2:.6f} MeV")  