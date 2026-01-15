import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Install with: pip install fipy
from fipy import CylindricalGrid1D, CellVariable, TransientTerm, DiffusionTerm, ImplicitSourceTerm

# Physical parameters (converted to mm units)
D = 3e-3  # Oxygen diffusion coefficient in water (fmols/mm^2/s)
ocr = 1  # Oxygen consumption rate (fmols/mm2/s)
C_air = 200.0  # Oxygen concentration at air interface (µM, ~8.5 mg/L)

# Geometry
L = 3.1  # Depth of cylindrical well (mm)
nr = 100  # Number of mesh points
dx = L / nr
# Create 1D cylindrical mesh (radial direction)
# For simplicity, we'll model the axial (depth) direction
# since this 1D grid uses only "left" and "right" boundaries
mesh = CylindricalGrid1D(dx=dx, nx=nr)

C_initial = C_air #well is initially saturated with O2
# Create the concentration variable
C = CellVariable(name="oxygen concentration", mesh=mesh, value=C_initial)

# Boundary conditions
#NOTE that CylindricalGrid1D has only "left" and "right" faces for boundary conditions
# we will use "right" (will last index in C values) as our air-liquid interface
# (top surface of media in the well)
#  we will use "left" (first index in C  values array) as the bottom of the
#  well - where we assume no diffusion can occur

# Top surface right: Contact with air - fixed concentration
C.constrain(C_air, mesh.facesRight)  # Air interface

# Bottom surface "left" (z=L): Sealed - zero flux (natural BC)
C.faceGrad.constrain((0,), where=mesh.facesLeft)

# Time parameters
dt = 0.9 * dx**2 / (2 * D)
steps = 2000  # Number of time steps
total_time = dt * steps

# Define the PDE: dC/dt = D * d²C/dz² - k*C
# TransientTerm for dC/dt
# DiffusionTerm for D * d²C/dz²
# ImplicitSourceTerm for -k*C (consumption)
k = D*ocr

eq = (TransientTerm() == DiffusionTerm(coeff=D) - ImplicitSourceTerm(coeff=k))

# Storage for plotting
times_to_plot = [0, 100, 500, 1000]
profiles = []
times = []

cs_at_z = []
# Solve
print("Solving oxygen diffusion with consumption...")
for step in range(steps + 1):
    profiles.append(C.value.copy())
    times.append(step * dt)
    print(f"Step {step}/{steps}, Time: {step*dt:.1f} s")

    if step < steps:
        eq.solve(var=C, dt=dt)
        print(f'c0 {C.value[0]} cN {C.value[-1]}')

# Print some statistics
print(f"\n=== Simulation Results ===")
print(f"Well depth: {L:.1f} mm")
print(f"Diffusion coefficient: {D:.2e} mm²/s")
print(f"Consumption rate: {k:.2e} 1/s")
print(f"Air O₂ concentration: {C_air:.1f} µM")
print(f"Final surface concentration: {C.value[0]:.1f} µM")
print(f"Final bottom concentration: {C.value[-1]:.1f} µM")
print(f"Penetration depth (90% depletion): {np.where(C.value < 0.1*C_air)[0][0] * L/nr:.2f} mm"
      if np.any(C.value < 0.1*C_air) else "Full penetration")

# Calculate characteristic length scale
L_char = np.sqrt(D / k)  # mm
print(f"Characteristic diffusion-reaction length: {L_char:.2f} mm")
print(f"Damköhler number (Da = L²k/D): {(L**2 * k / D):.3f}")
