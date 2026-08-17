from ase.build import bulk
from ase.mep import NEB, NEBTools
from ase.optimize import BFGS

from mace.calculators import mace_mp

calc = mace_mp(
    model="medium",
    dispersion=False,
    default_dtype="float64",
    device="cpu",
)

''' Test structures
bulk_cu = bulk("Cu", "fcc", a=3.615, cubic=True) * (3, 3, 3)
initial = bulk_cu.copy()
del initial[0]  # vacancy at site 0

final = bulk_cu.copy()
del final[1]    # vacancy hopped to neighbouring site 1
''' 
# Relax the two endpoints
# for atoms, name in [(initial, "initial"), (final, "final")]:
#     atoms.calc = calc
#     BFGS(atoms, logfile=f"{name}_relax.log").run(fmax=0.02)

# Build the band and interpolate a starting guess
n_images = 10  # number of intermediate images
ini
images = [initial] + [initial.copy() for _ in range(n_images)] + [final]

neb = NEB(images, method="improvedtangent")
neb.interpolate(method="idpp")  # image-dependent pair potential: better guess than linear

for image in images[1:-1]:
    image.calc = calc

opt = BFGS(neb, trajectory="neb.traj", logfile="neb.log")
opt.run(fmax=0.10)

neb.climb = True
opt.run(fmax=0.05)

# Read off the barrier
nebtools = NEBTools(images)
Ef, dE = nebtools.get_barrier()
print(f"Forward migration barrier: {Ef:.3f} eV")
print(f"Reaction energy (final - initial): {dE:.3f} eV")

fig = nebtools.plot_band()
fig.savefig("neb_band.png")
print("Band plot saved to neb_band.png")
