from ase.io import read, write
import sys

# Turns an ASE trajectory (.traj) into an animated GIF (for slides / GitHub).
# The atoms are shown from the side and slightly tilted. Hg = grey-blue, Se = orange.
#
# usage : python make_gif_01.py  <file.traj>  <output.gif>
# e.g.  : python make_gif_01.py ../03_Hg_on_tSe/Hg_tSe_site1.traj Hg_on_tSe_relax.gif

traj_file = sys.argv[1] if len(sys.argv) > 1 else '../03_Hg_on_tSe/Hg_tSe_site1.traj'
gif_file = sys.argv[2] if len(sys.argv) > 2 else 'Hg_on_tSe_relax.gif'

images = read(traj_file, index=':')
# (no repeat : a 2x2 repeat would also show the periodic copies of Hg)
for atoms in images:
    atoms.set_constraint()      # do not draw the crosses of fixed atoms

print(f"{traj_file}: {len(images)} frames -> {gif_file}")
write(gif_file, images, rotation='-75x,10y', interval=120, show_unit_cell=0)
print("done")
