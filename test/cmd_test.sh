# Divide the MoS2 bilatyer in the two single layers
# The plane cuts through z=5, perpenticular to the vertical axis, so the direction, i.e. family of planes, is 0 0 1 and the intercept, i.e. the specific plane we need, is at 0 0 5.
# This is one of the infinite possilbe choices.
../str_plane_cut.py -n 0 0 1 -p 0 0 5 --plot MoS2.vasp | add_elem2ase_vasp.sh
