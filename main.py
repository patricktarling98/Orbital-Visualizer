import visualizer
import atoms
"-----------------------------------------------PARAMS-------------------------------------------------------------------------------"
numPoints = 10000
atom_to_show = "Hydrogen"
"--------------------------------------------------------------------------------------------------------------------------------------"

if __name__ == "__main__":
    print("Starting")
    atom = atoms.Atom(atom_to_show)

    visualizer.render(numPoints=numPoints, atom=atom) #REAL WORK