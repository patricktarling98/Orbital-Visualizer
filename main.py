import visualizer
import orbitals
"-----------------------------------------------PARAMS-------------------------------------------------------------------------------"
numPoints = visualizer.numPoints
"--------------------------------------------------------------------------------------------------------------------------------------"

if __name__ == "__main__":
    print("Starting")
    orbital = orbitals.Orbital(n=1, l=0, m=0, Z=1)

    visualizer.render(orbital=orbital, numPoints=numPoints) #REAL WORK