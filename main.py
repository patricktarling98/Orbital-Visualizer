import visualizer
import orbitals

if __name__ == "__main__":
    print("Starting")
    orbital = orbitals.Orbital(n=2, l=1, m=0, Z=1)

    visualizer.render(orbital=orbital) #REAL WORK