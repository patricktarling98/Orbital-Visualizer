import visualizer
import orbitals

if __name__ == "__main__":
    letter_to_number = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4, 'h': 5, 'i': 6, 'j': 7, 'k': 8, 'l': 9}

    input = input("Enter orbital to show, enter response as nlm (e.g. 4f0): ")

    n = int(input[0])
    l = letter_to_number.get(input[1])
    m = int(input[2])

    orbital = orbitals.Orbital(n=n, l=l, m=m, Z=1)
    visualizer.render(orbital=orbital)