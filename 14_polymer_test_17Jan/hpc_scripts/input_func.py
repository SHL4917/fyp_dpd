import numpy as np

def create_poly_datafile(box_size, polymer_conc, polymer_structure, 
                         filename = "data.txt", mass_array = [1, 1.12, 0.866], 
                         rho_star = 3):
    def got_space(x, y, z, matrix, polymer_length):
        space = True
        for i in range(polymer_length):
            if matrix[x + i][y][z] != 0:
                space = False
        
        return space
    
    
    polymer_length = len(polymer_structure)
    polymer_structure = polymer_structure.astype(int)
    description = str(polymer_length) + "_" + str(polymer_conc)
    
    angles = 0
    dihedrals = 0
    impropers = 0
    
    # One more atom type to allow for srp between bonds
    atom_types = len(mass_array) + 1
    bond_types = 1
    angle_types = 0
    dihedral_types = 0
    improper_types = 0
    
    # Calculate simulation parameters
    particle_count = rho_star * box_size**3
    
    polymer_beadcount = np.round((particle_count * polymer_conc)).astype(int)
    num_chains = np.round((polymer_beadcount/polymer_length)).astype(int)
    polymer_beadcount = num_chains * polymer_length
    
    num_solvent = particle_count - polymer_beadcount
    
    gridsize = np.ceil((particle_count ** (1/3))).astype(int)
    spacing = box_size / (gridsize - 1)
    
    matrix = np.zeros((gridsize, gridsize, gridsize))
    polymer_coords = np.zeros((num_chains, 3))
    # Set the coordinates for the polymer heads and other chains
    # Polymer chains are initially all oriented in the x-direction
    
    """
    Sets the polymer bead type array such that the head has value of 69, such that 
    when filtering the matrix values in writing to input file, the polymer head 
    will not be mis-attributed
    
    e.g. when polymer_structure = [2, 2, 3, 3]
    
    """
    polymer_encoding = polymer_structure.copy()
    polymer_encoding[0] = 69
    
    
    i = 0 
    while i < num_chains:
        x = np.random.randint(0, gridsize - polymer_length)
        y = np.random.randint(0, gridsize - 1)
        z = np.random.randint(0, gridsize - 1)
        
        if got_space(x, y, z, matrix, polymer_length):
            for j in range(polymer_length):
                matrix[x + j][y][z] = polymer_encoding[j]
                
            polymer_coords[i][0] = x
            polymer_coords[i][1] = y
            polymer_coords[i][2] = z
            
            i = i + 1
            
    # Populate remaining with solvent beads
    i = 0
    for x in range(gridsize):
        for y in range(gridsize):
            for z in range(gridsize):
                if (matrix[x][y][z] == 0) and (i < num_solvent):
                    matrix[x][y][z] = 1
                    i = i + 1
    
    
    with open(filename, mode = "w") as line:
        line.write("LAMMPS " + description + "\n \n")
        
        line.write(str(particle_count) + " atoms\n")
        line.write(str(num_chains * (polymer_length - 1)) + " bonds\n")
        line.write(str(angles) + " angles\n")
        line.write(str(dihedrals) + " dihedrals\n")
        line.write(str(impropers) + " impropers\n \n")
        
        line.write(str(atom_types) + " atom types\n")
        line.write(str(bond_types) + " bond types\n")
        line.write(str(angle_types) + " angle types\n")
        line.write(str(dihedral_types) + " dihedral types\n")
        line.write(str(improper_types) + " improper types\n \n")
        
        line.write("0 " + str(box_size) + " xlo xhi\n")
        line.write("0 " + str(box_size) + " ylo yhi\n")
        line.write("0 " + str(box_size) + " zlo zhi\n \n")
        
        line.write("Masses\n \n")
        for index, val in enumerate(mass_array):
            line.write(str(index + 1) + " " + str(val) + "\n")
    
        line.write(str(atom_types) + " 1\n")
        line.write("\n")
        
        line.write("Atoms\n \n")
        i = 1
        j = 0
        polymer_number = np.zeros((num_chains, polymer_length))
        for index, val in np.ndenumerate(matrix):
            x = index[0] * spacing
            y = index[1] * spacing
            z = index[2] * spacing
            if val == 1:
                line.write(str(i) + " 0 1 " + 
                           "{:.3f}".format(x) + " " + "{:.3f}".format(y) + " " 
                           + "{:.3f}".format(z) + "\n")
                i = i + 1
            elif val == 69:
                for k in range(polymer_length):
                     line.write(str(i) + " 0 " + str(polymer_structure[k]) + 
                                " " + "{:.3f}".format(x + spacing * k) + " " + 
                           "{:.3f}".format(y) + " " + "{:.3f}".format(z) + 
                           "\n")
                     polymer_number[j][k] = i
                     i = i + 1
                j = j + 1
                
            
        line.write("\n")
        
        line.write("Bonds\n \n")
        i = 1
        bond_per_chain = polymer_length - 1
        polymer_number = polymer_number.astype(int)
        for row in polymer_number:
            for j in range(bond_per_chain):
                line.write(str(i) + " 1 " + str(row[j]) + " " + str(row[j + 1])
                           + "\n")
                i = i + 1
                
    return polymer_number
    
    

        

