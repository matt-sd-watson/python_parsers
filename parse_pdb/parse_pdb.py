import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# the file path will need to contain the location of the protein databse file to parse
filepath = "5kkk.pdb"

# Part 1

def parsePDB(file):

  """
  Input parameter:
  @file: a filepath to the pdb file to parse

  Objects returned:
  @residues: a list of all protein residues found in CA atoms
  @coordinates_array:a list of x,y,z coordinates for the list of residues found in CA atoms
  @list_heterogens: a list of all heterogens found in the parsed file
  @residue_seq_number: a list of all index positions of the residues
  """

  residues = []
  residue_seq_number = []
  heterogens = []
  coordinates = []
  with open(file) as pdbfile:
      for line in pdbfile:
          if line[:6] == "ATOM  " and line[12:16] == " CA " and not line[16] == "B":
              residues.append(line[17:20])
              residue_seq_number.append(line[22:26])
              coordinates_string = line[30:54]
              splitted = [float(i) for i in coordinates_string.split()]
              coordinates.append(splitted)
          elif line[:6] == "HETATM":
              heterogens.append(line)
  coordinates_array = np.array(coordinates)
  residue_seq_number_array = np.array(residue_seq_number)
  return residues, coordinates_array, heterogens, residue_seq_number_array


residues, coordinates_array, list_heterogens, residue_numbers = parsePDB(filepath)

# Part 2


def residue_hydrophobicity(residue_list):

   """
   Input parameter:
   @residue_list: a list of all protein residues found in CA atoms

   Objects returned:
   @residue_merged = a list of residues with a boolean indicator for hydrophobicity
   """

   hydrophobic_residues = ['ALA', 'CYS', 'PHE', 'ILE', 'LEU', 'MET', 'PRO', 'VAL', 'TRP']
   residue_flag = []
   for line in residue_list:
       if line in hydrophobic_residues:
           residue_flag.append(True)
       else:
           residue_flag.append(False)
   flags_array = np.array(residue_flag)
   residue_array = np.array(residue_list)
   residue_merged = np.column_stack((residue_array, flags_array))
   return residue_merged


def distance2Center(coordinates_array):

  """
   Input parameter:
   @coordinates_array: a list of x,y,z coordinates for the list of residues corresponding to CA atoms

   Objects returned:
   @distance_to_center = a list of respective distances from the mean for each of the
   residues input into the function
   """
  central_average = np.sum(coordinates_array, axis=0) / coordinates_array.shape[0]
  squared_difference = (central_average - coordinates_array) ** 2
  squared_dist = np.sum(squared_difference, axis=1)
  distance_to_center = np.sqrt(squared_dist)
  return (list(distance_to_center))

output = np.column_stack((residue_hydrophobicity(residues), distance2Center(coordinates_array)))
write_to_file = pd.DataFrame(output, columns=["Residue", "Hydrophobicity", "Distance to Centre"])
write_to_file.to_csv("output_part_ii.csv", index=False)




def drawGraph(residues, coordinates_array):

  """
   Input parameter:
   @residues: a list of all protein residues found in CA atoms
   @coordinates_array: a list of x,y,z coordinates for the list of residues

   Output:
   @plt.show() - a grouped boxplot showing the average distance to the centre for both
   hydrophobic and hydrophilic residues
   """

  hydrophobic_table = []
  hydrophilic_table = []
  table = np.column_stack((residue_hydrophobicity(residues), distance2Center(coordinates_array)))
  for i in range(0, table.shape[0]):
      if (table[i, 1] == 'True'):
          hydrophobic_table.append(table[i, 2])
          hydrophobic_distance = [float(i) for i in hydrophobic_table]
      else:
          hydrophilic_table.append(table[i, 2])
          hydrophilic_distance = [float(i) for i in hydrophilic_table]
  box_plot_data = [hydrophobic_distance, hydrophilic_distance]
  plt.boxplot(box_plot_data, patch_artist=True, labels=['Hydrophobic AA', 'Hydrophilic AA'])
  plt.show()


drawGraph(residues, coordinates_array)

# Part 3


def searchHeterogen(heterogen_list, atom_name):
  """

   Input parameter:
   @list_heterogens: a list of all heterogens found in the parsed file
   @atom_name: a string representing the name of an atom in the PDB file i.e. "FE"

   Output:
   @heterogen_coordinates: the x,y,z coordinates of the specific atom requested in the input
   """

  heterogen_position = []
  for line in heterogen_list:
      if "HETATM" and atom_name in line:
          heterogens = line[30:54]
          heterogen_float = [float(i) for i in heterogens.split()]
          heterogen_position.append(heterogen_float)

  heterogen_coordinates = np.array(heterogen_position)
  return heterogen_coordinates


# To get iron coordinates
fe_coordinates = searchHeterogen(list_heterogens, "FE")



def distance_from_heterogen(coordinate_array, heterogen_coordinates):
  """

   Input parameter:
   @coordinates_array:a list of x,y,z coordinates for the list of residues
   @heterogen_coordinates: the x,y,z coordinates of the specific atom requested in searchHeterogen
   i.e. "FE"

   Objects returned:
   @distance_to_center = a list of respective distances from the heterogen centre for each of the
   residues input into the function
   """

  squared_difference = (heterogen_coordinates - coordinate_array) ** 2
  squared_dist = np.sum(squared_difference, axis=1)
  dist = np.sqrt(squared_dist)
  distance_to_heterogen = np.array(dist)
  return distance_to_heterogen


# Create a csv output file containing the residue index, residue name, hydrophobicity, and distance to
# the FE atom
output_2 = np.column_stack((residue_numbers, residue_hydrophobicity(residues), distance_from_heterogen(coordinates_array, fe_coordinates)))
output_2_to_file = pd.DataFrame(output_2,columns=["Residue Index", "Residues", "Hydrophobicity", "Distance to FE Atom"])
output_2_to_file.to_csv("output_part_iii.csv", index=False)

# To print the 5 closest residues to the iron (FE) atom
data = pd.read_csv("output_part_iii.csv")
least5 = data.nsmallest(5, "Distance to FE Atom")
print(least5)


