Copyright 2023 Novartis Institutes for BioMedical Research Inc.
 
Licensed under the MIT License (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
 
https://www.mit.edu/~amini/LICENSE.md
 
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.


from rdkit import Chem
from rdkit.Chem import Descriptors,Draw
import base64


#check if there is N-nitroamine pattern in the molecule, if yes return with every pair of N index and alpha carbons' index
def is_nitro_pattern(mol):
    key_list = []
    nitro_pattern_smile = 'N(N=O)'
    nitro_pattern = Chem.MolFromSmarts(nitro_pattern_smile)
    matches = mol.GetSubstructMatches(nitro_pattern)
    if matches:
        for match in matches:
            nitrogen = mol.GetAtomWithIdx(match[0])
            if nitrogen.GetDegree() == 3:
                neighbors_symbol = [a.GetSymbol() for a in nitrogen.GetNeighbors()]
                neighbors =  [a for a in nitrogen.GetNeighbors()]
                if 'C' in neighbors_symbol and neighbors_symbol.count('C') == 2:
                    carbon_idx = []
                    for neighbor in neighbors:
                        if neighbor.GetSymbol() == 'C':
                            carbon_idx.append(neighbor.GetIdx())
                    key_list.append([match[0],carbon_idx])
    return key_list


# You could uncomment this part and commment out the previous part to enable NH and N(CH3)2 detecting

# def is_nitro_pattern(mol):
#     nitrogen = []
#     key_list = []
#     for atom in mol.GetAtoms():
#         if atom.GetSymbol() == 'N':
#             nitrogen.append(atom)
#     if len(nitrogen)>0:
#         for nitrogen in nitrogen:
#             #check if the nitrogen is in a aromatic ring, if yes, it's not a nitrosamine pattern
#             if not nitrogen.GetIsAromatic():
#                 carbons = []
#                 carbons_idx = []
#                 hydrogen = []
#                 methyl_c = []
#                 for neighbor in nitrogen.GetNeighbors():
#                     if neighbor.GetSymbol() == 'C':
#                         # find methyl group connected to nitrogen
#                         if neighbor.GetDegree() == 4:
#                             neighbors = [a.GetSymbol() for a in neighbor.GetNeighbors()]
#                             if 'H' in neighbors and neighbors.count('H') == 3:
#                                 methyl_c.append(neighbor)
#                         carbons.append(neighbor)
#                         carbons_idx.append(neighbor.GetIdx())
#                     elif neighbor.GetSymbol() == 'H':
#                         hydrogen.append(neighbor)
#                     # Detect Nitrosamine Group(N=0)
#                     elif neighbor.GetSymbol() == 'N':
#                         for a in neighbor.GetNeighbors():
#                             if a.GetSymbol() == 'O':
#                                  bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), a.GetIdx())
#                                  if  bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
#                                      hydrogen.append(neighbor)
#                 # if N connected with 2 CH3, then it's a nitrosamine pattern
#                 if len(methyl_c)==2:
#                     hydrogen.append(methyl_c[0])
#                     carbons.remove(methyl_c[0])
#                     carbons_idx.remove(methyl_c[0].GetIdx())
#                 # if a nitrosamine pattern is found
#                 if len(carbons) == 2 and len(hydrogen) == 1:
#                     double_bonded_to_heteroatom = False
#                     for carbon in carbons:
#                         for neighbor in carbon.GetNeighbors():
#                             if neighbor.GetSymbol() != 'H' and neighbor.GetSymbol() != 'C':
#                                 bond = mol.GetBondBetweenAtoms(carbon.GetIdx(), neighbor.GetIdx())
#                                 if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
#                                     double_bonded_to_heteroatom = True
#                                     break
#                         if double_bonded_to_heteroatom:
#                             break
#                     if not double_bonded_to_heteroatom:
#                         key_list.append([nitrogen.GetIdx(),carbons_idx])
#     return key_list


# check alpha carbon hydrogen
def hydrogen_num(mol,list):
    alpha_h = []
    for carbon in list:
        atom = mol.GetAtomWithIdx(carbon)
        h_num = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 1)
        alpha_h.append(h_num)
    return alpha_h

#check tertiary_carbon
def tertiary_carbon(mol,list):
    for carbon in list:
        atom = mol.GetAtomWithIdx(carbon)
        if atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3 and sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 6) == 3:
            # print('yes')
            return True
    return False

# Find caboxylic_acid group in the molecule
def caboxylic_acid(mol):
    acid = None
    carboxylic_acid = Chem.MolFromSmarts('C(=O)O[H]')
    match_pattern = mol.GetSubstructMatch(carboxylic_acid)
    print(match_pattern)
    if mol.HasSubstructMatch(carboxylic_acid):
        single_oxygen = mol.GetAtomWithIdx(match_pattern[2])
        print(single_oxygen.GetSymbol())
        for neighbor in single_oxygen.GetNeighbors():
            print(neighbor.GetSymbol())
            if neighbor.GetSymbol() == 'H':
                acid = True
        # return False
    return acid

# pyrrolidine ring is found or not
def is_in_pyrrolidine(mol,N_index):
    pyrrolidine = Chem.MolFromSmarts('C1CCNC1')
    matches = mol.GetSubstructMatches(pyrrolidine)
    is_pyrrolidine = False
    if matches:
        # print("Pyrrolidine ring is present in the molecule")
        for match in matches:
            if N_index in match:
                is_pyrrolidine = True
            else:
                continue
    return is_pyrrolidine


## Get A ring containing a specific atom
def get_ring_atoms(mol, atom_idx):
    for ring in mol.GetRingInfo().AtomRings():
        if atom_idx in ring:
            return ring
    return None

def is_in_morpholine(mol,N_index):
# create a substructure object for the morpholine ring
    substructure = Chem.MolFromSmarts('C1COCCN1')

# get the indices of all atoms in the molecule that match your substructure
    atom_indices = mol.GetSubstructMatches(substructure)
    is_morpholine = False
    if atom_indices:
        for atom_indice in atom_indices:
            if N_index in atom_indice:
                is_morpholine = True
    return is_morpholine

##  Check Ring Score (pyrrolidine ring, ring with sulfur, morphline ring, 7-membered ring)

def get_ring_related_score(mol,N_index):
    de_score = 0
    nitroso_6_sulfur_ring = False
    atom = mol.GetAtomWithIdx(N_index)
    if atom.IsInRingSize(5):
        de_score = 2
        is_pyrrolidine = is_in_pyrrolidine(mol,N_index)
        if is_pyrrolidine:
            de_score = 3
    elif atom.IsInRingSize(6):
        de_score = 2
        ring = get_ring_atoms(mol, N_index)
        for atom_idx in ring:
            r_atom = mol.GetAtomWithIdx(atom_idx)
            if r_atom.GetSymbol()=='S':
                de_score = 3
        if is_in_morpholine(mol,N_index):
            de_score = 1
    elif atom.IsInRingSize(7):
        de_score = 1
    return de_score

def find_chains(mol, atom, chain):
    # base case: if the chain is empty, add the current atom
    if not chain:
        chain.append(atom)
    # base case: if the chain has reached the desired length, return it
    if len(chain) >=6:
        return [chain]
    # recursive case: loop through the neighbors of the current atom
    chains = []
    for neighbor in mol.GetAtomWithIdx(atom).GetNeighbors():
        # check if the neighbor is not already in the chain
        if neighbor.GetSymbol() !='H' and neighbor.GetIdx() not in chain:
            # make a copy of the chain and append the neighbor
            new_chain = chain.copy()
            new_chain.append(neighbor.GetIdx())
            # recursively find the chains from the neighbor
            chains.extend(find_chains(mol, neighbor.GetIdx(), new_chain))
    return chains


def check_atoms_in_same_ring(mol,atom_indices):

    atoms = [mol.GetAtomWithIdx(i) for i in atom_indices]
    # check if all atoms are in a ring
    if all(atom.IsInRing() for atom in atoms):
        # check if all atoms are in the same ring
        ring_info = mol.GetRingInfo()
        if all(ring_info.AreAtomsInSameRing(atoms[0].GetIdx(), atom.GetIdx()) for atom in atoms):
            # check if ring size is smaller than 8
            # ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
            for ring in ring_info.AtomRings():
                if all(atom_indices[i] in ring for i in range(len(atom_indices))):
                    if len(ring) < 8:
                        return True
                    else:
                        return False
        else:
            return False
    else:
        return False


# make sure that the chains are on both side, do not share the same atom
def common_data(list1, list2):
    for x in list1:
        for y in list2:
            if x == y:
                return True
    return False

def check_sublist(list1, list2):
    for sublist1 in list1:
        for sublist2 in list2:
            if not common_data(sublist1, sublist2):
                return True
    return False


# Find the chains of non-hydrogen atoms that are connected to each carbon atom

def five_consec_atoms(mol, N_index, carbon_index, chain):
    de_score = 0
    ri_info = mol.GetRingInfo()
    if ri_info.NumAtomRings(N_index) > 0 and ri_info.AtomRingSizes(N_index)[0] < 8:
        return de_score
    result = find_chains(mol, N_index, chain)
    side_chain = []
    for list in result:
        list.remove(N_index)
        side_chain.append(list)
    carbon_1=[]
    carbon_2=[]
    carbon_1_ring = True
    carbon_2_ring = True
    for list in side_chain:
        if carbon_index[0] in list:
            carbon_1.append(list)
        elif carbon_index[1] in list:
            carbon_2.append(list)
    if carbon_1 and carbon_2 and check_sublist(carbon_1,carbon_2):
        for ring in carbon_1:
            if not check_atoms_in_same_ring(mol,ring):
                carbon_1_ring = False
        for ring in carbon_2:
            if not check_atoms_in_same_ring(mol,ring):
                carbon_2_ring = False
    if carbon_1_ring == False and carbon_2_ring == False:
        de_score = 1
    return de_score

str_EWG_group = ['CC(F)(F)(F)','CC#N','C[N+]','[O-][N+](C)=O','O=S(C)']
R2 = ['O','N','S']
R3 = ['C','F','Br','I','OC','N','N(C)','N(C)C']

# Initialize an empty list to store the smiles patterns
carb_de_smiles = []

# Loop through each element in R2 and R3
for r2 in R2:
  for r3 in R3:
    # Concatenate the strings to form a smiles pattern
    smile = r2 + "=C" +'(C)' + r3
    # Append the smile to the list
    carb_de_smiles.append(smile)
carb_de_smiles.remove('O=C(C)C')
# # Print the list of smiles patterns
# print(carb_de_smiles)

# Halogen group
halogen_smiles = ['C(F)','C(Cl)','C(Br)','C(I)']
aromatic_smiles = ['CO','CS','CN','[#6]/[#6]([H])=[#6]([H])','CC#C']
michael_smiles = ['[#6]/[#6]([H])=[#6]([H])/[#6](=O)','[#6]/[#6]([H])=[#6]([H])/[#6](=N)','[#6]/[#6]([H])=[#6]([H])/[#6](=S)','[#6]/[#6]#[#6]/[#6](=O)','[#6]/[#6]#[#6]/[#6](=N)','[#6]/[#6]#[#6]/[#6](=S)',
                  '[#6]/[#6]([H])=[#6]/C#N','CC#CC#N']


com_ewg_group = str_EWG_group + carb_de_smiles + halogen_smiles + aromatic_smiles + michael_smiles

def identify_ewg_alpha(mol,N_index,single_carbon_index):
    ri = mol.GetRingInfo()
    de_score = 0
    for smile in com_ewg_group:
        substructure = Chem.MolFromSmarts(smile)
        atom_indices = mol.GetSubstructMatches(substructure)
        if atom_indices:
            for atom_indice in atom_indices:
                if N_index not in atom_indice and single_carbon_index in atom_indice:
                    for k, index in enumerate(atom_indice):
                        if single_carbon_index == index:
                            if smile in aromatic_smiles:
                                if smile == '[#6]/[#6]([H])=[#6]([H])':
                                    last_atom = atom_indice[-2]
                                else:
                                    last_atom = atom_indice[-1]
                                last_atom = mol.GetAtomWithIdx(last_atom)
                                for neighbor in last_atom.GetNeighbors():
                                    if neighbor.GetSymbol() == 'C' and neighbor.GetIsAromatic() and not ri.AreAtomsInSameRing(single_carbon_index, neighbor.GetIdx()):
                                        de_score = 1
                            else:
                                de_score = 1
    return de_score

def hydroxyl_beta_carbon(mol,single_carbon_index):
      # Get the atom object from the molecule using the index
    de_score = 0
    alpha_carbon = mol.GetAtomWithIdx(single_carbon_index)
    for neighbor in alpha_carbon.GetNeighbors():
  # Check if the atom is carbon
        if neighbor.GetAtomicNum() == 6:
    # Initialize a variable for the oxygen atom
            oxygen_atom = None
    # Loop through the neighbors of the atom
            for neighbor_beta in neighbor.GetNeighbors():
      # Get the atomic number of the neighbor
                neighbor_beta_num = neighbor_beta.GetAtomicNum()
      # If the neighbor is oxygen, assign it to the oxygen atom variable
                if neighbor_beta_num == 8:
                    oxygen_atom = neighbor_beta
    # Check if the oxygen atom is not None

      # Loop through the neighbors of the oxygen atom
                    for neighbor_h in oxygen_atom.GetNeighbors():
                        if neighbor_h.GetAtomicNum() == 1 and neighbor.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                            de_score = 1
    return de_score

#check the atoms are in the same rings,
def are_atoms_in_same_ring(mol,atom_1, atom_2):
    ring_info = mol.GetRingInfo()
    atom_1_index = atom_1.GetIdx()
    atom_2_index = atom_2.GetIdx()
    ring_idx_list = ring_info.AtomRings()
    atom1_rings = []
    atom2_rings = []
    for ring in ring_idx_list:
        if atom_1_index in ring:
            atom1_rings.append(ring)
    for ring in ring_idx_list:
        if atom_2_index in ring:
            atom2_rings.append(ring)

    if not atom1_rings or not atom2_rings:
        return False
    for ring in atom1_rings:
        if ring in atom2_rings:
            return True
    return False

# check aromatci bonds
def find_aromatic_bonds(mol, atom_indices):
    aromatic_bonds = []
    for bond in mol.GetBonds():
        if bond.GetIsAromatic():
            atom1_idx = bond.GetBeginAtom().GetIdx()
            atom2_idx = bond.GetEndAtom().GetIdx()
            if atom1_idx in atom_indices and atom2_idx in atom_indices:
                aromatic_bonds.append(bond)
    return aromatic_bonds

# check double bounds
def has_double_bond(atom):
    for bond in atom.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            return True
    return False


def aryl_ring_atoms(mol, atom_idx, single_carbon_index):
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        if atom_idx in ring and single_carbon_index not in ring:
            return ring

def is_ring_connected(mol, atom_idx):
    neighbors = mol.GetAtomWithIdx(atom_idx).GetNeighbors()
    neighbors_list = []
    for neighbor in neighbors:
        # print(neighbor.GetIdx())
        if neighbor.IsInRing():
            neighbors_list.append(neighbor.GetIdx())
    return neighbors_list

#Find aryl_group
def aryl_group_check(mol,single_carbon_index):
    aromatic_num = 0
    Chem.Kekulize(mol)
    aromatic_ring = False
    if mol.GetAtomWithIdx(single_carbon_index).GetIsAromatic():
        aromatic_ring = False
    elif is_ring_connected(mol, single_carbon_index):
        connected_idx = is_ring_connected(mol, single_carbon_index)
        for k in connected_idx:
            ring_atoms = aryl_ring_atoms(mol, k, single_carbon_index)
            if ring_atoms:
                aromatic_ring = True
                for at_idx in ring_atoms:
                    if not mol.GetAtomWithIdx(at_idx).GetIsAromatic():
                        aromatic_ring = False
    if aromatic_ring == True:
        return True
    return False

# Define a function to check if an atom is connected to a methyl group

def methyl_group(mol, single_carbon_index):
  # Check if the atom is a carbon atom
  alpha_C = mol.GetAtomWithIdx(single_carbon_index)
  for neighbor in alpha_C.GetNeighbors():
    if neighbor.GetSymbol() == 'C':
      beta_h = sum(1 for sub_neigh in neighbor.GetNeighbors() if sub_neigh.GetAtomicNum() == 1)
     # This rule only applies when there is one hydrogen on beta-carbon
      if beta_h == 1:
        # Get the number of hydrogen atoms connected to the atom
       for sub_neighbor in neighbor.GetNeighbors():
          if sub_neighbor.GetSymbol() == 'C':
              num_h = sum(1 for neighbor in sub_neighbor.GetNeighbors() if neighbor.GetAtomicNum() == 1)
        # Check if the number of hydrogen atoms is 3
              if num_h == 3:
        # Return True if both conditions are met
                return True
  # Return False otherwise
  return False

def is_all_zero(lst):
  # Loop over the numbers in the list
  for num in lst:
    # Check if the number is not zero
    if num != 0:
      # Return False if any number is not zero
      return False
  # Return True if all numbers are zero
  return True


# combine the previous rules
def combine_all_rules_together(input_smiles):
    score = None
    mol =  Chem.MolFromSmiles(input_smiles)
    mol = Chem.AddHs(mol)
    #alpha-hydrogen score dict
    alpha_Hydrogen_Score ={'[0, 2]':3,'[0, 3]':2,'[1, 2]':3,'[1, 3]':3,'[2, 2]':1,'[2, 3]':1}
    # #initial_group_number
    # group_number = 5
    #generate the Nitrosamine-pattern related index N and alpha carbons
    pattern_index = is_nitro_pattern(mol)
    if not pattern_index:
        message = []
        message.append('No nitrosamine pattern is found in the molecule')
        return message, score
    elif len(pattern_index) == 1:
        N_index = pattern_index[0][0]
        message = []
        alpha_carbon = pattern_index[0][1]
        hydrogen_pattern = hydrogen_num(mol, alpha_carbon)
        #greater than zero alpha-hydrogen counts
        greater_zero_hydrogen_counts = 0
        #greater than one alpha-hydrogen counts
        greater_one_hydrogen_counts = 0
        for h_num in hydrogen_pattern:
                if h_num > 0:
                    greater_zero_hydrogen_counts += 1
                if h_num > 1:
                    greater_one_hydrogen_counts += 1
        if  greater_zero_hydrogen_counts == 0:
            message.append('No alpha hydrogen on its alpha carbon')
            return message,score
        elif greater_one_hydrogen_counts == 0:
            message.append('Do not have more than one alpha hydrogen on its alpha carbons')
            return message,score
        elif  greater_one_hydrogen_counts > 0:
            if tertiary_carbon(mol,alpha_carbon):
                message.append('Have a tertiary alpha carbon')
                return message,score
            else:
                message.append('Calculate Potency Score')
                for i in alpha_Hydrogen_Score.keys():
                    # print(str(sorted(hydrogen_pattern)),i)
                    if str(sorted(hydrogen_pattern)) == i:
                        if i == '[0, 2]':
                            for carbon_idx in alpha_carbon:
                                carbon = mol.GetAtomWithIdx(carbon_idx)
                                carbon_h = sum(1 for neighbor in carbon.GetNeighbors() if neighbor.GetAtomicNum() == 1)
                                if carbon_h == 2:
                                    for sub_carb in carbon.GetNeighbors():
                                        if sub_carb.GetSymbol() == 'C':
                                            if sum(1 for sub_h in sub_carb.GetNeighbors() if sub_h.GetAtomicNum() == 1) == 3:
                                                score = 2
                                                result = f"Count of hydrogen atoms on each α-carbon {i} and the methylene α-carbon is part of an ethyl group and corresponding α-Hydrogen Score is {score}"
                                            else:
                                                score = alpha_Hydrogen_Score[i]
                                                result = f"Count of hydrogen atoms on each α-carbon {i} and corresponding α-Hydrogen Score is {score}"
                        else:
                            score = alpha_Hydrogen_Score[i]
                            result = f"Count of hydrogen atoms on each α-carbon {i} and corresponding α-Hydrogen Score is {score}"
                        message.append(result)
                # Deactivating feature
                if caboxylic_acid(mol):
                    score += 3
                    message.append('Carboxylic acid group is found, Individual Deactivating Feature Score +3')
                chain = []
                if get_ring_related_score(mol,N_index)>0:
                    score += get_ring_related_score(mol,N_index)
                    message.append(f'Ring related is found, Individual Deactivating Feature Score is +{get_ring_related_score(mol,N_index)}')
                if five_consec_atoms(mol,N_index,alpha_carbon,chain) >0:
                    score += five_consec_atoms(mol,N_index,alpha_carbon,chain)
                    message.append('Chains of ≥5 consecutive non-hydrogen atoms (cyclic or acyclic) on both side,Individual Deactivating Feature Score +1')
                for k in alpha_carbon:
                    if identify_ewg_alpha(mol,N_index,k)>0:
                        score += identify_ewg_alpha(mol,N_index,k)

                        message.append('Electron-withdrawing group** bonded to α-carbon,Individual Deactivating Feature Score +1')
                    # print(hydroxyl_beta_carbon(mol,k))
                    if hydroxyl_beta_carbon(mol,k)>0:
                        score += hydroxyl_beta_carbon(mol,k)
                        message.append('Hydroxyl group bonded to β-carbon,Individual Deactivating Feature Score +1')
                for i in alpha_carbon:
                    if aryl_group_check(mol,i):
                        score -= 1
                        message.append('Aryl group bonded to α-carbon, Activating score -1')
                        break
                for i in alpha_carbon:
                    # print(i)
                    # print(methyl_group(mol,i))
                    if methyl_group(mol,i):
                        score -= 1
                        message.append('Methyl group bonded to beta-carbon, Activating score -1')
                        break
                return message, score
    else:
            # print('>>2')
            score_list =[]
            message_list =[]
            for sub_pattern in pattern_index:
                message = []
                # print(sub_pattern)
                N_index = sub_pattern[0]
                message = []
                alpha_carbon = sub_pattern[1]
                hydrogen_pattern = hydrogen_num(mol, alpha_carbon)
        #greater than zero alpha-hydrogen counts
                greater_zero_hydrogen_counts = 0
        #greater than one alpha-hydrogen counts
                greater_one_hydrogen_counts = 0
                for h_num in hydrogen_pattern:
                    if h_num > 0:
                        greater_zero_hydrogen_counts += 1
                    if h_num > 1:
                        greater_one_hydrogen_counts += 1
                if  greater_zero_hydrogen_counts == 0:
                    message.append('No alpha hydrogen on its alpha carbon')
                elif greater_one_hydrogen_counts == 0:
                    message.append('Do not have more than one alpha hydrogen on its alpha carbons')
                elif  greater_one_hydrogen_counts > 0:
                    if tertiary_carbon(mol,alpha_carbon):
                        message.append('Have a tertiary alpha carbon')
                    else:
                        message.append('Calculate Potency Score')
                        for i in alpha_Hydrogen_Score.keys():
                            # print(str(sorted(hydrogen_pattern)),i)
                            if str(sorted(hydrogen_pattern)) == i:
                                if i == '[0, 2]':
                                    for carbon_idx in alpha_carbon:
                                        carbon = mol.GetAtomWithIdx(carbon_idx)
                                        carbon_h = sum(1 for neighbor in carbon.GetNeighbors() if neighbor.GetAtomicNum() == 1)
                                        if carbon_h == 2:
                                            for sub_carb in carbon.GetNeighbors():
                                                if sub_carb.GetSymbol() == 'C':
                                                    if sum(1 for sub_h in sub_carb.GetNeighbors() if sub_h.GetAtomicNum() == 1) == 3:
                                                        score = 2
                                                        result = f"Count of hydrogen atoms on each α-carbon {i} and the methylene α-carbon is part of an ethyl group and corresponding α-Hydrogen Score is {score}"
                                                    else:
                                                        score = alpha_Hydrogen_Score[i]
                                                        result = f"Count of hydrogen atoms on each α-carbon {i} and corresponding α-Hydrogen Score is {score}"
                                else:
                                    score = alpha_Hydrogen_Score[i]
                                    result = f"Count of hydrogen atoms on each α-carbon {i} and corresponding α-Hydrogen Score is {score}"
                                message.append(result)
                # Deactivating feature
                        if caboxylic_acid(mol):
                            score += 3
                            message.append('Carboxylic acid group is found, Individual Deactivating Feature Score +3')
                        chain = []
                        if get_ring_related_score(mol,N_index)>0:
                            score += get_ring_related_score(mol,N_index)
                            message.append(f'Ring related is found, Individual Deactivating Feature Score is +{get_ring_related_score(mol,N_index)}')
                        if five_consec_atoms(mol,N_index,alpha_carbon,chain) >0:
                                score += five_consec_atoms(mol,N_index,alpha_carbon,chain)
                                message.append('Chains of ≥5 consecutive non-hydrogen atoms (cyclic or acyclic) on both side,Individual Deactivating Feature Score +1')
                        for k in alpha_carbon:
                            if identify_ewg_alpha(mol,N_index,k)>0:
                                score += identify_ewg_alpha(mol,N_index,k)
                                message.append('Electron-withdrawing group** bonded to α-carbon,Individual Deactivating Feature Score +1')
                            if hydroxyl_beta_carbon(mol,k)>0:
                                score += hydroxyl_beta_carbon(mol,k)
                                message.append('Hydroxyl group bonded to β-carbon,Individual Deactivating Feature Score +1')
                        for i in alpha_carbon:
                            if aryl_group_check(mol,i):
                                score -= 1
                                message.append('Aryl group bonded to α-carbon, Activating score -1')
                                break
                        for i in alpha_carbon:
                            if methyl_group(mol,i):
                                score -= 1
                                message.append('Methyl group bonded to beta-carbon, Activating score -1')
                                break
                message_list.append(message)
                if not score:
                    score_list.append(100)
                else:
                    score_list.append(score)
            final_score = min(score_list)
            min_index = score_list.index(final_score)
            final_message = message_list[min_index]
            return final_message, final_score

    return message,score
 
def score_to_category(score):
    if score == None or score == 100:
        category = 5
    elif score >= 4:
        category = 4
    elif score == 3:
        # category = 'Potency Category 3 : AI = 400 ng/day'
        category = 3
    elif score == 2:
        category = 2
    elif score <= 1:
        # category = 'Potency Category 1 : AI = 18 ng/day'
        category = 1
    return category
