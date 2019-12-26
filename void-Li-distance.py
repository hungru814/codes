import os
import copy
import numpy as np
import pickle
from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Element
import itertools


metal = ['Li','Be','B','N','Na','Mg','Al','Si','P','S','K','Ca','Sc','Ti','V',\
         'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Rb','Sr','Y','Zr',\
         'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Cs','Ba','La','Ce','Pr','Nd','Pm',\
         'Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Ir','Pt','Au','Pb','Bi']

metal_b = ['Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Y',\
         'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Hf','Ta','W','Re','Os','Ir','Pt','Au','Pb','Bi']

metal_large = ['Na','K','Rb','Cs','Ca','Sr','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu']

unwanted = ['H','He','N','F','NE','S','Cl','Ar','Br','Kr','Te','I','Xe','Hg','Tl',\
            'Ac','Th','Pa','U','Np','Pu']



#################################################################################################
def judge(structure,a):

#    print(structure)

    lattice = structure.lattice.matrix

    structure_orig = copy.deepcopy(structure)

    with open("POSCAR-TEMP-"+m_id,"w") as file:
        file.write('void\n')
        file.write('1.0\n')
        file.write('  '.join(map(str,lattice[0])))
        file.write('\n')
        file.write('  '.join(map(str,lattice[1])))
        file.write('\n')
        file.write('  '.join(map(str,lattice[2])))
        file.write('\n')
        file.write('Li\n')
        file.write('1\n')
        file.write('Direct\n')
        file.write('0 0 0')

    pos_l = Structure.from_file("POSCAR-TEMP-"+m_id)
    pos_l.remove_species("3")

    pos_m = copy.deepcopy(structure)
    pos_m.remove_species("3")
    pos_m.remove_species("8")

#    print(pos_m)

    pos_o = Structure.from_file("POSCAR-TEMP-"+m_id)
    pos_o.remove_species("3")
    for k in range(0,len(structure.sites)):
        if structure.sites[k].specie.symbol == 'O':
            pos_o.append(structure.sites[k].specie.symbol,structure.sites[k].frac_coords)

#    print(pos_o)

    pos_lm = Structure.from_file("POSCAR-TEMP-"+m_id)
    pos_lm.remove_species("3")
    for k in range(0,len(structure.sites)):
        if structure.sites[k].specie.symbol in metal_large:
            pos_lm.append(structure.sites[k].specie.symbol,structure.sites[k].frac_coords)

    dx=int((structure.lattice.a)/0.2)
    dy=int((structure.lattice.b)/0.2)
    dz=int((structure.lattice.c)/0.2)

    radius_O = 1.57
    radius_OO = 2.22
    radius_M = 1.84
    radius_LM = 2.15
#    radius_MM = 2.0
    number = 0

    structure.remove_species("3")

    for k in range(0,dx):
     for l in range(0,dy):
      for m in range(0,dz):
        pt = k*(1/dx)*lattice[0]+l*(1/dy)*lattice[1]+m*(1/dz)*lattice[2]

        if len(pos_lm) == 0:
#            print(pos_o.get_sites_in_sphere(pt,radius_O))
            if len(pos_o.get_sites_in_sphere(pt,radius_O)) == 0 and \
               len(pos_m.get_sites_in_sphere(pt,radius_M)) == 0 and \
               len(pos_o.get_sites_in_sphere(pt,radius_OO)) >= 2 :
#               len(pos_m.get_sites_in_sphere(pt,radius_MM)) < 3 :
                number = number + 1
                structure.append("3",(pt[0],pt[1],pt[2]),coords_are_cartesian=True)
                pos_l.append("3",(pt[0],pt[1],pt[2]),coords_are_cartesian=True)

        if len(pos_lm) > 0:
            if len(pos_o.get_sites_in_sphere(pt,radius_O)) == 0 and \
               len(pos_m.get_sites_in_sphere(pt,radius_M)) == 0 and \
               len(pos_lm.get_sites_in_sphere(pt,radius_LM)) == 0 and \
               len(pos_o.get_sites_in_sphere(pt,radius_OO)) >= 2 :
#               len(pos_m.get_sites_in_sphere(pt,radius_MM)) < 3 :
                number = number + 1
                structure.append("3",(pt[0],pt[1],pt[2]),coords_are_cartesian=True)
                pos_l.append("3",(pt[0],pt[1],pt[2]),coords_are_cartesian=True)

    if number/structure.volume > 14:
        structure.to("poscar",'POSCAR-'+m_id)
        pos_l.to("poscar",'POSCAR-Li-'+m_id)
        structure_orig.to("poscar",'POSCAR-orig-'+m_id)
        with open('entries.dat','at') as fff:
            fff.write(m_id+'   '+formula + '    ' + str(e_above_hull) + '\n')

    os.remove("POSCAR-TEMP-"+m_id)

##################################################################################


with open('/home/chen/MP_data_2019_Mar/mp_structure.pkl.pd_','rb') as f :
    mp_structure = pickle.load(f)


with open('/home/chen/MP_data_2019_Mar/mp_inorganic.pkl.pd_','rb') as g :
    mp_inorganic = pickle.load(g)


list1= list(itertools.combinations(metal,2))

for a in range(0,len(list1)):

    list2=list(list1[a])
    list2.append('O')

    with open(list1[a][0]+'-'+list1[a][1]+'-O',"w") as file:
        file.write('test')

    for i in range(0,len(mp_inorganic)):
        if sorted(mp_inorganic.iloc[i]['elements']) == sorted(list2) \
           and mp_inorganic.iloc[i]['e_above_hull'] < 0.025 :

            e_above_hull = mp_inorganic.iloc[i]['e_above_hull']
            formula = mp_inorganic.iloc[i]['pretty_formula']
            m_id = mp_inorganic.iloc[i].name
            ss = mp_structure.loc[m_id]['structure']
            judge(ss,a)

    os.remove(list1[a][0]+'-'+list1[a][1]+'-O')

