#iFerment156

from __future__ import print_function
from cobra import Model, Reaction, Metabolite
import pandas
import numpy
from cobra.flux_analysis import flux_variability_analysis
import cobra.test
import math

################################
###Part I: REACTOR PARAMETERS###
################################

model = Model("iFerment156")

#Set reactor conditions
T=35 #deg C
pH_out = 5.5 #s.u.

##Consider Transport Energetics

TransEnergetics = True

##Set Extracellular Concentrations of unprotonated forms (M)
S_Ethanol = 0.0183
S_Lactate = 0.000001
S_Formate = 0.0047
S_Acetate = 0.1073
S_Propionate = 0.00001
S_Butyrate = 0.0802
S_Valerate = 0.0019
S_Hexanoate = 0.0452
S_Heptanoate = 0.000001
S_Octanoate = 0.0032


##Set Intracellular Concentrations (M)
C_in_Formate = 0.001
C_in_Acetate = 0.001
C_in_Propionate = 0.001
C_in_Butyrate = 0.001
C_in_Valerate = 0.001
C_in_Hexanoate = 0.001
C_in_Heptanoate = 0.001
C_in_Octanoate = 0.001
C_in_Lactate = 0.001
C_in_Ethanol = 0.001

pH_in = 7

R = 8.3145*10**-3 #kj/mol-K

deltaG_ATP_Hydrolysis = -50

h_j = 1 #Number of protons translocated for all acid products
deltaG_pH = -2.3*h_j*R*T*(pH_out-pH_in)
print ('dG_pH: ', deltaG_pH)

delta_Sai = 33.33*(pH_in-pH_out)-143.33
print ('dSai: ', delta_Sai) # = mV
c_j = -1 #For Protons- Net charge transported from outside to inside the cell
F= .096485 # kJ/mV*mol

deltaG_Sai = 1*delta_Sai*c_j*F
print ('dG_Sai: ', deltaG_Sai)


print("Reactions: " + str(len(model.reactions)))
print("Metabolites: " + str(len(model.metabolites)))
print("Genes: " + str(len(model.genes)))


################################
###Part II: BUILD THE MODEL#####
################################

#Complex Carbohydrate Degradation

#Xylan Degradation

xyl4_e = Metabolite('xyl4_e',formula='C20H34O17',name='xyl4_e',compartment='e')

xyl4_c = Metabolite('xyl4_c',formula='C20H34O17',name='xyl4_c',compartment='c')

xyl__D_e = Metabolite('xyl__D_e',formula='C5H10O5',name='xylose-D',compartment='e')

xyl__D_c = Metabolite('xyl__D_c', formula='C5H10O5', name='xylose-D', compartment='c', charge=0)

h2o_c = Metabolite('h2o_c', formula='H2O', name='H2O', compartment='c', charge=0)

reaction = Reaction('EX_xyl4_e')
reaction.name = 'Xylan Exchange'
reaction.subsystem = 'Complex Carbohydrate Degradation'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({xyl4_e: -1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('xyl4t')
reaction.name = 'Xylan Transport'
reaction.subsystem = 'Complex Carbohydrate Degradation'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({xyl4_e: -1.0,
                          xyl4_c: 1.0})


model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('C5Hyd')
reaction.name = 'Xylan Hydrolysis'
reaction.subsystem = 'Complex Carbohydrate Degradation'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({xyl4_c: -1.0,
                          h2o_c: -3.0,
                          xyl__D_e: 4.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#Glucan Degradation

glc4_e = Metabolite('glc4_e',formula='C24H42O21',name='glucan',compartment='e')

glc4_c = Metabolite('glc4_c',formula='C24H42O21',name='glucan',compartment='c')

glc__D_e = Metabolite('glc__D_e',formula='C6H12O6',name='D-Glucose',compartment='e', charge=0)

glc__D_c = Metabolite('glc__D_c', formula='C6H12O6',name='D-Glucose',compartment='c', charge=0)

reaction = Reaction('EX_glc4_e')
reaction.name = 'Glucan Exchange'
reaction.subsystem = 'Complex Carbohydrate Degradation'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

model.add_reactions([reaction])

reaction.add_metabolites({glc4_e: -1.0})

print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('glc4t')
reaction.name = 'Glucan Transport'
reaction.subsystem = 'Complex Carbohydrate Degradation'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({glc4_e: -1.0,
                          glc4_c: 1.0})


model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('C6Hyd')
reaction.name = 'Glucan Hydrolysis'
reaction.subsystem = 'Complex Carbohydrate Degradation'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({glc4_c: -1.0,
                          h2o_c: -3.0,
                          glc__D_e: 4.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

##Pentose Utilization

reaction = Reaction('EX_xyl__D_e')
reaction.name = 'D-Xylose exchange'
reaction.subsystem = 'Transport'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({xyl__D_e: -1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#xyl__D_e <-> xyl__D_c

reaction = Reaction('XYLt')
reaction.name = 'D xylose reversible transport'
reaction.subsystem = 'Pentose Utilization'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({xyl__D_e: -1.0,
                          xyl__D_c: 1.0})

model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

#atp_c + h2o_c + xyl__D_e <-> adp_c + h_c + pi_c + xyl__D_c

atp_c = Metabolite('atp_c', formula='C10H12N5O13P3', name='ATP', compartment='c', charge=-4)
adp_c = Metabolite('adp_c', formula='C10H12N5O10P2', name='ADP', compartment='c', charge=-3)
h_c = Metabolite('h_c', formula='H', name='H+', compartment='c', charge=1)
pi_c = Metabolite('pi_c', formula='HO4P', name='xylose-D', compartment='c', charge=-2)

reaction = Reaction('XYLabc')
reaction.name = 'D-xylose transport via ABC system'
reaction.subsystem = 'Pentose Utilization'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({atp_c: -1.0,
                          h2o_c: -1.0,
                          xyl__D_e: -1.0,
                          adp_c: 1.0,
                          h_c: 1.0,
                          pi_c: 1.0,
                          xyl__D_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#xyl__D_c <-> xylu__D_c

xylu__D_c = Metabolite('xylu__D_c', formula='C5H10O5', name='D-xylulose', compartment='c', charge=0)

reaction = Reaction('XYLI1')
reaction.name = 'Xylose isomerase'
reaction.subsystem = 'Pentose Utilization'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({xyl__D_c: -1.0,
                          xylu__D_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#atp_c + xylu__D_c <-> adp_c + h_c + xu5p__D_c

xu5p__D_c = Metabolite('xu5p__D_c', formula='C5H9O8P', name='D-Xylulose 5-phosphate', compartment='c', charge=-2)

reaction = Reaction('XYLK')
reaction.name = 'Xylulokinase'
reaction.subsystem = 'Pentose Utilization'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({atp_c: -1.0,
                          xylu__D_c: -1.0,
                          adp_c: 1.0,
                          h_c: 1.0,
                          xu5p__D_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

###Phosphoketolase
#pi_c + xu5p__D_c <-> actp_c + g3p_c + h2o_c

actp_c = Metabolite('actp_c', formula='C2H3O5P', name='Acetyl phosphate', compartment='c', charge=-2)
g3p_c = Metabolite('g3p_c', formula='C3H5O6P', name='Glyceraldehyde 3-phosphate', compartment='c', charge=-2)

reaction = Reaction('PKETX')
reaction.name = 'Phosphoketolase (xylulose-5-phosphate utilizing)'
reaction.subsystem = 'Pentose Utilization'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({pi_c: -1.0,
                          xu5p__D_c: -1.0,
                          actp_c: 1.0,
                          g3p_c: 1.0,
                          h2o_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

##Pentose Phosphate Pathway

#ru5p__D_c <-> xu5p__D_c

ru5p__D_c = Metabolite('ru5p__D_c', formula='C5H9O8P', name='D-Ribulose 5-phosphate', compartment='c', charge=-2)

reaction = Reaction('RPE')
reaction.name = 'Ribulose 5-phosphate 3-epimerase'
reaction.subsystem = 'Pentose Utilization'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({ru5p__D_c: -1.0,
                          xu5p__D_c: 1.0})


model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#r5p_c <-> ru5p__D_c

r5p_c = Metabolite('r5p_c', formula='C5H9O8P', name='Alpha-D-Ribose 5-phosphate', compartment='c', charge=-2)

reaction = Reaction('RPI')
reaction.name = 'Ribose-5-phosphate isomerase'
reaction.subsystem = 'Pentose Utilization'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({r5p_c: -1.0,
                          ru5p__D_c: 1.0})


model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#r5p_c + xu5p__D_c <-> g3p_c + s7p_c

s7p_c = Metabolite('s7p_c', formula='C7H13O10P', name='Sedoheptulose 7-phosphate', compartment='c', charge=-2)

reaction = Reaction('TKT1')
reaction.name = 'Transketolase 1'
reaction.subsystem = 'Pentose Utilization'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({r5p_c: -1.0,
                          xu5p__D_c: -1.0,
                          g3p_c: 1.0,
                          s7p_c: 1.0})


model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#g3p_c + s7p_c <-> e4p_c + f6p_c

f6p_c = Metabolite('f6p_c', formula='C6H11O9P', name='D-Fructose 6-phosphate', compartment='c', charge=-2)
e4p_c = Metabolite('e4p_c', formula='C4H7O7P', name='D-Erythrose 4-phosphate', compartment='c', charge=-2)

reaction = Reaction('TALA')
reaction.name = 'Transaldolase'
reaction.subsystem = 'Pentose Utilization'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({g3p_c: -1.0,
                          s7p_c: -1.0,
                          e4p_c: 1.0,
                          f6p_c: 1.0})


model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#e4p_c + xu5p__D_c <-> f6p_c + g3p_c

reaction = Reaction('TKT2')
reaction.name = 'Transketolase 2'
reaction.subsystem = 'Pentose Utilization'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({e4p_c: -1.0,
                          xu5p__D_c: -1.0,
                          f6p_c: 1.0,
                          g3p_c: 1.0})


model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#Glucose Utilization

reaction = Reaction('EX_glc__D_e')
reaction.name = 'D-Glucose exchange'
reaction.subsystem = 'Transport'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({glc__D_e: -1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#glc__D_e <-> glc__D_c

reaction = Reaction('GLCt')
reaction.name = 'D glucose reversible transport'
reaction.subsystem = 'Hexose Utilization'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({glc__D_e: -1.0,
                          glc__D_c: 1.0})

model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

#atp_c + h2o_c + glc__D_e <-> adp_c + h_c + pi_c + glc__D_c


reaction = Reaction('GLCabc')
reaction.name = 'D-glucose transport via ABC system'
reaction.subsystem = 'Hexose Utilization'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({atp_c: -1.0,
                          h2o_c: -1.0,
                          glc__D_e: -1.0,
                          adp_c: 1.0,
                          h_c: 1.0,
                          pi_c: 1.0,
                          glc__D_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#atp_c + glc__D_c <-> adp_c + g6p_c + h_c

g6p_c = Metabolite('g6p_c', formula='C6H11O9P', name='D-Glucose 6-phosphate', compartment='c', charge=-2)

reaction = Reaction('HEX1')
reaction.name = 'Hexokinase (D-glucose:ATP)'
reaction.subsystem = 'Hexose Utilization'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({atp_c: -1.0,
                          glc__D_c: -1.0,
                          adp_c: 1.0,
                          g6p_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#atp_c + glc__D_c <-> adp_c + g6p_c + h_c

reaction = Reaction('PGI')
reaction.name = 'Glucose-6-phosphate isomerase'
reaction.subsystem = 'Hexose Utilization'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({g6p_c: -1.0,
                          f6p_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))


#Phosphoketolase

#f6p_c + pi_c <-> actp_c + e4p_c + h2o_c

reaction = Reaction('PKETF')
reaction.name = 'Phosphoketolase (fructose-6-phosphate utilizing)'
reaction.subsystem = 'Phosphoketolase'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({f6p_c: -1.0,
                          pi_c: -1.0,
                          actp_c: 1.0,
                          e4p_c: 1.0,
                          h2o_c: 1.0})


model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))


##Upper Glycolysis

#atp_c + f6p_c <-> adp_c + fdp_c + h_c

fdp_c = Metabolite('fdp_c', formula='C6H10O12P2', name='D-Fructose 1,6-bisphosphate', compartment='c', charge=-4)

reaction = Reaction('PFK')
reaction.name = 'Phosphofructokinase'
reaction.subsystem = 'Upper Glycolysis'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({atp_c: -1.0,
                          f6p_c: -1.0,
                          adp_c: 1.0,
                          fdp_c: 1.0,
                          h_c: 1.0})


model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#fdp_c <-> dhap_c + g3p_c

dhap_c = Metabolite('dhap_c', formula='C3H5O6P', name='Dihydroxyacetone phosphate', compartment='c', charge=-2)

reaction = Reaction('FBA')
reaction.name = 'Fructose-bisphosphate aldolase'
reaction.subsystem = 'Upper Glycolysis'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({fdp_c: -1.0,
                          dhap_c: 1.0,
                          g3p_c: 1.0})


model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#dhap_c <-> g3p_c

dhap_c = Metabolite('dhap_c', formula='C3H5O6P', name='Dihydroxyacetone phosphate', compartment='c', charge=-2)

reaction = Reaction('TPI')
reaction.name = 'Triose-phosphate isomerase'
reaction.subsystem = 'Upper Glycolysis'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({dhap_c: -1.0,
                          g3p_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

##Lower Glycolysis

#g3p_c + nad_c + pi_c <-> 13dpg_c + h_c + nadh_c

nad_c = Metabolite('nad_c', formula='C21H26N7O14P2', name='Nicotinamide adenine dinucleotide', compartment='c', charge=-1)
nadh_c = Metabolite('nadh_c', formula='C21H27N7O14P2', name='Nicotinamide adenine dinucleotide - reduced', compartment='c', charge=-2)
_13dpg_c = Metabolite('_13dpg_c', formula='C3H4O10P2', name='3-Phospho-D-glyceroyl phosphate', compartment='c', charge=-4)

reaction = Reaction('GAPD')
reaction.name = 'Glyceraldehyde-3-phosphate dehydrogenase'
reaction.subsystem = 'Lower Glycolysis'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({g3p_c: -1.0,
                          nad_c: -1.0,
                          pi_c: -1.0,
                          _13dpg_c: 1.0,
                          h_c: 1.0,
                          nadh_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#3pg_c + atp_c <-> 13dpg_c + adp_c

_3pg_c = Metabolite('_3pg_c', formula='C3H4O7P', name='3-Phospho-D-glycerate', compartment='c', charge=-3)

reaction = Reaction('PGK')
reaction.name = 'Phosphoglycerate kinase'
reaction.subsystem = 'Lower Glycolysis'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({_3pg_c: -1.0,
                          atp_c: -1.0,
                          _13dpg_c: 1.0,
                          adp_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#2pg_c <-> 3pg_c

_2pg_c = Metabolite('_2pg_c', formula='C3H4O7P', name='2-Phospho-D-glycerate', compartment='c', charge=-3)

reaction = Reaction('PGM')
reaction.name = 'Phosphoglycerate mutase'
reaction.subsystem = 'Lower Glycolysis'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({_2pg_c: -1.0,
                          _3pg_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#2pg_c <-> h2o_c + pep_c

pep_c = Metabolite('pep_c', formula='C3H2O6P', name='Phosphoenolpyruvate', compartment='c', charge=-3)

reaction = Reaction('ENO')
reaction.name = 'Enolase'
reaction.subsystem = 'Lower Glycolysis'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({_2pg_c: -1.0,
                          h2o_c: 1.0,
                          pep_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#adp_c + h_c + pep_c <-> atp_c + pyr_c

pyr_c = Metabolite('pyr_c', formula='C3H3O3', name='Pyruvate', compartment='c', charge=-1)

reaction = Reaction('PYK')
reaction.name = 'Pyruvate kinase'
reaction.subsystem = 'Lower Glycolysis'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({adp_c: -1.0,
                          h_c: -1.0,
                          pep_c: -1.0,
                          atp_c: 1.0,
                          pyr_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

##Lactate Metabolism

#lac__D_c + nad_c <-> h_c + nadh_c + pyr_c

lac__D_c = Metabolite('lac__D_c', formula='C3H5O3', name='D-Lactate', compartment='c', charge=-1)

reaction = Reaction('LDH-D')
reaction.name = 'D-lactate dehydrogenase'
reaction.subsystem = 'Lactate metabolism'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({lac__D_c: -1.0,
                          nad_c: -1.0,
                          h_c: 1.0,
                          nadh_c: 1.0,
                          pyr_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#lac__D_e <-> lac__D_c

lac__D_e = Metabolite('lac__D_e', formula='C3H5O3', name='D-Lactate', compartment='e', charge=-1)

reaction = Reaction('LACDt')
reaction.name = 'D-Lactate transport'
reaction.subsystem = 'Transport'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({lac__D_e: -1.0,
                          lac__D_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#lac_D_c + 2 nad_c + fdred_c -> pyr_c + 2 nadh + fdox_c + h_c

fdred_c = Metabolite('fdred_c', formula='Fe8S8X', name='Ferredoxin (reduced) 2[4Fe-4S]', compartment='c', charge= -2)
fdox_c = Metabolite('fdox_c', formula='Fe8S8X', name='Ferredoxin (oxidized) 2[4Fe-4S]', compartment='c', charge= 0)

reaction = Reaction('EB-iLDH-D')
reaction.name = 'Electron bifurcating lactate dehydrogenase'
reaction.subsystem = 'Lactate metabolism'
reaction.lower_bound = 0  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({lac__D_c: -1.0,
                          fdred_c: -1.0,
                          nad_c: -2.0,
                          pyr_c: 1.0,
                          nadh_c: 2.0,
                          fdox_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#lac__D_e <->

reaction = Reaction('EX_lac__D_e')
reaction.name = 'D-Lactate exchange'
reaction.subsystem = 'Transport'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({lac__D_e: -1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

##Formate metabolism

#coa_c + pyr_c <-> accoa_c + for_c

for_c = Metabolite('for_c', formula='CHO2', name='Formate', compartment='c', charge= -1)
accoa_c = Metabolite('accoa_c', formula='C23H34N7O17P3S', name='Acetyl-CoA', compartment='c', charge=-4)
coa_c = Metabolite('coa_c', formula='C21H32N7O16P3S', name='Coenzyme A', compartment='c', charge=-4)

reaction = Reaction('PFL')
reaction.name = 'Pyruvate formate lyase'
reaction.subsystem = 'Transport'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({pyr_c: -1.0,
                          coa_c: -1.0,
                          accoa_c: 1.0,
                          for_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

# for_c <-> for_e

for_e = Metabolite('for_e', formula='CHO2', name='Formate', compartment='e', charge= -1)

#for_e <->

reaction = Reaction('EX_for_e')
reaction.name = 'Formate exchange'
reaction.subsystem = 'Transport'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({for_e: -1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

##Acetate Metabolism

#ac_c + atp_c <-> actp_c + adp_c

ac_c = Metabolite('ac_c', formula='C2H3O2', name='Acetate', compartment='c', charge=-1)

reaction = Reaction('ACKr')
reaction.name = 'Acetate kinase'
reaction.subsystem = 'Acetate Metabolism'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({ac_c: -1.0,
                         atp_c: -1.0,
                         actp_c: 1.0,
                         adp_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#ac_e <-> ac_c

ac_e = Metabolite('ac_e', formula='C2H3O2', name='Acetate', compartment='e', charge=-1)

#ac_e <->

ac_e = Metabolite('ac_e', formula='C2H3O2', name='Acetate', compartment='e', charge=-1)

reaction = Reaction('EX_ac_e')
reaction.name = 'Acetate exchange'
reaction.subsystem = 'Transport'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({ac_e: -1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#accoa_c + pi_c <-> actp_c + coa_c

reaction = Reaction('PTAr')
reaction.name = 'Phosphotransacetylase'
reaction.subsystem = 'Acetate Metabolism'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({accoa_c: -1.0,
                          pi_c: -1.0,
                          actp_c: 1.0,
                          coa_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#accoa_c + h20_c -> ac_c + coa_c + h_c

reaction = Reaction('ACOAH')
reaction.name = 'Acteyl-CoA hydrolase'
reaction.subsystem = 'Acetate Metabolism'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({accoa_c: -1.0,
                          h2o_c: -1.0,
                          ac_c: 1.0,
                          coa_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))


#Pyruvate Oxidation

#coa_c + pyr_c + fdox_c <-> accoa_c + co2_c + fdred_c + h_c

fdred_c = Metabolite('fdred_c', formula='Fe8S8X', name='Ferredoxin (reduced) 2[4Fe-4S]', compartment='c', charge= -2)
fdox_c = Metabolite('fdox_c', formula='Fe8S8X', name='Ferredoxin (oxidized) 2[4Fe-4S]', compartment='c', charge= 0)
co2_c = Metabolite('co2_c', formula='CO2', name='CO2', compartment='c', charge= 0)

reaction = Reaction('PFOR')
#This reaction differs from BiGG database because a different ferredoxin is used and H+ is a product for mass and charge balance
reaction.name = '*Pyruvate flavodoxin oxidoreductase'
reaction.subsystem = 'Pyruvate Oxidation'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({coa_c: -1.0,
                          pyr_c: -1.0,
                          fdox_c: -1.0,
                          accoa_c: 1.0,
                          co2_c: 1.0,
                          fdred_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#coa_c + nad_c + pyr_c <-> accoa_c + co2_c + nadh_c
reaction = Reaction('PDH')

#This reaction differs from BiGG database because a different ferredoxin is used and H+ is a product for mass and charge balance
reaction.name = 'Pyruvate dehdyrogenase'
reaction.subsystem = 'Pyruvate Oxidation'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({coa_c: -1.0,
                          pyr_c: -1.0,
                          nad_c: -1.0,
                          accoa_c: 1.0,
                          co2_c: 1.0,
                          nadh_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

##Reverse Beta Oxidation

#Butyrate Production (Cycle 1)

#2.0 accoa_c <-> aacoa_c + coa_c

aacoa_c = Metabolite('aacoa_c', formula='C25H36N7O18P3S', name='Acetoacetyl-CoA', compartment='c', charge=-4)

reaction = Reaction('ACACT1r')
reaction.name = 'Acetyl-CoA C-acetyltransferase'
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({accoa_c: -2.0,
                          aacoa_c: 1.0,
                          coa_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#aacoa_c + h_c + nadh_c <-> 3hbcoa_c + nad_c

_3hbcoa_c = Metabolite('_3hbcoa_c', formula='C25H38N7O18P3S', name='(S)-3-Hydroxybutanoyl-CoA', compartment='c', charge=-4)

reaction = Reaction('HACD1')
reaction.name = '3-hydroxyacyl-CoA dehydrogenase (acetoacetyl-CoA)'
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({aacoa_c: -1.0,
                          h_c: -1.0,
                          nadh_c: -1.0,
                          _3hbcoa_c: 1.0,
                          nad_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#3hbcoa_c <-> b2coa_c + h2o_c

b2coa_c = Metabolite('b2coa_c', formula='C25H36N7O17P3S', name='Crotonoyl-CoA', compartment='c', charge=-4)

reaction = Reaction('ECOAH1')
reaction.name = '3-hydroxyacyl-CoA dehydratase (3-hydroxybutanoyl-CoA)'
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({_3hbcoa_c: -1.0,
                          b2coa_c: 1.0,
                          h2o_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#b2coa_c + 2 nadh_c + fdox_c <-> btcoa_c + 2 nad_c + fdred_c

btcoa_c = Metabolite('btcoa_c', formula='C25H38N7O17P3S', name='Butanoyl-CoA', compartment='c', charge= -4)

reaction = Reaction('EBACD1')
#BiGG does not have an electron bifurcating acyl-CoA dehydrogenase reaction
reaction.name = '*Electron Bifurcating Acyl-CoA Dehydrogenase (C4)'
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({b2coa_c: -1.0,
                          nadh_c: -2.0,
                          fdox_c: -1.0,
                          btcoa_c: 1.0,
                          nad_c: 2.0,
                          fdred_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#b2coa_c + h_c + nadh_c <-> btcoa_c + nad_c

reaction = Reaction('ACOAD1')
reaction.name = "Acyl-CoA dehydrogenase (butanoyl-CoA)"
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({b2coa_c: -1.0,
                          nadh_c: -1.0,
                          h_c: -1.0,
                          btcoa_c: 1.0,
                          nad_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#btcoa_c + h2o_c <-> but_c + coa_c

but_c = Metabolite('but_c', formula='C4H7O2', name='Butyrate (n-C4:0)', compartment='c', charge= -1)

reaction = Reaction('ACHC4')
#BiGG does not have this specific acyl-CoA hydrolase reaction
reaction.name = '*Acyl-CoA Hydrolase (C4:0)'
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({btcoa_c: -1.0,
                          h2o_c: -1.0,
                          but_c: 1.0,
                          coa_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#btcoa_c + ac_c <-> but_c + accoa_c


reaction = Reaction('CoATC4')
#BiGG does not have this specific CoAT hydrolase reaction
reaction.name = '*CoA Transferase (C4:0-C2:0)'
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({btcoa_c: -1.0,
                          ac_c: -1.0,
                          but_c: 1.0,
                          accoa_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#but_e <-> but_c

but_e = Metabolite('but_e', formula='C4H7O2', name='Butyrate (n-C4:0)', compartment='e', charge= -1)

#but_e <->

reaction = Reaction('EX_but_e')
reaction.name = 'Butyrate exchange'
reaction.subsystem = 'Exchange'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({but_e: -1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#Hexanoate Production (Cycle 2)

# accoa_x + btcoa_x <-> coa_x + 3ohcoa_x

_3ohcoa_c = Metabolite('3ohcoa_c', formula='C27H40N7O18P3S', name='3-Oxohexanoyl-CoA', compartment='c', charge=-4)

reaction = Reaction('ACACT2')
reaction.name = 'Butanoyl-CoA:acetyl-CoA C-butanoyltransferase'
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({accoa_c: -1.0,
                          btcoa_c: -1.0,
                          _3ohcoa_c: 1.0,
                          coa_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#_3ohcoa_c + h_c + nadh_c <-> _3hhcoa_c + nad_c

_3hhcoa_c = Metabolite('_3hhcoa_c', formula='C27H42N7O18P3S', name='(S)-3-Hydroxyhexanoyl-CoA', compartment='c', charge=-4)

reaction = Reaction('HACD2')
reaction.name = '3-hydroxyacyl-CoA dehydrogenase (3-oxohexanoyl-CoA)'
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({_3ohcoa_c: -1.0,
                          h_c: -1.0,
                          nadh_c: -1.0,
                          _3hhcoa_c: 1.0,
                          nad_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#_3hhcoa_c <-> h2o_c + hx2coa_c

hx2coa_c = Metabolite('hx2coa_c', formula='C27H40N7O17P3S', name='Trans-Hex-2-enoyl-CoA', compartment='c', charge=-4)

reaction = Reaction('ECOAH2')
reaction.name = '3-hydroxyacyl-CoA dehydratase (3-hydroxyhexanoyl-CoA)'
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({_3hhcoa_c: -1.0,
                          hx2coa_c: 1.0,
                          h2o_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#hx2coa_c + 2 nadh_c + fdox_c <-> hxcoa_c + 2 nad_c + fdred_c

hxcoa_c = Metabolite('hxcoa_c', formula='C27H42N7O17P3S', name='Hexanoyl-CoA (n-C6:0CoA)', compartment='c', charge= -4)

reaction = Reaction('EBACD2')
#BiGG does not have an electron bifurcating acyl-CoA dehydrogenase reaction
reaction.name = '*Electron Bifurcating Acyl-CoA Dehydrogenase (C6)'
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({hx2coa_c: -1.0,
                          nadh_c: -2.0,
                          fdox_c: -1.0,
                          hxcoa_c: 1.0,
                          nad_c: 2.0,
                          fdred_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#h_c + hx2coa_c + nadh_c <-> hxcoa_c + nad_c

reaction = Reaction('ACOAD2')
reaction.name = "Acyl-CoA dehydrogenase (hexanoyl-CoA)"
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({hx2coa_c: -1.0,
                          nadh_c: -1.0,
                          h_c: -1.0,
                          hxcoa_c: 1.0,
                          nad_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#hxcoa_c + h2o_c <-> but_c + coa_c

hxa_c = Metabolite('hxa_c', formula='C6H11O2', name='Hexanoate (n-C6:0)', compartment='c', charge= -1)

reaction = Reaction('ACH-C6')
#BiGG does not have this specific acyl-CoA hydrolase reaction
reaction.name = '*Acyl-CoA Hydrolase (C6:0)'
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({hxcoa_c: -1.0,
                          h2o_c: -1.0,
                          hxa_c: 1.0,
                          coa_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#hxcoa_c + ac_c <-> hxa_c + accoa_c

reaction = Reaction('CoATC6')
#BiGG does not have this specific acyl-CoA hydrolase reaction
reaction.name = '*CoA Transferase (C4:0-C2:0)'
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({hxcoa_c: -1.0,
                          ac_c: -1.0,
                          hxa_c: 1.0,
                          accoa_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#hex_e <-> but_c

hxa_e = Metabolite('hxa_e', formula='C6H11O2', name='Hexanoate (n-C6:0)', compartment='e', charge= -1)

#hex_e <->

reaction = Reaction('EX_hxa_e')
reaction.name = 'Hexanoate exchange'
reaction.subsystem = 'Exchange'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({hxa_e: -1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))


#Octanoate Production (Cycle 3)

#accoa_c + hxcoa_c <-> coa_c + 3oocoa_c

_3oocoa_c = Metabolite('_3oocoa_c', formula='C29H44N7O18P3S', name='3-Oxooctanoyl-CoA', compartment='c', charge=-4)

reaction = Reaction('ACACT3')
reaction.name = 'Hexanoyl-CoA:acetyl-CoA C-acyltransferase'
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({accoa_c: -1.0,
                          hxcoa_c: -1.0,
                          _3oocoa_c: 1.0,
                          coa_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#_3oocoa_c + h_c + nadh_c <-> _3hocoa_c + nad_c

_3hocoa_c = Metabolite('_3hocoa_c', formula='C29H46N7O18P3S', name='(S)-3-Hydroxyoctanoyl-CoA', compartment='c', charge=-4)

reaction = Reaction('HACD3')
reaction.name = '3-hydroxyacyl-CoA dehydrogenase (3-oxooctanoyl-CoA)'
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({_3oocoa_c: -1.0,
                          h_c: -1.0,
                          nadh_c: -1.0,
                          _3hocoa_c: 1.0,
                          nad_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#_3hocoa_c <-> h2o_c + oc2coa_c

oc2coa_c = Metabolite('oc2coa_c', formula='C29H44N7O17P3S', name='Trans-Oct-2-enoyl-CoA', compartment='c', charge=-4)

reaction = Reaction('ECOAH3')
reaction.name = '3-hydroxyacyl-CoA dehydratase (3-hydroxyoctanoyl-CoA)'
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({_3hocoa_c: -1.0,
                          oc2coa_c: 1.0,
                          h2o_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#oc2coa_c + 2 nadh_c + fdox_c <-> occoa_c + 2 nad_c + fdred_c

occoa_c = Metabolite('occoa_c', formula='C29H46N7O17P3S', name='Octanoyl-CoA (n-C8:0CoA)', compartment='c', charge= -4)

reaction = Reaction('EBACD3')
#BiGG does not have an electron bifurcating acyl-CoA dehydrogenase reaction
reaction.name = '*Electron Bifurcating Acyl-CoA Dehydrogenase (C8)'
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({oc2coa_c: -1.0,
                          nadh_c: -2.0,
                          fdox_c: -1.0,
                          occoa_c: 1.0,
                          nad_c: 2.0,
                          fdred_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#h_c + oc2coa_c + nadh_c <-> occoa_c + nad_c

reaction = Reaction('ACOAD3')
reaction.name = "Acyl-CoA dehydrogenase (octanoyl-CoA)"
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({oc2coa_c: -1.0,
                          nadh_c: -1.0,
                          h_c: -1.0,
                          occoa_c: 1.0,
                          nad_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#occoa_c + h2o_c <-> octa_c + coa_c

octa_c = Metabolite('octa_c', formula='C8H15O2', name='Octanoate (n-C8:0)', compartment='c', charge= -1)

reaction = Reaction('ACH-C8')
#BiGG does not have this specific acyl-CoA hydrolase reaction
reaction.name = '*Acyl-CoA Hydrolase (C8:0)'
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({occoa_c: -1.0,
                          h2o_c: -1.0,
                          octa_c: 1.0,
                          coa_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#occoa_c + ac_c <-> octa_c + accoa_c

reaction = Reaction('CoATC8')
#BiGG does not have this specific acyl-CoA hydrolase reaction
reaction.name = 'CoA Transferase (C8:0-C2:0)'
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({occoa_c: -1.0,
                          ac_c: -1.0,
                          octa_c: 1.0,
                          accoa_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#octa_e <-> octa_c

octa_e = Metabolite('octa_e', formula='C8H15O2', name='Octanoate (n-C8:0)', compartment='e', charge= -1)

#octa_e <->

reaction = Reaction('EX_octa_e')
reaction.name = 'Octanoate exchange'
reaction.subsystem = 'Exchange'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({octa_e: -1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

##Propionate Production

#Acryloyl-CoA Pathway (Metacyc: Pyruvate fermentation to propanoate II)

#lactate CoA transferase
ppa_c = Metabolite('ppa_c', formula='C3H5O2', name='Propionate (n-C3:0)', compartment='c', charge=-1)
ppcoa_c = Metabolite('ppcoa_c', formula='C24H36N7O17P3S', name='Propanoyl-CoA', compartment='c', charge=-4)
laccoa_c = Metabolite('laccoa_c', formula='C24H36N7O18P3S', name='Lactoyl-CoA', compartment='c', charge=-4)

#Reaction not in BIGG database.

reaction = Reaction('LCT')
reaction.name = 'Lactate CoA transferase'
reaction.subsystem = 'Propionate Production'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({lac__D_c: -1.0,
                          ppcoa_c: -1.0,
                          ppa_c: 1.0,
                          laccoa_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#ppa_e <-> ppa_c

ppa_e = Metabolite('ppa_e', formula='C3H5O2', name='Propionate (n-C3:0)', compartment='e', charge=-1)

#ppa_e <->

reaction = Reaction('EX_ppa_e')
reaction.name = 'Propionate exchange'
reaction.subsystem = 'Exchange'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({ppa_e: -1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))


#lactoyl-CoA dehydratase

laccoa_c = Metabolite('laccoa_c', formula='C24H36N7O18P3S', name='Lactoyl-CoA', compartment='c', charge=-4)
pp2coa_c = Metabolite('pp2coa_c', formula='C24H34N7O17P3S', name='Acrylyl-CoA', compartment='c', charge=-4)

reaction = Reaction('LCD')
reaction.name = 'Lactoyl coA dehydratase'
reaction.subsystem = 'Propionate Production'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({laccoa_c: -1.0,
                          pp2coa_c: 1.0,
                          h2o_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))


#acryloyl-CoA reductase

ppcoa_c = Metabolite('ppcoa_c', formula='C24H36N7O17P3S', name='Propanoyl-CoA', compartment='c', charge=-4)

#Reaction not in BIGG Database

reaction = Reaction('ACR')
reaction.name = 'Acryloyl-CoA reductase'
reaction.subsystem = 'Propionate Production'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({pp2coa_c: -1.0,
                          h_c: -1.0,
                          nadh_c: -1.0,
                          nad_c: 1.0,
                          ppcoa_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#propionate CoA transferase

#Reaction not in BIGG Database

reaction = Reaction('PCT')
reaction.name = 'Propionate CoA Transferase'
reaction.subsystem = 'Propionate Production'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({ppcoa_c: -1.0,
                          ac_c: -1.0,
                          accoa_c: 1.0,
                          ppa_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#Methylmalonyl-CoA Pathway (Metacyc: Pyruvate fermentation to propanoate I)"""

#methylmalonyl-CoA carboxyltransferase
mmcoa__S_c = Metabolite('mmcoa__S_c', formula='C25H35N7O19P3S', name='(S)-Methylmalonyl-CoA', compartment='c', charge=-5)
oaa_c = Metabolite('oaa_c', formula='C4H2O5', name='Oxaloacetate', compartment='c', charge=-2)

#Reaction not in BIGG database.

reaction = Reaction('MCC')
reaction.name = 'Methylmalonyl-CoA Carboxyltransferase'
reaction.subsystem = 'Propionate Production'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({mmcoa__S_c: -1.0,
                          pyr_c: -1.0,
                          ppcoa_c: 1.0,
                          oaa_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))


#malate dehydrogenase
mal__L_c = Metabolite('mal__L_c', formula='C4H4O5', name='L-Malate', compartment='c', charge=-2)

reaction = Reaction('MDH')
reaction.name = 'Malate dehydrogenase'
reaction.subsystem = 'Propionate Production'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({oaa_c: -1.0,
                          nadh_c: -1.0,
                          h_c: -1.0,
                          nad_c: 1.0,
                          mal__L_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#fumarase
fum_c = Metabolite('fum_c', formula='C4H2O4', name='Fumarate', compartment='c', charge=-2)

reaction = Reaction('FUM')
reaction.name = 'Fumarase'
reaction.subsystem = 'Propionate Production'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({mal__L_c: -1.0,
                          h2o_c: 1.0,
                          fum_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))


#fumarate reductase NADH
succ_c = Metabolite('succ_c', formula='C4H4O4', name='Succinate', compartment='c', charge=-2)

reaction = Reaction('FRDx')
reaction.name = 'Fumarate Reductase NADH'
reaction.subsystem = 'Propionate Production'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({fum_c: -1.0,
                          nadh_c: -1.0,
                          h_c: -1.0,
                          nad_c: 1.0,
                          succ_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))


#Propanoyl-CoA: succinate CoA-transferase
succoa_c = Metabolite('succoa_c', formula='C25H35N7O19P3S', name='Succinyl-CoA', compartment='c', charge=-5)

reaction = Reaction('PPCSCT')
reaction.name = 'Propanoyl-CoA: succinate CoA-transferase'
reaction.subsystem = 'Propionate Production'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({succ_c: -1.0,
                          ppcoa_c: -1.0,
                          ppa_c: 1.0,
                          succoa_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))


#Methylmalonyl-CoA mutase
mmcoa__R_c = Metabolite('mmcoa__R_c', formula='C25H35N7O19P3S', name='(R)-Methylmalonyl-CoA', compartment='c', charge=-5)

reaction = Reaction('MMM2')
reaction.name = 'Methylmalonyl-CoA mutase'
reaction.subsystem = 'Propionate Production'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({succoa_c: -1.0,
                          mmcoa__R_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))


#Methylmalonyl-CoA epimerase

reaction = Reaction('MME')
reaction.name = 'Methylmalonyl-CoA epimerase'
reaction.subsystem = 'Propionate Production'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({mmcoa__R_c: -1.0,
                          mmcoa__S_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

##Odd-chain reverse beta-oxidation

#Pentanoate production (Cycle 1)

#accoa_x + ppcoa_x <-> coa_x + 3optcoa_x

_3optcoa_c = Metabolite('_3optcoa_c', formula='C26H38N7O18P3S', name='3-Ocopentanoyl-CoA', compartment='c', charge=-4)

reaction = Reaction('VCACT')
reaction.name = 'Acetyl-CoA C-acyltransferase (3-oxovaleryl-CoA)'
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({accoa_c: -1.0,
                          ppcoa_c: -1.0,
                          _3optcoa_c: 1.0,
                          coa_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#_h_c + nadh_c + 3optcoa_c <-> nad_c + 3hptcoa_c

_3hptcoa_c = Metabolite('_3hptcoa_c', formula='C26H40N7O18P3S', name='3-Hydroxypentoyl-CoA', compartment='c', charge=-4)

reaction = Reaction('HVCD')
reaction.name = '3-hydroxyacyl-CoA dehydrogenase (3-hydroxyacyl-CoA dehydrogenase (3-oxovaleryl-CoA))'
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({_3optcoa_c: -1.0,
                          h_c: -1.0,
                          nadh_c: -1.0,
                          _3hptcoa_c: 1.0,
                          nad_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#_3hptcoa_c <-> h2o_c + pt2coa_c

pt2coa_c = Metabolite('pt2coa_c', formula='C26H38N7O17P3S', name='Pent-2-enoyl-CoA', compartment='c', charge=-4)

reaction = Reaction('VECOAH')
reaction.name = '3-hydroxyacyl-CoA dehydratase (3-hydroxypentanoyl-CoA)'
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({_3hptcoa_c: -1.0,
                          pt2coa_c: 1.0,
                          h2o_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#pt2coa_c + 2 nadh_c + fdox_c <-> ptcoa_c + 2 nad_c + fdred_c

ptcoa_c = Metabolite('ptcoa_c', formula='C26H40N7O17P3S', name='Pentanoyl-CoA', compartment='c', charge=-4)

reaction = Reaction('EBVCD')
#BiGG does not have an electron bifurcating acyl-CoA dehydrogenase reaction
reaction.name = '*Electron Bifurcating Acyl-CoA dehydrogenase (pentanoyl-CoA)'
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({pt2coa_c: -1.0,
                          nadh_c: -2.0,
                          fdox_c: -1.0,
                          ptcoa_c: 1.0,
                          nad_c: 2.0,
                          fdred_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#h_c + pt2coa_c + nadh_c <-> ptcoa_c + nad_c

reaction = Reaction('VCOAD')
reaction.name = "Acyl-CoA dehydrogenase (pentanoyl-CoA)"
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({pt2coa_c: -1.0,
                          nadh_c: -1.0,
                          h_c: -1.0,
                          ptcoa_c: 1.0,
                          nad_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#ptcoa_c + h2o_c <-> pta_c + coa_c

pta_c = Metabolite('pta_c', formula='C5H9O2', name='Pentanoate', compartment='c', charge= -1)

reaction = Reaction('ACHC5')
#BiGG does not have this specific acyl-CoA hydrolase reaction
reaction.name = '*Acyl-CoA Hydrolase (C5:0)'
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({ptcoa_c: -1.0,
                          h2o_c: -1.0,
                          pta_c: 1.0,
                          coa_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#ptcoa_c + ac_c <-> pta_c + accoa_c

reaction = Reaction('CoATC5')
#BiGG does not have this specific acyl-CoA hydrolase reaction
reaction.name = '*CoA Transferase (C5:0-C2:0)'
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({ptcoa_c: -1.0,
                          ac_c: -1.0,
                          pta_c: 1.0,
                          accoa_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#pta_e <-> pta_c
pta_e = Metabolite('pta_e', formula='C5H9O2', name='Pentanoate', compartment='e', charge= -1)

#pta_e <->

reaction = Reaction('EX_pta_e')
reaction.name = 'Pentanoate exchange'
reaction.subsystem = 'Exchange'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({pta_e: -1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#Heptanoate production (Cycle 2)

#accoa_x + ppcoa_x <-> coa_x + 3optcoa_x

#3-Oxoheptanoyl-CoA is only in BiGG as M00877. Will define as 3ohtcoa_c

_3ohtcoa_c = Metabolite('_3ohtcoa_c', formula='C28H42N7O18P3S', name='3-Oxoheptanoyl-CoA', compartment='c', charge=-4)

#Reaction not in BiGG Database
reaction = Reaction('VCACT2')
reaction.name = 'Acetyl-CoA C-acyltransferase (3-oxoheptanoyl-CoA)'
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({accoa_c: -1.0,
                          ptcoa_c: -1.0,
                          _3ohtcoa_c: 1.0,
                          coa_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#_h_c + nadh_c + 3ohtcoa_c <-> nad_c + 3hptcoa_c


_3hhtcoa_c = Metabolite('_3hhtcoa_c', formula='C28H44N7O18P3S', name='3-Hydroxyheptanoyl-CoA', compartment='c', charge=-4)

#Reaction is not in BiGG Database
reaction = Reaction('HVCD2')
reaction.name = '3-hydroxyacyl-CoA dehydrogenase (3-oxoheptanoyl-CoA)'
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({_3ohtcoa_c: -1.0,
                          h_c: -1.0,
                          nadh_c: -1.0,
                          _3hhtcoa_c: 1.0,
                          nad_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#_3hhtcoa_c <-> h2o_c + ht2coa_c

ht2coa_c = Metabolite('ht2coa_c', formula='C28H42N7O17P3S', name='Hept-2-enoyl-CoA', compartment='c', charge=-4)

reaction = Reaction('VECOAH2')
reaction.name = '3-hydroxyacyl-CoA dehydratase (3-hydroxyheptanoyl-CoA)'
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({_3hhtcoa_c: -1.0,
                          ht2coa_c: 1.0,
                          h2o_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#ht2coa_c + 2 nadh_c + fdox_c <-> hptcoa_c + 2 nad_c + fdred_c

hptcoa_c = Metabolite('hptcoa_c', formula='C28H44N7O17P3S', name='Heptanoyl-CoA', compartment='c', charge=-4)

reaction = Reaction('EBVCD2')
#BiGG does not have an electron bifurcating acyl-CoA dehydrogenase reaction
reaction.name = '*Electron Bifurcating Acyl-CoA dehydrogenase (heptanoyl-CoA)'
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({ht2coa_c: -1.0,
                          nadh_c: -2.0,
                          fdox_c: -1.0,
                          hptcoa_c: 1.0,
                          nad_c: 2.0,
                          fdred_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#h_c + pt2coa_c + nadh_c <-> ptcoa_c + nad_c

reaction = Reaction('VCOAD2')
reaction.name = "Acyl-CoA dehydrogenase (heptanoyl-CoA)"
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({ht2coa_c: -1.0,
                          nadh_c: -1.0,
                          h_c: -1.0,
                          hptcoa_c: 1.0,
                          nad_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#ptcoa_c + h2o_c <-> pta_c + coa_c

hpta_c = Metabolite('hpta_c', formula='C7H13O2', name='Pentanoate', compartment='c', charge= -1)

reaction = Reaction('ACH-C7')
#BiGG does not have this specific acyl-CoA hydrolase reaction
reaction.name = '*Acyl-CoA Hydrolase (C7:0)'
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({hptcoa_c: -1.0,
                          h2o_c: -1.0,
                          hpta_c: 1.0,
                          coa_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#hptcoa_c + ac_c <-> hpta_c + accoa_c

reaction = Reaction('CoATC7')
#BiGG does not have this specific acyl-CoA hydrolase reaction
reaction.name = '*CoA Transferase (C7:0-C2:0)'
reaction.subsystem = 'Reverse Beta Oxidation'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({hptcoa_c: -1.0,
                          ac_c: -1.0,
                          hpta_c: 1.0,
                          accoa_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#hpta_e <-> hpta_c

hpta_e = Metabolite('hpta_e', formula='C7H13O2', name='Heptanoate', compartment='e', charge= -1)

#hpta_e <->

reaction = Reaction('EX_hpta_e')
reaction.name = 'Heptanoate exchange'
reaction.subsystem = 'Exchange'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({hpta_e: -1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

##Ethanol Utilization

etoh_e = Metabolite('etoh_e',formula='C2H6O', name='Ethanol',compartment='e')

reaction = Reaction('EX_etoh_e')
reaction.name = 'Ethanol exchange'
reaction.subsystem = 'Transport'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({etoh_e: -1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#etoh_e <-> etoh_c

etoh_c = Metabolite('etoh_c',formula='C2H6O', name='Ethanol',compartment='c')

#etoh_c + nad_c <-> acald_c + h_c + nadh_c

acald_c = Metabolite('acald_c',formula='C2H4O', name='Acetaldehyde',compartment='c')

reaction = Reaction('ALCD2x')
reaction.name = 'Alcohol dehydrogenase (ethanol)'
reaction.subsystem = 'Ethanol Utilization'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({etoh_c: -1.0,
                          nad_c: -1.0,
                          acald_c: 1.0,
                          h_c: 1.0,
                          nadh_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#acald_c + coa_c + nad_c <-> accoa_c + h_c + nadh_c

reaction = Reaction('ACALD')
reaction.name = 'Acetaldehyde dehydrogenase (acetylating)'
reaction.subsystem = 'Ethanol Utilization'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({acald_c: -1.0,
                          coa_c: -1.0,
                          nad_c: -1.0,
                          accoa_c: 1.0,
                          h_c: 1.0,
                          nadh_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#Glycerol Utilization"""

glyc_e = Metabolite('glyc_e', formula='C3H8O3', name='Glycerol', compartment='e', charge= 0)
glyc_c = Metabolite('glyc_c', formula='C3H8O3', name='Glycerol', compartment='c', charge= 0)

reaction = Reaction('EX_glyc_e')
reaction.name = 'Glycerol Exchange'
reaction.subsystem = 'Transport'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({glyc_e: -1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('glyct')
reaction.name = 'Glycerol transport'
reaction.subsystem = 'Transport'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({glyc_e: -1.0,
                          glyc_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#atp_c + glyc_c <-> adp_c + glyc3p_c + h_c

glyc3p_c = Metabolite('glyc3p_c', formula='C3H7O6P', name='Glycerol 3-phosphate', compartment='c', charge= -2)

reaction = Reaction('GLYK')
reaction.name = 'Glycerol kinase'
reaction.subsystem = 'Glycerol utilization'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({glyc_c: -1.0,
                          atp_c: -1.0,
                          adp_c: 1.0,
                          glyc3p_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])


print(reaction.name + ": " + str(reaction.check_mass_balance()))

#dhap_c + h_c + nadh_c <-> glyc3p_c + nad_c

reaction = Reaction('G3PD1')
reaction.name = 'Glycerol-3-phosphate dehydrogenase (NAD)'
reaction.subsystem = 'Glycerol utilization'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({dhap_c: -1.0,
                          h_c: -1.0,
                          nadh_c: -1.0,
                          glyc3p_c: 1.0,
                          nad_c: 1.0})

model.add_reactions([reaction])


print(reaction.name + ": " + str(reaction.check_mass_balance()))

##Hydrogen Generation

h2_c = Metabolite('h2_c', formula='H2', name='Hydrogen', compartment='c', charge= 0)

#fdred_c + 2.0 h_c <->  h2_c + fdox_c

reaction = Reaction('HYD1')
#The reaction in BiGG uses a different ferredoxin
#BiGG reaction is not balanced for H
reaction.name = '(FeFe)-hydrogenase, cytoplasm'
reaction.subsystem = 'Hydrogen Generation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({fdred_c: -1.0,
                          h_c: -2.0,
                          h2_c: 1.0,
                          fdox_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#fdred_c + 3.0 h_c <->  h2_c + fdox_c + h_i_

h_i = Metabolite('h_i', formula='H', name='H+', compartment='i', charge=1)

reaction = Reaction('ECH')
#The reaction in BiGG uses a different ferredoxin
#BiGG reaction is not balanced for H
reaction.name = 'Energy conserving hydrogenase'
reaction.subsystem = 'Hydrogen Generation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({fdred_c: -1.0,
                          h_c: -3.0,
                          h2_c: 1.0,
                          fdox_c: 1.0,
                          h_i: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#fdred_c + nadh_c + 4.0 h_c <-> 2 h2_c + fdox_c + nad_c

reaction = Reaction('HYDABC')
#The reaction in BiGG uses a different ferredoxin
#BiGG reaction is not balanced for H
reaction.name = 'Electron confurcating hydrogenase'
reaction.subsystem = 'Hydrogen Generation'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({fdred_c: -1.0,
                          nadh_c: -1.0,
                          h_c: -3.0,
                          h2_c: 2.0,
                          fdox_c: 1.0,
                          nad_c: 1.0})

model.add_reactions([reaction])
#Adding this reaction with the ferredoxin hydrogenase reaction creates a loop in the model

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#h2_e <-> h2_c

h2_e = Metabolite('h2_e', formula='H2', name='Hydrogen', compartment='e', charge= 0)

reaction = Reaction('H2t')
reaction.name = 'Hydrogen transport'
reaction.subsystem = 'Transport'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({h2_e: -1.0,
                          h2_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#h2_e <->

reaction = Reaction('EX_h2_e')
reaction.name = 'H2 exchange'
reaction.subsystem = 'Exchange'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({h2_e: -1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

###Homoacetogensis

#for_c + nad_c <-> co2_c + nadh_c

reaction = Reaction('FDH')
reaction.name = 'Formate dehydrogenase'
reaction.subsystem = 'Wood Ljungadhl Pathway'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({for_c: -1.0,
                          nad_c: -1.0,
                          co2_c: 1.0,
                          nadh_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#atp_c + for_c + thf_c -> 10fthf_c + adp_c + pi_c

thf_c = Metabolite('thf_c', formula='C19H21N7O6', name='5,6,7,8-Tetrahydrofolate', compartment='c', charge= -2)
_10fthf_c = Metabolite('_10fthf_c', formula='C20H21N7O7', name='10-Formyltetrahydrofolate', compartment='c', charge= -2)

reaction = Reaction('FTHFLi')
reaction.name = 'Formate-tetrahydrofolate ligase'
reaction.subsystem = 'Wood Ljungadhl Pathway'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({for_c: -1.0,
                          atp_c: -1.0,
                          thf_c: -1.0,
                          _10fthf_c: 1.0,
                          adp_c: 1.0,
                          pi_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

# 10fthf_c + h_c <-> h2o_c + methf_c

methf_c = Metabolite('methf_c', formula='C20H20N7O6', name='5,10-Methenyltetrahydrofolate', compartment='c', charge= -1)

reaction = Reaction('MTHFC')
reaction.name = 'Methenyltetrahydrofolate cyclohydrolase'
reaction.subsystem = 'Wood Ljungadhl Pathway'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({ _10fthf_c: -1.0,
                          h_c: -1.0,
                          h2o_c: 1.0,
                          methf_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

# methf_c + nadh_c <-> mlthf_c + nad_c

mlthf_c = Metabolite('mlthf_c', formula='C20H21N7O6', name='5,10-Methylenetetrahydrofolate', compartment='c', charge= -2)

reaction = Reaction('MTHFD2i')
reaction.name = 'Methylenetetrahydrofolate dehydrogenase NAD'
reaction.subsystem = 'Wood Ljungadhl Pathway'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({methf_c: -1.0,
                          nadh_c: -1.0,
                          mlthf_c: 1.0,
                          nad_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#2.0 h_c + mlthf_c + nadh_c -> 5mthf_c + nad_c

_5mthf_c = Metabolite('_5mthf_c', formula='C20H24N7O6', name='5-Methyltetrahydrofolate', compartment='c', charge= -1)

reaction = Reaction('MTHFR2')
reaction.name = '5,10-methylenetetrahydrofolate reductase (NADH)'
reaction.subsystem = 'Wood Ljungadhl Pathway'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({mlthf_c: -1.0,
                          nadh_c: -1.0,
                          h_c: -2.0,
                          _5mthf_c: 1.0,
                          nad_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

# 5mthf_c + cfesp_c -> thf_c + mecfsp_c

cfesp_c = Metabolite('cfesp_c', formula='C19CoN4R21', name='Corrinoid Iron sulfur protein', compartment='c', charge=-1)
mecfsp_c = Metabolite('mecfsp_c', formula='C20H3CoN4R21', name='Methylcorrinoid iron sulfur protein', compartment='c', charge=0)

reaction = Reaction('METR')
reaction.name = 'Methyltetrahydrofolate:corrinoid/iron-sulfur protein methyltransferase (MeTr)'
reaction.subsystem = 'Wood Ljungadhl Pathway'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({_5mthf_c: -1.0,
                          cfesp_c: -1.0,
                          thf_c: 1.0,
                          mecfsp_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#co2_c + 2.0 h_c + fdred_c <-> h2o_c + co_c + fdox__c
#BIGG uses a differnt form of ferredoxin

co_c = Metabolite('co_c', formula='CO', name='Carbon monoxide', compartment='c', charge=0)

reaction = Reaction('CODH4')
reaction.name = 'Carbon monoxide dehydrogenase / acetyl-CoA synthase 2'
reaction.subsystem = 'Wood Ljungadhl Pathway'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({co2_c: -1.0,
                          h_c: -2.0,
                          fdred_c: -1.0,
                          h2o_c: 1.0,
                          co_c: 1.0,
                          fdox_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))


#co_c + coa_c  + mecfsp_c -> accoa_c + cfesp_c + h_c

reaction = Reaction('*ACSWL')
reaction.name = 'Acetyl-CoA synthase'
reaction.subsystem = 'Wood Ljungadhl Pathway'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({co_c: -1.0,
                          coa_c: -1.0,
                          mecfsp_c: -1.0,
                          accoa_c: 1.0,
                          cfesp_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

###Energy Generation

#3.0 h_c + nad_c + fdred_c <-> nadh_c + 2.0 h_i + fdox_c


reaction = Reaction('RNF1')
#This reaction differs from the BiGG reaction because a different type of ferredoxin is used.
reaction.name = '*Ferredoxin:NAD oxidoreductase (2 protons translocated)'
reaction.subsystem = 'Energy Generation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({h_c: -3.0,
                          nad_c: -1.0,
                          fdred_c: -1.0,
                          nadh_c: 1.0,
                          h_i: 2.0,
                          fdox_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#adp_c + pi_c + 4.0 h_i <-> atp_c + 3.0 h_c + h2o_c

reaction = Reaction('ATPS4r')
#This reaction differs from the BiGG reaction because this model assumes a different compartment for ion motive force generation
reaction.name = '*ATP Synthase'
reaction.subsystem = 'Energy Generation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({adp_c: -1.0,
                          pi_c: -1.0,
                          h_i: -4.0,
                          atp_c: 1.0,
                          h_c: 3.0,
                          h2o_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

##Other
#h2o_e <-> h2o_c

h2o_e = Metabolite('h2o_e', formula='H2O', name='H2O', compartment='e', charge=0)

reaction = Reaction('H2Ot')
reaction.name = 'H2O transport'
reaction.subsystem = 'Transport'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({h2o_e: -1.0,
                          h2o_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#h2o_e <->

reaction = Reaction('EX_h2o_e')
reaction.name = 'H2O exchange'
reaction.subsystem = 'Transport'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({h2o_e: -1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#h_e <-> h_c

h_e = Metabolite('h_e', formula='H', name='H+', compartment='e', charge=1)

#h_e <->

reaction = Reaction('EX_h_e')
reaction.name = 'H+ exchange'
reaction.subsystem = 'Transport'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({h_e: -1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#co2_e <-> co2_c

co2_e = Metabolite('co2_e', formula='CO2', name='CO2', compartment='e', charge=0)

reaction = Reaction('co2t')
reaction.name = 'CO2 transport'
reaction.subsystem = 'Transport'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({co2_e: -1.0,
                          co2_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#h_e <->

reaction = Reaction('EX_co2_e')
reaction.name = 'CO2 exchange'
reaction.subsystem = 'Exchange'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({co2_e: -1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#ATP Hydrolysis

reaction = Reaction('ATP_Hydrolysis')
reaction.name = 'ATP Hydrolysis'
reaction.subsystem = 'ATP Hydrolysis'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({atp_c: -1.0,
                          h2o_c: -1.0,
                          adp_c: 1.0,
                          pi_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

##Import and Export Reactions For Energy Calculations

#Formate Transport

reaction = Reaction('Formate_import')
reaction.name = 'Formate import'
reaction.subsystem = 'Transport'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({for_e: -1.0,
                          h_e: -1.0,
                          for_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('Formate_export')
reaction.name = 'Formate_export'
reaction.subsystem = 'Transport'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({for_c: -1.0,
                          h_c: -1.0,
                          for_e: 1.0,
                          h_e: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#Acetate Transport

reaction = Reaction('Acetate_import')
reaction.name = 'Acetate import'
reaction.subsystem = 'Transport'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({ac_e: -1.0,
                          h_e: -1.0,
                          ac_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('Acetate_export')
reaction.name = 'Acetate export'
reaction.subsystem = 'Transport'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({ac_c: -1.0,
                          h_c: -1.0,
                          ac_e: 1.0,
                          h_e: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#Propionate Transport

reaction = Reaction('Propionate_import')
reaction.name = 'Propionate import'
reaction.subsystem = 'Transport'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({ppa_e: -1.0,
                          h_e: -1.0,
                          ppa_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('Propionate_export')
reaction.name = 'Propionate export'
reaction.subsystem = 'Transport'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({ppa_c: -1.0,
                          h_c: -1.0,
                          ppa_e: 1.0,
                          h_e: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#Butyrate Transport

reaction = Reaction('Butyrate_import')
reaction.name = 'Butyrate import'
reaction.subsystem = 'Transport'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({but_e: -1.0,
                          h_e: -1.0,
                          but_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('Butyrate_export')
reaction.name = 'Butyrate export'
reaction.subsystem = 'Transport'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({but_c: -1.0,
                          h_c: -1.0,
                          but_e: 1.0,
                          h_e: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#Valerate Transport

reaction = Reaction('Valerate_import')
reaction.name = 'Valerate import'
reaction.subsystem = 'Transport'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({pta_e: -1.0,
                          h_e: -1.0,
                          pta_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('Valerate_export')
reaction.name = 'Valerate export'
reaction.subsystem = 'Transport'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({pta_c: -1.0,
                          h_c: -1.0,
                          pta_e: 1.0,
                          h_e: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#Hexanoate Transport

reaction = Reaction('Hexanoate_import')
reaction.name = 'Hexanoate import'
reaction.subsystem = 'Transport'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({hxa_e: -1.0,
                          h_e: -1.0,
                          hxa_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('Hexanoate_export')
reaction.name = 'Hexanote export'
reaction.subsystem = 'Transport'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({hxa_c: -1.0,
                          h_c: -1.0,
                          hxa_e: 1.0,
                          h_e: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#Heptanoate Transport

reaction = Reaction('Heptanoate_import')
reaction.name = 'Heptanoate import'
reaction.subsystem = 'Transport'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({hpta_e: -1.0,
                          h_e: -1.0,
                          hpta_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('Heptanoate_export')
reaction.name = 'Heptanote export'
reaction.subsystem = 'Transport'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({hpta_c: -1.0,
                          h_c: -1.0,
                          hpta_e: 1.0,
                          h_e: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#Octanoate Transport

reaction = Reaction('Octanoate_import')
reaction.name = 'Octanoate import'
reaction.subsystem = 'Transport'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({octa_e: -1.0,
                          h_e: -1.0,
                          octa_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('Octanoate_export')
reaction.name = 'Octanote export'
reaction.subsystem = 'Transport'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({octa_c: -1.0,
                          h_c: -1.0,
                          octa_e: 1.0,
                          h_e: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#Lactate Transport

reaction = Reaction('Lactate_import')
reaction.name = 'Lactate import'
reaction.subsystem = 'Transport'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 0.  # This is the default

reaction.add_metabolites({lac__D_e: -1.0,
                          h_e: -1.0,
                          lac__D_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('Lactate_export')
reaction.name = 'Lactate export'
reaction.subsystem = 'Transport'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({lac__D_c: -1.0,
                          h_c: -1.0,
                          lac__D_e: 1.0,
                          h_e: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#Ethanol Transport

reaction = Reaction('Ethanol_import')
reaction.name = 'Ethanol import'
reaction.subsystem = 'Transport'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({etoh_e: -1.0,
                          etoh_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('Ethanol_export')
reaction.name = 'Ethanol export'
reaction.subsystem = 'Transport'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({etoh_c: -1.0,
                          etoh_e: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#Proton transport

reaction = Reaction('H_import')
reaction.name = 'H+ import'
reaction.subsystem = 'Transport'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({h_e: -1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('H_export')
reaction.name = 'H+ export'
reaction.subsystem = 'Transport'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({h_c: -1.0,
                          h_e: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#ATP for Transport

#Formate_Transport_ATP
reaction = Reaction('Formate_Transport_ATP')
reaction.name = 'Formate Transport ATP'
reaction.subsystem = 'ATP Hydrolysis'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({atp_c: -1.0,
                          h2o_c: -1.0,
                          adp_c: 1.0,
                          pi_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#Acetate_Transport_ATP

reaction = Reaction('Acetate_Transport_ATP')
reaction.name = 'Acetate Transport ATP Hydrolysis'
reaction.subsystem = 'ATP Hydrolysis'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({atp_c: -1.0,
                          h2o_c: -1.0,
                          adp_c: 1.0,
                          pi_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#Propionate_Transport_ATP

reaction = Reaction('Propionate_Transport_ATP')
reaction.name = 'Propionate Transport ATP Hydrolysis'
reaction.subsystem = 'ATP Hydrolysis'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({atp_c: -1.0,
                          h2o_c: -1.0,
                          adp_c: 1.0,
                          pi_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#Butyrate_Transport_ATP

reaction = Reaction('Butyrate_Transport_ATP')
reaction.name = 'Butyrate Transport ATP Hydrolysis'
reaction.subsystem = 'ATP Hydrolysis'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({atp_c: -1.0,
                          h2o_c: -1.0,
                          adp_c: 1.0,
                          pi_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#Valerate_Transport_ATP

reaction = Reaction('Valerate_Transport_ATP')
reaction.name = 'Valerate Transport ATP Hydrolysis'
reaction.subsystem = 'ATP Hydrolysis'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({atp_c: -1.0,
                          h2o_c: -1.0,
                          adp_c: 1.0,
                          pi_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#Hexanoate_Transport_ATP

reaction = Reaction('Hexanoate_Transport_ATP')
reaction.name = 'Hexanoate Transport ATP Hydrolysis'
reaction.subsystem = 'ATP Hydrolysis'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({atp_c: -1.0,
                          h2o_c: -1.0,
                          adp_c: 1.0,
                          pi_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#Heptanoate_Transport_ATP

reaction = Reaction('Heptanoate_Transport_ATP')
reaction.name = 'Heptanoate Transport ATP Hydrolysis'
reaction.subsystem = 'ATP Hydrolysis'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({atp_c: -1.0,
                          h2o_c: -1.0,
                          adp_c: 1.0,
                          pi_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#Octanoate_Transport_ATP

reaction = Reaction('Octanoate_Transport_ATP')
reaction.name = 'Octanoate Transport ATP Hydrolysis'
reaction.subsystem = 'ATP Hydrolysis'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({atp_c: -1.0,
                          h2o_c: -1.0,
                          adp_c: 1.0,
                          pi_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#Lactate_Transport_ATP

reaction = Reaction('Lactate_Transport_ATP')
reaction.name = 'ATP Transport'
reaction.subsystem = 'ATP Hydrolysis'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({atp_c: -1.0,
                          h2o_c: -1.0,
                          adp_c: 1.0,
                          pi_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#Ethanol_Transport_ATP

reaction = Reaction('Ethanol_Transport_ATP')
reaction.name = 'Ethanol Transport ATP Hydrolysis'
reaction.subsystem = 'ATP Hydrolysis'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({atp_c: -1.0,
                          h2o_c: -1.0,
                          adp_c: 1.0,
                          pi_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

#Proton_Transport_ATP
reaction = Reaction('Proton_Transport_ATP')
reaction.name = 'ATP Transport'
reaction.subsystem = 'ATP Hydrolysis'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({atp_c: -1.0,
                          h2o_c: -1.0,
                          adp_c: 1.0,
                          pi_c: 1.0,
                          h_c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))


#Summarize Model Reactions and Metabolites
print("Reactions: " + str(len(model.reactions)))
print("Metabolites: " + str(len(model.metabolites)))
print("Genes: " + str(len(model.genes)))

## Transport Energy

if TransEnergetics == True:

##Formate Transport Energy
    deltaG_trans_grad_Formate = R*(T+273.15)*(math.log(S_Formate/C_in_Formate))
    ATP_trans_Formate = -1*(deltaG_trans_grad_Formate + deltaG_pH)/ deltaG_ATP_Hydrolysis
    if ATP_trans_Formate > 0:
        Constraint_trans_Formate = model.problem.Constraint(model.reactions.Formate_Transport_ATP.flux_expression - ATP_trans_Formate* model.reactions.Formate_export.flux_expression, lb=0, ub=0)
        model.add_cons_vars(Constraint_trans_Formate)

##Acetate Transport Energy
    deltaG_trans_grad_Acetate = R*(T+273.15)*(math.log(S_Acetate/C_in_Acetate))
    ATP_trans_Acetate = -1*(deltaG_trans_grad_Acetate + deltaG_pH)/ deltaG_ATP_Hydrolysis
    if ATP_trans_Acetate > 0:
        Constraint_trans_Acetate = model.problem.Constraint(model.reactions.Acetate_Transport_ATP.flux_expression - ATP_trans_Acetate * model.reactions.Acetate_export.flux_expression, lb=0, ub=0)
        model.add_cons_vars(Constraint_trans_Acetate)

##Propionate Transport Energy
    deltaG_trans_grad_Propionate = R*(T+273.15)*(math.log(S_Propionate/C_in_Propionate))
    ATP_trans_Propionate = -1*(deltaG_trans_grad_Propionate + deltaG_pH)/ deltaG_ATP_Hydrolysis
    if ATP_trans_Propionate > 0:
        Constraint_trans_Propionate = model.problem.Constraint(model.reactions.Propionate_Transport_ATP.flux_expression - ATP_trans_Propionate* model.reactions.Propionate_export.flux_expression, lb=0, ub=0)
        model.add_cons_vars(Constraint_trans_Propionate)

##Butyrate Transport Energy
    deltaG_trans_grad_Butyrate = R*(T+273.15)*(math.log(S_Butyrate/C_in_Butyrate))
    ATP_trans_Butyrate = -1*(deltaG_trans_grad_Butyrate + deltaG_pH)/ deltaG_ATP_Hydrolysis
    if ATP_trans_Butyrate > 0:
        Constraint_trans_Butyrate = model.problem.Constraint(model.reactions.Butyrate_Transport_ATP.flux_expression - ATP_trans_Butyrate* model.reactions.Butyrate_export.flux_expression, lb=0, ub=0)
        model.add_cons_vars(Constraint_trans_Butyrate)

##Valerate Transport Energy
    deltaG_trans_grad_Valerate = R*(T+273.15)*(math.log(S_Valerate/C_in_Valerate))
    ATP_trans_Valerate = -1*(deltaG_trans_grad_Valerate + deltaG_pH)/ deltaG_ATP_Hydrolysis
    if ATP_trans_Valerate > 0:
        Constraint_trans_Valerate = model.problem.Constraint(model.reactions.Valerate_Transport_ATP.flux_expression - ATP_trans_Valerate* model.reactions.Valerate_export.flux_expression, lb=0, ub=0)
        model.add_cons_vars(Constraint_trans_Valerate)

##Hexanoate Transport Energy
    deltaG_trans_grad_Hexanoate = R*(T+273.15)*(math.log(S_Hexanoate/C_in_Hexanoate))
    ATP_trans_Hexanoate = -1*(deltaG_trans_grad_Hexanoate + deltaG_pH)/ deltaG_ATP_Hydrolysis
    if ATP_trans_Hexanoate > 0:
        Constraint_trans_Hexanoate = model.problem.Constraint(model.reactions.Hexanoate_Transport_ATP.flux_expression - ATP_trans_Hexanoate* model.reactions.Hexanoate_export.flux_expression, lb=0, ub=0)
        model.add_cons_vars(Constraint_trans_Hexanoate)

##Heptanoate Transport Energy
    deltaG_trans_grad_Heptanoate = R*(T+273.15)*(math.log(S_Heptanoate/C_in_Heptanoate))
    ATP_trans_Heptanoate = -1*(deltaG_trans_grad_Heptanoate + deltaG_pH)/ deltaG_ATP_Hydrolysis
    if ATP_trans_Heptanoate > 0:
        Constraint_trans_Heptanoate = model.problem.Constraint(model.reactions.Heptanoate_Transport_ATP.flux_expression - ATP_trans_Heptanoate* model.reactions.Heptanoate_export.flux_expression, lb=0, ub=0)
        model.add_cons_vars(Constraint_trans_Heptanoate)

##Octanoate Transport Energy
    deltaG_trans_grad_Octanoate = R*(T+273.15)*(math.log(S_Octanoate/C_in_Octanoate))
    ATP_trans_Octanoate = -1*(deltaG_trans_grad_Octanoate + deltaG_pH)/ deltaG_ATP_Hydrolysis
    if ATP_trans_Octanoate > 0:
        Constraint_trans_Octanoate = model.problem.Constraint(model.reactions.Octanoate_Transport_ATP.flux_expression - ATP_trans_Octanoate* model.reactions.Octanoate_export.flux_expression, lb=0, ub=0)
        model.add_cons_vars(Constraint_trans_Octanoate)

##Lactate Transport Energy
    deltaG_trans_grad_Lactate = R*(T+273.15)*(math.log(S_Lactate/C_in_Lactate))
    ATP_trans_Lactate = -1*(deltaG_trans_grad_Lactate + deltaG_pH)/ deltaG_ATP_Hydrolysis
    if ATP_trans_Lactate > 0:
        Constraint_trans_Lactate = model.problem.Constraint(model.reactions.Lactate_Transport_ATP.flux_expression - ATP_trans_Lactate* model.reactions.Lactate_export.flux_expression, lb=0, ub=0)
        model.add_cons_vars(Constraint_trans_Lactate)

##Proton Transport Energy
    S_H = 10*math.exp(-pH_out)
    C_in_H = 10*math.exp(-pH_in)
    deltaG_trans_grad_Proton = R*(T+273.15)*(math.log(S_H/C_in_H))
    ATP_trans_Proton = 1*(deltaG_trans_grad_Proton + deltaG_Sai)/ deltaG_ATP_Hydrolysis
    if ATP_trans_Proton > 0:
        Constraint_trans_Proton = model.problem.Constraint(model.reactions.Proton_Transport_ATP.flux_expression - ATP_trans_Proton* model.reactions.H_export.flux_expression, lb=0, ub=0)
        model.add_cons_vars(Constraint_trans_Proton)

##Ethanol Transport Energy
    deltaG_trans_grad_Ethanol = R*(T+273.15)*(math.log(S_Ethanol/C_in_Ethanol))
    ATP_trans_Ethanol = -1*(deltaG_trans_grad_Ethanol)/ deltaG_ATP_Hydrolysis
    if ATP_trans_Ethanol > 0:
        Constraint_trans_Ethanol = model.problem.Constraint(model.reactions.Ethanol_Transport_ATP.flux_expression - ATP_trans_Ethanol* model.reactions.Ethanol_export.flux_expression, lb=0, ub=0)
        model.add_cons_vars(Constraint_trans_Ethanol)


#Summarize Model Reactions and Metabolites
print("Reactions: " + str(len(model.reactions)))
print("Metabolites: " + str(len(model.metabolites)))
print("Genes: " + str(len(model.genes)))

##################################
###Part III: SUBSTRATE UPTAKE#####
##################################

print (model.medium)
medium = model.medium

#Medium for reactor simulations
medium["EX_xyl__D_e"] = 0.1317 #mmol/hr
medium["EX_xyl4_e"] = 0.0081 #mmol/hr
medium["EX_glc4_e"] = 0.0081 #mmol/hr
medium["EX_glc__D_e"] = 0.0125 #mmol/hr
medium["EX_glyc_e"] = 0.0360 #mmol/hr
medium["EX_lac__D_e"] = 0.0005 #mmol/hr
medium["EX_etoh_e"] = 0
medium["EX_ac_e"] = 0
medium["EX_but_e"] = 0
medium["EX_hxa_e"] = 0
medium["EX_h2_e"] = 0
medium["EX_octa_e"] = 0
medium["EX_ppa_e"] = 0
medium["EX_pta_e"] = 0
medium["EX_hpta_e"] = 0
medium["EX_co2_e"] = 0
medium["EX_for_e"] = 0

model.medium = medium
print (model.medium)


#######################################
##PART IV: SET ADDITIONAL CONSTRAINTS##
#######################################

#Turn off end products
#model.reactions.EX_octa_e.knock_out()
#model.reactions.EX_hxa_e.knock_out()
#model.reactions.EX_but_e.knock_out()
#model.reactions.EX_ac_e.knock_out()
#model.reactions.EX_lac__D_e.knock_out()
#model.reactions.EX_h2_e.knock_out()
#model.reactions.EX_etoh_e.knock_out()
#model.reactions.EX_for_e.knock_out()
#model.reactions.EX_ppa_e.knock_out()
#model.reactions.EX_pta_e.knock_out()
#model.reactions.EX_hpta_e.knock_out()

#To allow acetate uptake only
#model.reactions.EX_ac_e.upper_bound = 0
#model.reactions.ACKr.knock_out()
#model.reactions.ACKr.lower_bound = 0

#Turn off electron bifurcating acyl-CoA dehydrogenase
#model.reactions.EBACD1.knock_out()
#model.reactions.EBACD2.knock_out()
#model.reactions.EBACD3.knock_out()
#model.reactions.EBACD3.knock_out()
#model.reactions.EBVCD.knock_out()
#model.reactions.EBVCD2.knock_out()

#Turn off non-electron bifurcating acyl-CoA dehydrogenase
model.reactions.ACOAD1.knock_out()
model.reactions.ACOAD2.knock_out()
model.reactions.ACOAD3.knock_out()
model.reactions.VCOAD.knock_out()
model.reactions.VCOAD2.knock_out()

#Turn off CoaT
#model.reactions.CoATC4.knock_out()
#model.reactions.CoATC6.knock_out()
#model.reactions.CoATC8.knock_out()
#model.reactions.CoATC5.knock_out()
#model.reactions.CoATC7.knock_out()

#Turn of RNF Complex
#model.reactions.RNF1.knock_out()

#Turn of PFOR
#model.reactions.PFOR.knock_out()

#Constrain Hydrogenases to Avoid Loops
if medium["EX_h2_e"] > 0:
    model.reactions.HYD1.upper_bound = 0
    model.reactions.ECH.upper_bound = 0
    model.reactions.HYDABC.upper_bound = 0
else:
    model.reactions.HYD1.lower_bound = 0
    model.reactions.ECH.lower_bound = 0
    model.reactions.HYDABC.lower_bound = 0

#model.reactions.HYD1.knock_out()
#model.reactions.ECH.knock_out()
#model.reactions.HYDABC.knock_out()

#Turn off Phosphoketolase
#model.reactions.PKETX.knock_out()
#model.reactions.RPE.knock_out()
#model.reactions.GLCabc.knock_out()
#model.reactions.PDH.knock_out()
#model.reactions.PFL.upper_bound = 1

#model.reactions.PKETF.knock_out()

#Turn of Odd-chain Production

#Turn off propionate production via acryloyl-CoA pathway
model.reactions.LCD.knock_out()

#Turn off propionate production via methylmalonyl-CoA pathway
model.reactions.MCC.knock_out()

#Turn off homoacetogensis
model.reactions.MCC.knock_out()

#Turn off homoacetogensis
model.reactions.CODH4.knock_out()
model.reactions.FDH.knock_out()

#This is where we set the objective function
model.objective = 'ATP_Hydrolysis'
#model.objective = 'EX_octa_e'
#model.objective = 'EX_lac__D_e'
#model.objective = 'EX_ppa_e'
#model.objective = 'EX_pta_e'
#model.objective = 'EX_hpta_e'
#model.objective = 'EX_ac_e'

#model.reactions.EX_octa_e.upper_bound = 0.0034
#model.reactions.EX_octa_e.lower_bound = 0.0034

#model.reactions.EX_hxa_e.upper_bound = 0.0478
#model.reactions.EX_hxa_e.lower_bound = 0.0478

"""model.reactions.EX_but_e.upper_bound = 0.0827
model.reactions.EX_but_e.lower_bound = 0.0827

model.reactions.EX_ac_e.lower_bound = 0.0740
model.reactions.EX_ac_e.upper_bound = 0.0740

model.reactions.EX_etoh_e.lower_bound = 0.0032
model.reactions.EX_etoh_e.upper_bound = 0.0032

model.reactions.EX_for_e.lower_bound = 0.0001
model.reactions.EX_for_e.upper_bound = 0.0001

#model.reactions.EX_h2_e.lower_bound = 0.437
#model.reactions.EX_h2_e.upper_bound = 0.437

model.reactions.EX_lac__D_e.upper_bound = 0
model.reactions.EX_lac__D_e.lower_bound = 0"""

"""model.reactions.EX_glyc_e.upper_bound = -0.0346708
model.reactions.EX_glyc_e.lower_bound = -0.0346708

model.reactions.EX_glc__D_e.upper_bound = -0.012016617
model.reactions.EX_glc__D_e.lower_bound = -0.012016617

model.reactions.EX_xyl__D_e.upper_bound = -0.1264797
model.reactions.EX_xyl__D_e.lower_bound = -0.1264797

model.reactions.EX_xyl4_e.upper_bound = -0.01994533
model.reactions.EX_xyl4_e.lower_bound = -0.01994533

model.reactions.EX_glc4_e.upper_bound = -0.01994533
model.reactions.EX_glc4_e.lower_bound = -0.01994533"""

#Run pFBA

pfba_solution = cobra.flux_analysis.pfba(model)
model.summary()
print (pfba_solution.fluxes)

#Write pFBA results to excel
writer = pandas.ExcelWriter('output.xlsx')
pfba_solution.fluxes.to_excel(writer,'Sheet1')
writer.save()

model.summary()


#Calculate overall reaction thermodynamics
XYL = pfba_solution["EX_xyl__D_e"]
GLC = pfba_solution["EX_glc__D_e"]
GLYC = pfba_solution["EX_glyc_e"]
LAC = pfba_solution["EX_lac__D_e"]
ETOH = pfba_solution["EX_etoh_e"]
H2 = pfba_solution["EX_h2_e"]
H2O = pfba_solution["EX_h2o_e"]
CO2 = pfba_solution["EX_co2_e"]
H = pfba_solution["EX_h_e"]
C1 = pfba_solution["EX_for_e"]
C2 = pfba_solution["EX_ac_e"]
C4 = pfba_solution["EX_but_e"]
C6 = pfba_solution["EX_hxa_e"]
C8 = pfba_solution["EX_octa_e"]
C3 = pfba_solution["EX_ppa_e"]
C5 = pfba_solution["EX_pta_e"]
C7 = pfba_solution["EX_hpta_e"]

G_XYL = pfba_solution["EX_xyl__D_e"]*-753.37
G_GLC = pfba_solution["EX_glc__D_e"]*-913.28
G_XYL4 = pfba_solution["EX_xyl4_e"]*-547.08*4.184 #Used dG for tetra-arabinofuranoside
G_GLC4 = pfba_solution["EX_glc4_e"]*-694.61*4.184 #Used dG for stachyose
G_GLYC = pfba_solution["EX_glyc_e"]*-116.18*4.184
G_LAC = pfba_solution["EX_lac__D_e"]*-515.34
G_ETOH = pfba_solution["EX_etoh_e"]*-181.75
G_H2 = pfba_solution["EX_h2_e"]*0
G_H2O = pfba_solution["EX_h2o_e"]*-237.18
G_CO2 = pfba_solution["EX_co2_e"]*-386.02
G_H = pfba_solution["EX_h_e"]*0
G_C1 = pfba_solution["EX_for_e"]*-351.0376
G_C2 = pfba_solution["EX_ac_e"]*-369.41
G_C4 = pfba_solution["EX_but_e"]*-352.63
G_C6 = pfba_solution["EX_hxa_e"]*-335.85
G_C8 = pfba_solution["EX_octa_e"]*-322.29
G_C3 = pfba_solution["EX_ppa_e"]*-356.18
G_C5 = pfba_solution["EX_pta_e"]*-342.6
G_C7 = pfba_solution["EX_hpta_e"]*-329.2

dG0 = G_XYL + G_GLC + G_XYL4 + G_GLC4 + G_GLYC + G_LAC + G_ETOH + G_H2 + G_H2O + G_CO2 + G_H + G_C1 + G_C2 + G_C4 + G_C6 + G_C8 + G_C3 + G_C5 + G_C7

print ("dG0: ", dG0)

dG0_prime = dG0 + ((8.3145*10**-3)*298)*numpy.log((10**(-7))**H)

print ("dG0_prime: ", dG0_prime)

G_Per_ATP = dG0_prime/pfba_solution["ATP_Hydrolysis"]

print ("G_Per_ATP: ", G_Per_ATP)

if G_Per_ATP > -50:
    print ("Free Energy per mol ATP < 60 kJ")

rxn_count = 0
for n in pfba_solution.fluxes:
    if n > 0.0000001 or n < -0.0000001:
        rxn_count += 1

print ("Reactions carrying flux: ",rxn_count)

print ("XYL: ",XYL)
print ("GLC: ",GLC)
print ("GLYC: ",GLYC)
print ("LAC: ",LAC)
print ("ETOH: ", ETOH)
print ("H2: ",H2)
print ("H2O: ",H2O)
print ("CO2: ",CO2)
print ("H+: ",H)
print ("C1: ",C1)
print ("C2: ",C2)
print ("C3: ",C3)
print ("C4: ",C4)
print ("C5: ",C5)
print ("C6: ",C6)
print ("C7: ",C7)
print ("C8: ",C8)
print ("ATP: ",pfba_solution["ATP_Hydrolysis"])

if G_Per_ATP > -50:
    print ("Free Energy per mol ATP < 50 kJ")

print (type(model.solver))

sol= model.optimize()
model.summary(fva=1.00)

#Run FVA
fva = flux_variability_analysis(model, loopless=True, fraction_of_optimum=1.00)
print (fva)

#Print FVA results to excel
writer = pandas.ExcelWriter('FVA_iFerment156_Simulation1_NoProdConstr.xlsx')
fva.to_excel(writer,'Sheet1')
writer.save()

#Create .SBML model for use in other modeling platforms

#cobra.io.write_sbml_model(model, "iFerment156.xml")