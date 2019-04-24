# Mixed Culture Fermentation Metabolic Models

Copyright (c) 2019, Matthew J. Scarborough
Contact: scarborough.matthew@gmail.com

This readme describes the use of two mixed culture fermentation metabolic models: a unicellular model and a mixed community model

## Prerequisites
* [Python 3.5 or later] (https://www.python.org/downloads/) -  Both metabolic models are written with python. Matlab versions can be created using CobraPy. 
 
* [CobraPy] (https://cobrapy.readthedocs.io/en/latest/) - Python package for metabolic modeling using the COBRA framework.

## Recommended Software
* [PyCharm] (https://www.jetbrains.com/pycharm/) - Integrated development environment used for constructing the models.

## Model Descriptions

**iFerment156**
iFerment156 is a unicellular metabolic model with 156 reactions.

**iFermGuilds548**
iFermGuids548 is a community metbaolic model with 548 reactions distributed amongst six guilds.

## Getting started

**Constraining the reactor environment**
In the model, see **Part I: REACTOR PARAMETERS**

Here, you can tell the model if it should consider transport energetics related to extracellular concentrations of end-products. You can also set the temperature, the biomass concentration (volatile suspended solids, VSS) concentraion, the intracellular and extracellular pH, and the concentrations of products inside and outside the cell. To ignore transport energetics, you can simply set ```TransEnergetics = False```. If you are not interested in modeling a bioreactor, and are interested in flux distributions, you can just set the volume and VSS concentrationeach to ```1``` which will result in all reported fluxes being in units of mmol gDCW<sup>-1</sup> hr<sup>-1</sup>  

**Building the metabolic networks**
In the model, see **Part II: BUILDING THE MODEL**

**Feeding the model**

In the model, see **Part III: SUBSTRATE UPTAKE**

To set substrate uptake, in mmol/hr, set the "medium" value to the desired uptake rate. This is the maximum uptake rate, so if you specify multiple medium componenents, one may be limiting. For all exchange metabolites that you do not want to test,set the value to 0. This will ensure the model does not consume these compounds, but it still allows the model to produce these compounds.  As an example, to test co-utilization of xylose and hydrogen gas with xylose being the limiting substrate, I could set the media constraints as follows:

```
medium = model.medium

medium["EX_xyl__D_e"] = 1 
medium["EX_xyl4_e"] = 0
medium["EX_glc4_e"] = 0
medium["EX_glc__D_e"] = 0
medium["EX_glyc_e"] = 0
medium["EX_lac__D_e"] = 0 
medium["EX_etoh_e"] = 0
medium["EX_ac_e"] = 0
medium["EX_but_e"] = 0
medium["EX_hxa_e"] = 0
medium["EX_h2_e"] = 1000
medium["EX_octa_e"] = 0
medium["EX_ppa_e"] = 0
medium["EX_pta_e"] = 0
medium["EX_hpta_e"] = 0
medium["EX_co2_e"] = 0
medium["EX_for_e"] = 0

model.medium = medium

```

**Setting reaction constraints**

In the model, see **PART IV: SET ADDITIONAL CONSTRAINTS**

To turn off end-products, use the ```.knockout()``` function from CobraPy to turnoff exhcange reactions. If we wanted to prevent the model from producing formate, for instance, we could use the following:

```
#Turn off end products
model.reactions.EX_for.knock_out()
```

To turn off a reaction, use the ```.knockout()``` function from CobraPy. As an example, if we want to turn of Coenzyme A Transferase for all products from 4-8 carbons in length, we could use the following:

```
#Turn off CoaT
model.reactions.CoATC4.knock_out()
model.reactions.CoATC6.knock_out()
model.reactions.CoATC8.knock_out()
model.reactions.CoATC5.knock_out()
model.reactions.CoATC7.knock_out()
```

**Constraining products**



**Constraining "growth rates"**

**Considering transport energetics**


