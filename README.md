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

**Feeding the model**

In the model, see **Part III: Substrate Uptake**

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

**Turning reactions on and off**

```
#Turn off CoaT
model.reactions.LEO_CoATC4.knock_out()
model.reactions.LEO_CoATC6.knock_out()
model.reactions.LEO_CoATC8.knock_out()
model.reactions.LEO_CoATC5.knock_out()
model.reactions.LEO_CoATC7.knock_out()
```

**Constraining products**

**Constraining "growth rates"**

**Considering transport energetics**


