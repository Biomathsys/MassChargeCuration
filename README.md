MassChargeCuration
========================================

Python module to automatically curate the mass and charge assignments for metabolites in a metabolic model.


[![License (MIT Licence)](https://img.shields.io/badge/license-MIT-blue.svg?style=plastic)](https://opensource.org/licenses/MIT)
![Lines of Code](https://img.shields.io/tokei/lines/github/Finnem/MassChargeCuration?color=orange&style=plastic)
![Download count](https://img.shields.io/github/downloads/Finnem/MassChargeCuration/total.svg?style=plastic)

----
*Authors*: [Finn Mier](https://uni-tuebingen.de/en/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/systems-biology/team/),
[Reihaneh Mostolizadeh](https://uni-tuebingen.de/en/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/systems-biology/team/dr-reihaneh-mostolizadeh/)


*Repository*: [https://github.com/Finnem/MassChargeCuration](https://github.com/Finnem/MassChargeCuration)



► Getting started with MassChargeCuration
----------------------------

Requires Microsoft Research's Z3 to run. You can download it from https://github.com/Z3Prover/z3.

In order to use this module, you need to install it first, e.g. by running `pip install -e .` in the folder you downloaded this repository to.

To apply it to a model, simply load your model via cobrapy and instantiate a MassChargeCuration class with your model as parameter.
If you are unsure about your models [identifiers.org](http://identifiers.org/) annotation, we also recommend to pass `update_ids = True` to the constructor. This will update any id used in the balancing effort, but will take longer (~15 minutes for 1200 Metabolites).
```
from MCC import MassChargeCuration
import cobra
model = cobra.io.read_sbml_model(model_path)
balancer = MassChargeCuration(model, update_ids = True)
```
If you have downloaded the BioCyc database, you can also pass the folder which contains it as an argument. In addition to that, the balancer will attempt to download other databases if they are not present, unless `no_local = True` is passed to it. You can pass a directory path for these databases as well via `data_path = database_path`.
```
balancer = MassChargeCuration(model, update_ids = True, data_path = "../database_path", biocyc_path = "../data/25.1/data")
```

Once the balancer has finished you can use it to generate some reports about its results.
A visual comparison of the original model and the new model can be generated with `balancer.generate_visual_report()`.
A report about unbalanced or noticable reactions can be generated with `balancer.generate_reaction_report(f"{model.id}_reactions")`.
A report about the chosen assignments for every metabolite can be generated with balancer.generate_metabolite_report(f"{model.id}_metabolites").

As an example we can look at all metabolites for which no or incomplete information was available with:
```
import pandas as pd
df = balancer.generate_metabolite_report(f"{model.id}_metaoblites")
pd.set_option('display.max_rows', None) # displays entire DF, takes a while
df[df["Inferrence Type"] != "Clean"]
```

Or at those metabolites where the assignment now differs from the original model with:
```
df = balancer.generate_metabolite_report(f"{model.id}_metaoblites")
pd.set_option('display.max_rows', None) # displays entire DF, takes a while
df[df["Similarity"] != "Same"]
```


☮ Licensing and distribution
----------------------------

MassChargeCuration is Copyright (C) 2021-2022 by the following organization:

The University of Tuebingen, Germany

MassChargeCuration is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation.
