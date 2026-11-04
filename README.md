The repository has been created to include a python script whose purpose is to provide an example of working fluid selection in a binary ORC power cycle operating with a typical medium-enthalpy geothermal source. The script is developed monolithically in a single file given the short length.

The algorithm computes a condensing Rankine Cycle with no superheating for each of the fluids that are present in CoolProp 8.0.0. (with adequate critical temperature for a subcritical cycle to be possibile) at optimized maximum cycle pressure, limited by the minimum temperature difference inside the evaporator. The data and hypothesis can be deduced from the script.

The script has been developed for my bachelor's thesis in Energy Engineering with title "ORC systems for the exploitation of medium‑to‑low enthalpy geothermal sources" and uploaded to github on November 2026.

## Requirements and installation

```
conda env create -f environment.yml
conda activate orc-wfs
python sub_1plv.py
```

The script has been written using Spyder as IDE and Anaconda as the Python distribution and package manager. Tne above lines are needed to reproduce the minimal correct setup to run the code and should be run in the Anaconda Prompt terminal.

