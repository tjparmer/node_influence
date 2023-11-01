Dynamical Methods for Target Control of Biological Networks
=======================================================

Code to analyze node influence in Boolean and discrete-state networks, in particular to find the domain of influence (DOI) of a seed set for the purpose of target control  


**Tutorials:**

- Tutorial_Dynamical_Methods_for_Target_Control.ipynb shows how to create generalized threshold networks (GTNs) with various representations of the transfer functions, calculate DOI on these networks, and compare different methods (GTNs, IBMFA) for target control in biological networks (e.g., the _drosophila_ single-cell SPN).


**Scripts:**

- gtn_construction.py (code to create GTNs)
- mean_field_computations.py (code to run IBMFA)
- simulations.py (code to run and analyze simulations)
- brute_force_computations.py (code to run and analyze brute-force calculations for small networks)
- modules.py (utility functions related to influence pathways and finding DOI for generalized threshold networks)
- utils.py (utility functions for node influence)

**Original notebooks:**

Note: the scripts above were created from the functions originally used in the jupyter notebooks below.  If there's a bug in the above code, you may refer to the original function defined in the notebook to help troubleshoot.
- See network_influence.ipynb for analysis and comparison of dynamical methods for finding node influence in biological networks from the Cell Collective repository (http://cellcollective.org/) [1]
- See network_influence_results.ipynb for results of the analysis.  This notebook was intended to have limited dependencies on scripts or python packages but requires rather saved data files on which to make calculations. 
- See Comparing Influence in Boolean Networks.ipynb for figure results used in the paper.


The corresponding paper has been published in Royal Society Open Science [2] with corresponding Zenodo repository https://zenodo.org/records/10064204

References:
---------

[1] Tom´aˇs Helikar, Bryan Kowal, Sean McClenathan, Mitchell Bruckner, Thaine Rowley, Alex Madrahimov, Ben Wicks, Manish Shrestha, Kahani Limbu, and Jim A
Rogers. The cell collective: toward an open and collaborative approach to systems biology. BMC systems biology, 6(1):1–14, 2012.

[2] Parmer Thomas and Radicchi Filippo. 2023. Dynamical methods for target control of biological networks. R. Soc. open sci.10230542230542. https://doi.org/10.1098/rsos.230542
