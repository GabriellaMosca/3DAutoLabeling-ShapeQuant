# 3DAutoLabeling-ShapeQuant
Plugin for MorphoMechanX to classify in semi-automated way cell type in the ovule (L1 and pSMC/SMC  need to be selected) and to quantify principal axes of elongations in each cell

To run this model on the test mesh provided, get MorphoMechanX first (you can find it here: www.morphomechanx.org). 
After installing MorphoMechnX(part of MorphoDynamX modeling environment) this plugin won't need any further installation.

Download the plugin from this git repository and put it into a folder.

From a terminal move into the folder where the plugin is places and type:
- make clean
- make run

A window will open showing already the 3D mesh with cell type classification active. 
To test how cell classification works, select the wedge simbol from the left panel of the window (prolonged clicked of the left  mouse button will allow to select between triangle and wedge) and with that select manyally all the L1 cells in the area of interest. 
The run: 
Model/CCF/50 Set Cell Type 
with the option "L1" selected.
For the Ovule it is necessary to specify the pSMC/SMC manually as well (to see inside use the Clippin Plane option), run again: 
Model/CCF/50 Set Cell Type 
with the option "pSMC" selected. 

After this select the region of ovule you are interested to for your analysis with the cube symbol on the left panel of the window (prolonged click of the mouse will allow to select between rectangle and cube - all shown with a diagonal)

and run: 
Model/CCF/60 Shape and Connectivity Quantifiers/00 Global Shape and Connectivity Quantifier Process

This will process the cell and assigne them a type based on connectivity (L2/L3, Companion Cell) and measure the principal elongation axes of the selected cells.
The results will be reported in a csv file. Its name can be specified into the "File Name" field in the process :

Model/CCF/60 Shape and Connectivity Quantifiers/05 Write Shape and Connectivity Quantifiers  (run this process to save the fields to a file).

Within this Git repo, there is also a video Tutorial available. 

For further information refer to the Appendix 1 provided with the paper (Hernandez-Lagana, Mosca, Mendocilla-Sato et al.). The Appendix can be found as pre-print also here https://www.biorxiv.org/content/10.1101/2020.07.30.226670v3.supplementary-material under "Supplemental File 1"

For further questions, contact: gabriella.mosca@gmail.com



 
