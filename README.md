# 3DsemiAutoLabelingOvule
Plugin for MorphoMechanX to classify in semi-automated way cell type in the ovule (L1 and pSMC/SMC  need to be selected) and to quantify principal axes of elongations in each cell

To run this model on the test mesh provided, get MorphoMechanX first (release upon request, https://www.mpipz.mpg.de/MorphoGraphX/MorphoMechanX) and install it, this plugin works with git revision of MorphoMechanX e618301477af8a0bda2c344a7312195ba5acf379 and this info should be provided together with the software request). 
After installing MorphoMechnX( in turn plugin of MorphoDynamX) this plugin won't need any further installation.
Download the plugin from this git repository and put it into a folder.
From a terminal move into the folder where the plugin is places and type:
make run

A window will be opened showing already the 3D mesh with cell type classification active. 
To test how cell classification works, select the wedge simbol from the left panel of the window (prolonged clicked of the left  mouse button will allow to select between triangle and wedge) and with that select manyally all the L1 cells in the area of interest. 
The run: 
Model/CCF/50 Set Cell Type 
with the option "L1" selected.
For the Ovule it is necessary to specify the pSMC/SMC manyally as well (to see inside use the Clippin Plane option), run again: 
Model/CCF/50 Set Cell Type 
with the option "pSMC" selected. 

After this select the region of ovule you are interested to for your analysis with the cube symbol on the left panel of the window (prolonged click of the mouse will allow to select between rectangle and cube - all shown with a diagonal)
and run: 
Model/CCF/60 Shape Quantifier/00 Global Shape Quantifier Process

This will process the cell and assigned them a type based on connectivity (L2/L3, Companion Cell, measure the principal elongation axes of the selected.
The results can be reported in a csv file by running:
Model/CCF/70 Write Shape Quantifiers (specify the file name in the filed on the right).

A video tutorial will be provided soon

 
