# Detecting-Demyelination-with-Anomalous-Diffusion-NMR
The project was part of the research done for my Master's Degree in Physics. It aims at demonstrating how Nuclear Magnetic 
Resonance - NMR and anomalous diffusion theory can be used for the detection of demyelination in white matter axonal tracts. 
Particularly, the anomalous diffusion theory parameters D and gamma were evaluated and were shown to be potential biomarkers 
for the early diagnosis of Multiple Sclerosis. For further details, see the paper attached to this project (Modeling biological
water diffusion in complex tissues.pdf).

# Built with
Matlab R2016b

# Executing the code
All code files used in this project are included in the folder modules. 

"Geometry" is the main module used for creating the simulation space. It lets the user choose between two kinds of space. If a 
model of an axonal tract is desired, the module calls submodule "axons_256". The second possibility is a space consisting of
spheres of equal diameter (for which the submodule "monospheres" is called). As is apparent, one can implement submodules 
to create any thinkable kind of geometry and the new submodules can be easily incorporated to the program's structure.

There are two modules whose names begin with "random_walk_function": one was devised for use with simulations of axonal tracts
and the other for simulations of spaces consisting of spheres of equal diameter. In both cases, the module is responsible for
the creation of random walk trajectories inside the respective spaces.

Finally, the modules whose names begin with "Signal_X_Experiments" simulate the acquisition of the NMR signal. Analogously, 
there are two modules: one for use with models of white matter axons; one for use with the simulations of spheres.

# Author
Tiago Monteiro Cardoso

# Acknowledgments
For the Mittag-Leffler function implementation I used the code provided by Igor Podlubny and Martin Kacenak in 
<https://www.mathworks.com/matlabcentral/fileexchange/8738-mittag-leffler-function>



