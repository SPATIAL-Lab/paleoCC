# Fern Spike Box Model

Box model and driver scripts to simulate global carbon cycle response to landscape 
de-vegetation during the PETM.

## Model

The model is written in C++ and was compiled using Microsoft Visual Studio 2022 
Community Edition. The solution file fernSpike.sln should open the full solution
in VS and allow local compilation. fernSpike.cpp contains the main function and 
most of the model code.

The executable fernspike.exe was compiled on a 64-bit PC running Windows 11 and 
may be functional on similar systems.

The model reads configuration information and initial conditions from the following
text files:

- config.txt
- tmp00.txt
- tmp01.txt

Model output is saved in subdirectories corresponding to the model case name (e.g.,
base3k). This consists of a single comma-separated text file (exogenic.txt) 
containing a selection of state variables and fluxes saved at specified time steps 
throughout the simulation and a set of 25 text files (sediment00.txt, ...) storing
seafloor sediment properties. Each of the sediment files corresponds to a different
subsurface depth layer (from the sediment surface = 00 to 0.5 meters = 24), and
records the mass of carbonate and silicate, and the carbonate carbon isotope value
in each of 15 bathymetric intervals defined in the model (shallowest = 1).

## Driver

paleoCC.Rproj is a RStudio project file that will open a set of scripts used to
run and parse and plot output from the model. These scripts can be found in the
fernSpike subdirectory.

- experiments.R includes code that conducts and saves summary output from the 
experiments conducted for the Nelissen et al. manuscript. 
- helpers.R is source code for functions used in experiments.R. Of note, 
write.config() saves model parameter values provided as arguments in config.txt
where they can be read by the model.

The fernSpike subdirectory also includes .csv versions of the model output for 
each experiment, which were used in generating manuscript figures.