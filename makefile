# Makefile

sourceVoronoi  = infnan.o constants_mod.o vector_mod.o common_mod.o \
	VoronoiVoropp.o VoronoiFortranInterface.o \
	interpolation_mod.o \
	set_input_mod.o hydro_mod.o ph_mod.o composition_mod.o \
	continuum_mod.o ionization_mod.o pathIntegration_mod.o \
	grid_mod.o voronoi_grid_mod.o dust_mod.o emission_mod.o photon_mod.o \
	voronoi_photon_mod.o \
	update_mod.o output_mod.o iteration_mod.o \
        mocassin.o

sourcePlot  = infnan.o constants_mod.o vector_mod.o common_mod.o \
	VoronoiVoropp.o VoronoiFortranInterface.o \
	interpolation_mod.o \
	set_input_mod.o hydro_mod.o ph_mod.o composition_mod.o \
	continuum_mod.o ionization_mod.o pathIntegration_mod.o \
	grid_mod.o voronoi_grid_mod.o dust_mod.o emission_mod.o photon_mod.o \
	voronoi_photon_mod.o \
	update_mod.o output_mod.o iteration_mod.o \
        mocassinPlot.o

sourceWarm  = infnan.o constants_mod.o vector_mod.o common_mod.o \
	VoronoiVoropp.o VoronoiFortranInterface.o \
	interpolation_mod.o \
	set_input_mod.o hydro_mod.o ph_mod.o composition_mod.o \
	continuum_mod.o ionization_mod.o pathIntegration_mod.o \
	grid_mod.o voronoi_grid_mod.o dust_mod.o emission_mod.o photon_mod.o \
	voronoi_photon_mod.o \
	update_mod.o output_mod.o iteration_mod.o \
        mocassinWarm.o

source1  = source/infnan.f90 source/constants_mod.f90 source/vector_mod.f90 source/common_mod.f90 \
	source/VoronoiVoropp.f90 source/VoronoiFortranInterface.f90 source/interpolation_mod.f90 \
	source/set_input_mod.f90 source/hydro_mod.f90 source/ph_mod.f90 source/composition_mod.f90 \
	source/continuum_mod.f90 source/ionization_mod.f90 source/pathIntegration_mod.f90 \
	source/grid_mod.f90 source/voronoi_grid_mod.f90 source/dust_mod.f90 source/emission_mod.f90 source/photon_mod.f90  \
	source/voronoi_photon_mod.f90 source/update_mod.f90 \
	source/output_mod.f90 source/iteration_mod.f90 source/mocassin.f90 

source2  = source/infnan.f90 source/constants_mod.f90 source/vector_mod.f90 source/common_mod.f90 source/interpolation_mod.f90 \
	source/set_input_mod.f90 source/hydro_mod.f90 source/ph_mod.f90 source/composition_mod.f90 \
	source/continuum_mod.f90 source/ionization_mod.f90 source/pathIntegration_mod.f90 \
	source/grid_mod.f90 source/dust_mod.f90 source/emission_mod.f90 source/photon_mod.f90  \
	source/update_mod.f90 \
	source/output_mod.f90 source/iteration_mod.f90 source/mocassinWarm.f90 

sourceOutput  = infnan.o constants_mod.o vector_mod.o common_mod.o \
	VoronoiVoropp.o VoronoiFortranInterface.o \
	interpolation_mod.o \
	set_input_mod.o hydro_mod.o ph_mod.o composition_mod.o \
	continuum_mod.o ionization_mod.o pathIntegration_mod.o \
	grid_mod.o voronoi_grid_mod.o dust_mod.o emission_mod.o photon_mod.o \
	voronoi_photon_mod.o \
	update_mod.o output_mod.o iteration_mod.o \
        mocassinOutput.o

source4  = source/infnan.f90 source/constants_mod.f90 source/vector_mod.f90 source/common_mod.f90 source/interpolation_mod.f90 \
	source/set_input_mod.f90 source/hydro_mod.f90 source/ph_mod.f90 source/composition_mod.f90 \
	source/continuum_mod.f90 source/ionization_mod.f90 source/pathIntegration_mod.f90 \
	source/grid_mod.f90 source/dust_mod.f90 source/emission_mod.f90 source/photon_mod.f90  \
	source/update_mod.f90 \
	source/output_mod.f90 source/iteration_mod.f90 source/mocassinPlot.f90 
 
F90  = mpif90
CPP = g++
LIBS =	-lm
OPT1 = -O2 -fno-range-check
OPT2 = -fno-range-check -g -traceback
CFLAGS = -O2 -DNDIM=3
VORO++ = /Users/barbara/voronoi_mocassin_V0/voro++-0.4.6/src
VOROHEADERS = $(VORO++)
VOROLIB = $(VORO++)
VPATH = $(PWD)/source

.SUFFIXES: .cpp .f90 .c .o

%.o: %.f90
	$(F90) $(OPT1) -c $<

%.o: %.cpp
	$(CPP) $(CFLAGS) -c $< -I$(VORO++) 

mocassinVoronoi: $(sourceVoronoi)
	$(F90) -lc -lstdc++ -o mocassinVoronoi $(sourceVoronoi) -L$(VORO++) -lvoro++ -lm

mocassin:
	$(F90) $(OPT1) -o mocassin $(source1) $(LIBS)

mocassinPlot: $(sourcePlot)
	$(F90) -lc -lstdc++ -o mocassinPlot $(sourcePlot) -L$(VORO++) -lvoro++ -lm

mocassinWarm: $(sourceWarm)
	$(F90) -lc -lstdc++ -o mocassinWarm $(sourceWarm) -L$(VORO++) -lvoro++ -lm

mocassinOutput: $(sourceOutput)
	$(F90) -lc -lstdc++ -o mocassinOutput $(sourceOutput) -L$(VORO++) -lvoro++ -lm


clean:
	/bin/rm *.o *~ *.mod mocassin mocassinWarm mocassinOutput mocassinPlot mocassinVoronoi
