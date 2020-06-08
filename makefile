# Compiler
F2PY = f2py

# Compile flags
FLAGS = -c

# Macros, Current Path 
PWD = $(CURDIR)
# Macros, not useful now
SPICELIB_Path = /Users/woodywu/Desktop/Research/Software/toolkit_libspicelib/lib \
SPICELIB_Name = spicelib \
SOFALIB_Path = /Users/woodywu/Desktop/Research/Proper_Element/thalassaMex/librarySetup/sofa/20180130/f77/src/ \
SOFALIB_Name = sofa \
SUBTHALASSALIB_Path = /Users/woodywu/Desktop/Research/Software/thalassa_f2py/thalassa/ \
SUBTHALASSALIB_Name = subthalassa

# Source and object files
OBJECTS = cart2coe4c.o kinds.o phys_const.o settings.o io.o kepler.o cart_coe.o nsgrav.o \
sun_moon.o drag_exponential.o US76_PATRIUS.o $(J77) nrlmsise00_sub.o srp.o \
perturbations.o initialize.o integrate.o $(SLSODAR) $(DLSODAR) $(FORMUL) \
propagate.o auxiliaries.o

# create dynamic library
#pytha.co: ./thalassa/libsubthalassa.a
#	$(F2PY)  -L/Users/woodywu/Desktop/Research/Software/thalassa_python/thalassa -lsubthalassa -L/Users/woodywu/Desktop/Research/Software/toolkit_libspicelib/lib -lspicelib -L/Users/woodywu/Desktop/Research/Proper_Element/thalassaMex/librarySetup/sofa/20180130/f77/src/ -lsofa -c ./thalassa/thalassaSub.f90 -m pytha

pytha.co: ./thalassa/libsubthalassa.a
	$(F2PY)  -L/$(PWD)/thalassa -lsubthalassa -L/$(PWD)/Lib/ -lspicelib -lsofa -c ./thalassa/thalassaSub.f90 -m pytha
	
./thalassa/libsubthalassa.a:
	make -C ./thalassa libsubthalassa.a

.PHONY: clean
clean:
	rm -rf *.so
	
clean_thalassa:
	cd ./thalassa && rm -rf *.o *.a
