# Compiler
F2PY = f2py

# Compile flags
FLAGS = -c
# Link flags
SPICELIB_Path = /Users/woodywu/Desktop/Research/Software/toolkit_libspicelib/lib \
SPICELIB_Name = spicelib \
SOFALIB_Path = /Users/woodywu/Desktop/Research/Proper_Element/thalassaMex/librarySetup/sofa/20180130/f77/src/ \
SOFALIB_Name = sofa \
SUBTHALASSALIB_Path = /Users/woodywu/Desktop/Research/Software/thalassa_f2py/thalassa/ \
SUBTHALASSALIB_Name = subthalassa

# create dynamic library
pytha.co: ./thalassa/libsubthalassa.a
	$(F2PY)  -L/Users/woodywu/Desktop/Research/Software/thalassa_f2py/thalassa -lsubthalassa -L/Users/woodywu/Desktop/Research/Software/toolkit_libspicelib/lib -lspicelib -L/Users/woodywu/Desktop/Research/Proper_Element/thalassaMex/librarySetup/sofa/20180130/f77/src/ -lsofa -c ./thalassa/thalassaSub.f90 -m pytha

.PHONY: clean
clean:
	rm -rf *.so
