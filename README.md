# ISR_vol_img
Developers: Juan Araújo, Francisco López

## Supporting code for the manuscript:
Efficient Computation and visualization of Ionospheric Volumetric Images for the enhanced interpretation of Incoherent Scatter Radar Data,
Juan Araújo<sup>c</sup>, Francisco López, Stefan Johansson, Assar Westman, Madelen Bodin
<sup>c</sup> Corresponding author. 

## Software requirements (Linux systems)
+ g++ version 9.4.0, https://gcc.gnu.org/
+ HDF5 C++ API (1.12.2), https://docs.hdfgroup.org/archive/support/HDF5/doc1.8/cpplus_RM/index.html
+ OpenGL Mathematics (GLM), https://www.opengl.org/sdk/libs/GLM/

## ISR data
The data used for this application is in .hdf5 format and must be stored under the directory data/
Download AMISR data from http://cedar.openmadrigal.org/, following the instructions provided there.

Examples of data files used in the manuscript:
pfa150317.004.hdf5, pfa170128.002.hdf5, pfa200116.002.hdf5, ras161121.002.hdf5, ras190510.004.hdf5, ras200113.004.hdf5

## Compilation  (Linux systems):
+ g++ -std=c++17 -Wall -O3 -mavx -mfma -mavx2 -march=native -c interpolation_multi_parallel.cpp -o interpolation_multi_parallel.o -fopenmp;
+ g++ -std=c++17 -Wall -I/usr/lib/x86_64-linux-gnu/hdf5/serial/lib/include -I/src interpolation_multi_parallel.o main.cpp -L/usr/lib/x86_64-linux-gnu/hdf5/serial/lib -o main -lhdf5_cpp -lhdf5 -fopenmp;

## Execution:
./main -F [filename] -[options]

### Options:
-n: number of 1D interpolation points. Total interpolation points is the cube.
-v: verbose on interpolation process
-t: number of threads in parallel interpolation
-p: Shepard parameter

### Examples:
+ ./main -F ras161121.002
+ ./main -F ras161121.002 -n 64
+ ./main -F ras161121.002 -n 64 -p 3 
+ ./main -F ras161121.002 -n 64 -p 3 -t 8

## Output:
The process of the computation of the volumetric images is output on the terminal prompt.
The result of the computation of the volumetric images is written in .vtk files and stored under the directory output/ 
Output .vtk files can be visualized on the Paraview software (https://www.paraview.org/) using volume rendering.
