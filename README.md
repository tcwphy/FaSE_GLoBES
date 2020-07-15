# FaSE_GLoBES
This package serves to provide a toolkit for GLoBES to test general flavor symmetry models in neutrino oscillation experiments.


*INSIDE FaSE*
-----------------
Two files are provided in FaSE package: 'FaSE_GLoBES.c'  and 'model-input.c'. To use the toolkit, the user only needs to include their binary files to compile the code.

FaSE_GLoBES: The function of 'FaSE_GLoBES.c' is to simulate the probability profile and implement the prior values for standard neutrino oscillation parameters, which are translated from the model parameters given by the user.

model-input: The user-defined input, which includes the model set up and any restrictions on the model parameters, are defined in 'model-input.c'.



*DOWLOAD AND INSTALLATION*
---------------------------------------
This toolkit is not a standalone program, but it requires an installation of GLoBES version 3.0.0 or higher. For those who are using GLoBES for the first time, please visit the GLoBES website: 'https://www.mpi-hd.mpg.de/personalhomes/globes/'.

The user can download FaSE simply from https://github.com/tcwphy/FaSE_GLoBES. git directly. Or, the user can also download it via:
git clone https://github.com/tcwphy/FaSE_GLoBES.git

As FaSE works based on GLoBES, the user should use the GLoBES Makefile, but include two binary files of FaSE ('FaSE_GLoBES.o'  and 'model-input.o'). To do so, we include the script in the Makefile:
{makefile_execution}:  my_program.o FaSE_GLoBES.o model-input.o
gcc my_program.o FaSE_GLoBES.o model-input.o -o my_executable $(LDFLAGS) $(local_LDFLAGS)

where {compile_execusion} is the makefile command for the program my_program.o and the execution my_executable is the output. After giving this script in the GLoBES Makefile, the user needs to makefile the program my_program on the terminal:
makefile {makefile_execution}

And, the user can execute the program with:
./my_executable



*CONTACT US*
------------------
More details can be found in the user manual in '/doc/FaSE_GLoBES.pdf'.
If the user has any problems, suggestions or comments, please contact us:
Jian Tang 'tangjian5@mail.sysu.edu.cn' & TseChun Wang 'tcwphy@gmail.com'.
