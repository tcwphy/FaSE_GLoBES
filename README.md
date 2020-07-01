# FaSE_GLoBES
This package serves to provide a toolkit for GLoBES to test general flavor symmetry models in neutrino oscillation experiments. Two files are provided in this package: 'FaSE_GLoBES.c'  and 'model-input.c'. To use the toolkit, the user only needs to include their binary files to compile the code. The function of 'FaSE_GLoBES.c' is to simulate the probability profile and implement the prior values for standard neutrino oscillation parameters, which are translated from the model parameters given by the user. The user-defined input, which includes the model set up and any restrictions on the model parameters, are defined in 'model-input.c'. More details can be found in the user manual in '/doc/FaSE_GLoBES.pdf'.

This toolkit is not a standalone program, but it requires an installation of GLoBES version 3.0.0 or higher. For those who are using GLoBES for the first time, please visit the GLoBES website: 'https://www.mpi-hd.mpg.de/personalhomes/globes/'.

If the user has any problems, suggestions or comments, please contact us:
Jian Tang 'tangjian5@mail.sysu.edu.cn' & TseChun Wang 'tcwphy@gmail.com'.
