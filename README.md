# GetLaplaceCoefficient

Program to calculate the Laplace coefficient b_s^(j)(alpha) at the value alpha, which is used in dynamical calculations involving secular interactions or mean-motion resonances.

The program returns b_s^(j)(alpha) as defined by Eq. 7 in Pearce & Wyatt (2015), which was taken from Murray & Durmott (1999). This is

b_s^(j)(alpha) = 1/pi integral_0^{2pi} cos(j*y) / (1 - 2*alpha*cos(y) + alpha^2)^s dy.

Note that the superscript (j) does *not* indicate b_s(alpha) to the power of j, but is rather the commonly used nomenclature to define the Laplace coefficient evaluated at a particular s, j and alpha.

To access the program, first download and unpack the ZIP file (press the green 'Code' button on GitHub, then 'Download ZIP', then unzip the file on your computer).

To use the program, simply edit the values in the 'User Inputs' section, then run the code. This will print out the value of b_s^(j) evaluated at alpha, and plot b_s^(j)(alpha) from alpha=0 to 1. 

This code was written by Tim D. Pearce. Feel free to use it, and if the results go into a publication, then please cite Costa, Pearce & Krivov (2023, submitted). Also, please let Tim know if you find any bugs or have any requests.
