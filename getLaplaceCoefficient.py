
'''Program to calculate the Laplace coefficient b_s^(j)(alpha) at the
value alpha, which is used in dynamical calculations involving secular
interactions or mean-motion resonances.


The program returns b_s^(j)(alpha) as defined by Eq. 7 in Pearce & Wyatt
(2015), which was taken from Murray & Durmott (1999). This is

b_s^(j)(alpha) = 1/pi integral_0^{2pi} cos(j*y) / (1 - 2*alpha*cos(y) + alpha^2)^s dy.

Note that the superscript (j) does *not* indicate b_s(alpha) to the power
of j, but is rather the commonly used nomenclature to define the Laplace 
coefficient evaluated at a particular s, j and alpha.


To use the program, simply edit the values in the 'User Inputs' section
below, then run the code. This will print out the value of b_s^(j)
evaluated at alpha, and plot b_s^(j)(alpha) from alpha=0 to 1. 


This code was written by Tim D. Pearce. Feel free to use it, and if the
results go into a publication, then please cite Costa, Pearce & Krivov
(2023, submitted). Also, please let Tim know if you find any bugs or have
any requests.'''

############################### Libraries ###############################
import sys
import numpy as np
import math
from scipy import integrate
import matplotlib.pyplot as plt

############################## User Inputs ##############################
'''Parameters to change. Don't change anything outside of this section'''

# The values s and j in the Laplace coefficient denoted b_s^(j)
s = 1.5
j = 1.

# The value of alpha at which the Laplace coefficient is to be evaluated,
# where 0 <= alpha < 1
alpha = 0.5

# The number of alpha values to plot on the graph
numberOfAlphasToPlot = 1000

############################### Functions ###############################
def GetLaplaceCoefficient(j, s, alpha):
	'''Returns the Laplace coefficient defined b_s^j(alpha)'''
	
	laplaceCoefficientIntegrand = lambda psi: np.cos(j*psi) / (1.0 - 2*alpha*np.cos(psi) + alpha**2)**s
	
	laplaceCoefficient = 1.0 / np.pi * integrate.quad(laplaceCoefficientIntegrand, 0, 2*np.pi)[0]
	
	return laplaceCoefficient

#------------------------------------------------------------------------
def MakePlot(laplaceCoefficientAtInputAlpha, sString, jString):
	'''Plot the Laplace coefficient vs. alpha, for alpha between 0 and 
	1.'''

	print('Making plot...')
	PrintEmptyLine()
		
	# Define the alpha values
	alphaStep = 1.0 / (numberOfAlphasToPlot + 1)
	alphas = np.arange(alphaStep, 1, alphaStep)

	# Initialise the Laplace coefficients
	laplaceCoefficients = []

	# Calculate the Laplace coefficient at each alpha
	for alphaForLine in alphas:
		laplaceCoefficient = GetLaplaceCoefficient(j, s, alphaForLine)
		laplaceCoefficients.append(laplaceCoefficient)

	# If all values of the Laplace coefficient are negative, modify the
	# plot 
	maxLaplaceCoeff = -1e99
	
	for laplaceCoefficient in laplaceCoefficients:
		maxLaplaceCoeff = max(laplaceCoefficient, maxLaplaceCoeff)

	if maxLaplaceCoeff < 0:
		absoluteValueString = 'Absolute value of '
		
		for i in range(len(laplaceCoefficients)):
			laplaceCoefficients[i] = abs(laplaceCoefficients[i])

	else:
		absoluteValueString = ''
		
	# Plot the Laplace coefficents
	plt.semilogy(alphas, laplaceCoefficients)

	# Mark the Laplace coefficient at the user-inputted alpha value
	plt.plot(alpha, laplaceCoefficientAtInputAlpha, 'kx')
	
	# Set various plot parameters
	plt.xlim(0,1)

	plt.xlabel(r'$\alpha$')
	plt.ylabel(r'%sLaplace coefficient $b^{(%s)}_{%s}(\alpha)$' % (absoluteValueString, jString, sString))

	# Display the plot
	plt.show()

#------------------------------------------------------------------------
def PrintEmptyLine():
	'''Print an empty line (done this way to enable compatibility across 
	Python versions)'''
	
	if sys.version_info[0] < 3: print
	else: print()

#------------------------------------------------------------------------
def CheckUserInputsOK():
	'''Check the user inputs are OK'''
	
	# Initialise the flag for whether inputs are OK, and the list of
	# reasons why the inputs are not OK
	areUserInputsOK = True
	reasonsUnputsAreBad = []
	
	# All system parameters must be floats or ints, and non-nan
	for value in [s, j, alpha, numberOfAlphasToPlot]:
		if isinstance(value, (float, int)) == False or math.isnan(value):
			areUserInputsOK = False
			reasonsUnputsAreBad.append('All system parameters must be floats or ints, and not nan')
			break
			
	# Alpha must be between 0 and 1
	if isinstance(alpha, (float, int)):
		if alpha < 0 or alpha >= 1:
			areUserInputsOK = False
			reasonsUnputsAreBad.append('alpha must be in the range 0 <= alpha < 1')
	
	# Number of alphas to plot must be at least two
	if isinstance(numberOfAlphasToPlot, (float, int)):
		if numberOfAlphasToPlot < 2:
			areUserInputsOK = False
			reasonsUnputsAreBad.append('numberOfAlphasToPlot must be at least 2')

	# Warn the user if the inputs are bad
	if areUserInputsOK == False:
		print('***ERROR*** Problem(s) with user inputs:')
		for reasonUnputsAreBad in reasonsUnputsAreBad:
			print('     -%s' % reasonUnputsAreBad)
		PrintEmptyLine()
		
	return areUserInputsOK	

#------------------------------------------------------------------------
def PrintUserInputs():
	'''Print the user inputs'''
	
	print('Evaluating the Laplace coefficient b_s^(j) evaluated at alpha, where:')
	print('     s = %s' % s)
	print('     j = %s' % j)
	print('     alpha = %s' % alpha)
	PrintEmptyLine()
	
#------------------------------------------------------------------------
def PrintOutput(s, j, alpha, laplaceCoefficient):
	'''Print the output'''
	
	print('Result:')
	
	# Write s and j as ints, unless must be floats
	if int(s) == s: sString = '%s' % int(s)
	else: sString = '%s' % float(s)
		
	if int(j) == j: jString = '%s' % int(j)
	else: jString = '%s' % float(j)
	
	print('     b_%s^(%s)(%s) = %s' % (sString, jString, alpha, laplaceCoefficient))
	PrintEmptyLine()
		
	return sString, jString
		
################################ Program ################################
PrintEmptyLine()
	
# Check user inputs fine
areUserInputsOK = CheckUserInputsOK()

# Continue if the user inputs are OK
if areUserInputsOK:

	# Print the user inputs
	PrintUserInputs()
	
	# Evalute Laplace coefficient and print
	laplaceCoefficientAtInputAlpha = GetLaplaceCoefficient(j, s, alpha)
	sString, jString = PrintOutput(s, j, alpha, laplaceCoefficientAtInputAlpha)
	
	# Make the plot
	MakePlot(laplaceCoefficientAtInputAlpha, sString, jString)

	print('Complete')
	PrintEmptyLine()

#########################################################################

