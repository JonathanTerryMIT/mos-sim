# Numerical methods libaries and plotting
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import newton

# Definitions of constants
from definitions import *

##-------------------Source Code (as opposed to Drain Code)------------------

'''
MOSFET Class

sourceBody: 	float 	representing the sourceBody bias
drainSource: 	float	representing drainSource voltage bias
gateSource: 	float 	representing gateSource votlage bias

'''

class MOSFET(object):

	def __init__(self, flatband, acceptorDoping, oxideThickness, width, length, mobilityParameters, verbose = False, dataPoints = 250):

		# Given parameters
		self.flatband = flatband
		self.width = width
		self.length = length
		self.mobilityParameters = mobilityParameters

		# Derived parameters
		self.oxideCapacitance = epsilonSi / (oxideThickness*(1e-7))
		self.fermiPotential = thermalVoltage*np.log(acceptorDoping/intrinsicCarriers)
		self.bodyCoeff = 0.53*(oxideThickness/10.0)*np.sqrt(acceptorDoping/(1e17))
		self.threshold = self.flatband + 2*self.fermiPotential + 6*thermalVoltage + self.bodyCoeff*np.sqrt(2*self.fermiPotential + 6*thermalVoltage)

		# Default simulation parameters
		self.verbose = verbose
		self.dataPoints = dataPoints

		# Storage parameters
		self.drainSourceCurrents = []

		if self.verbose:
			self.printSimulationParameters()


	def model(self, sourceBody, drainSource, gateSource):

		# Create the range needed
		self.drainSourceRange = np.linspace(drainSource[0], drainSource[1], self.dataPoints)
		psiRange = np.linspace(drainSource[0] + sourceBody, drainSource[1] + sourceBody, self.dataPoints)

		# Invoke vector processing
		diffVectorized = np.vectorize(self.diffusionCurrent)
		driftVectorized = np.vectorize(self.driftCurrent)

		for i in range(len(gateSource)):

			# Redefine voltages
			self.sourceBody = sourceBody
			self.gateBody = gateSource[i] + sourceBody

			# Print relevant information
			if self.verbose:
				print "[+] Running simulation sweep for gate-body voltage of: " + str(self.gateBody)

			# Calculate and memoize surface potentials and mobility
			self.potentials = np.vectorize(self.surfacePotential)(psiRange)
			self.mobility = self.calculateMobility()

			# Compute current components and graph
			diff = diffVectorized(self.potentials)
			drift = driftVectorized(self.potentials)
			drainSourceCurrent = np.add(diff, drift)

			self.drainSourceCurrents.append(drainSourceCurrent)
	
		return

	def diffusionCurrent(self, psiL):

		expansion = (psiL + self.potentials[0])/2.0
		alpha = 1.0 + (self.bodyCoeff/(2.0*np.sqrt(expansion)))
		
		mobilityComponent = (self.width/self.length)*self.mobility*self.oxideCapacitance*alpha*thermalVoltage
		drainComponent = psiL - self.potentials[0]

		return milliampScale*mobilityComponent*drainComponent

	def driftCurrent(self, psiL):

		expansion = (psiL + self.potentials[0])/2.0
		
		mobilityComponent = (self.width/self.length)*self.mobility*self.oxideCapacitance
		thresholdComponent = self.gateBody - self.flatband - expansion - self.bodyCoeff*np.sqrt(expansion)
		drainComponent = psiL - self.potentials[0]

		return milliampScale*mobilityComponent*thresholdComponent*drainComponent

	# Direct calculation of surface potential via newton-raphson
	def surfacePotential(self, contactBody):

		self.contactBody = contactBody
		
		# Attempt to speed up calculations
		try:
			guess = self.potentials[-1]
		
		except:
			guess = 1

		return newton(self.implicitPotential, guess, maxiter=500)

	# Implicit anonymous function
	def implicitPotential(self, psi):

		return self.gateBody - self.flatband - psi - self.bodyCoeff*np.sqrt(psi + thermalVoltage*np.exp((psi - 2*self.fermiPotential - self.contactBody)/thermalVoltage))


	def calculateMobility(self):

		# Mobility parameters of form (mu, [theta, thetaB])

		if len(self.mobilityParameters) == 1:
			return self.mobilityParameters[0]

		elif len(self.mobilityParameters) == 2:

			return self.mobilityParameters[0] / (1 + self.mobilityParameters[1]*(self.gateBody - self.sourceBody - self.threshold))

		else:
			return self.mobilityParameters[0] / (1 + self.mobilityParameters[1]*(self.gateBody - self.sourceBody - self.threshold) + self.mobilityParameters[2]*self.sourceBody)


	'''
	Function to plot MOSFET data.

	sourceBody: 	float 	representing the sourceBody bias
	drainSource: 	2-tuple representing drainSource voltage range
	gateSource: 	2-tuple representing gateSource votlage range

	'''
	def plotModel(self):

		for i in range(len(self.drainSourceCurrents)):
			plt.plot(self.drainSourceRange, self.drainSourceCurrents[i])

		plt.show()

		return

	def getSimulationData(self):
		return (self.drainSourceRange, self.drainSourceCurrents)

	def printSimulationParameters(self):
		
		# Print them
		print "[+] MOSFET initialized with the following parameters:"
		print "\t[-] Oxide Capacitance: " + str(self.oxideCapacitance)
		print "\t[-] Fermi Potential: " + str(self.fermiPotential)
		print "\t[-] Body Effect Coefficient: " + str(self.bodyCoeff)
		print "\t[-] Approximate threshold for strong inversion: " + str(self.threshold)