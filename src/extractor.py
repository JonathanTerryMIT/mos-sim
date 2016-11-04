'''
Parameter Extractor

dataFile:	String	Name of the file for which you wish to extract the MOSFET parameters. Should be located in the "../data" directory


'''

from definitions import *
import numpy as np
from matplotlib import pyplot as plt


class Extractor(object):

	def __init__(self, dataFile, thickness = 10.5, verbose = True):

		# Parse and load the datafile
		self.dataFile = dataFile
		self.thickness = thickness
		self.verbose = verbose


		self.info()
		self.parse()

		if verbose:
			print extractorASCII
			print "[+] Running extraction on: " + self.dataFile

			if (len(self.idvgTests) == 0): 
				print u"\t[FATAL] Unable to load ID-VGS data."
			if (len(self.idvdTests) == 0): 
				print u"\t[FATAL] Unable to load ID-VDS data."

			if (len(self.idvdTests) == 0 or len(self.idvgTests) == 0):
				return

			print "\n"
			print "\tCOMPILING EXPERIMENTAL DATA..."
			print u"\t[\u2713] Found " + str(len(self.idvgTests) + len(self.idvdTests)) + " test cases."
			print "\n"
			print "\tEXTRACTING PARAMETERS..."
			print u"\t[\u2713] Assuming channel width of " + str(self.width) + " nanometers."
			print u"\t[\u2713] Assuming channel length of " + str(self.length) + " nanometers."
			print u"\t[\u2713] Assuming oxide thickness of " + str(thickness) + " nanometers."

		# Extract the parameters
		self.runExtraction()

	def info(self):

		lengthWidth = self.dataFile.split("W")
		lengthWidth = lengthWidth[1].split("L")

		self.length = float(lengthWidth[1])
		self.width = float(lengthWidth[0])

		self.widthLength = self.width/self.length

		return

	def parse(self):

		try:
			# Load the experimental data (ID-VGS)
			data = np.loadtxt('../data/' + self.dataFile + '.idvg', skiprows=2)		

			## Splits the data into test cases based on VGS and VBS
			parsedTests = []

			tests = np.split(data, np.where(np.diff(data[:,0]))[0] + 1)

			for test in tests:

				parsedTests += np.split(test, np.where(np.diff(test[:,2]))[0] + 1)

			self.idvgTests = parsedTests

		except:
			self.idvgTests = []

		try:
			# Load the experimental data (ID-VDS)
			data = np.loadtxt('../data/' + self.dataFile + '.idvd', skiprows=2)		

			## Splits the data into test cases based on VDS and VBS
			parsedTests = []

			tests = np.split(data, np.where(np.diff(data[:,1]))[0] + 1)

			for test in tests:

				parsedTests += np.split(test, np.where(np.diff(test[:,2]))[0] + 1)

			self.idvdTests = parsedTests

		except:
			self.idvdTests = []

		return

	def extractCap(self):
		self.cap = epsilonSi / (self.thickness*10**(-7))
		if self.verbose: print u"\t[\u2713] Calculated oxide thickness of " + str(self.cap) + " Farad."
		return

	def extractMobility(self):

		# Find an appropriate part of the curve
		lowVoltageTest = self.idvgTests[0]

		for i in range(len(lowVoltageTest)):

			if lowVoltageTest[i][1] > 0.9:

				idsLower = lowVoltageTest[i][3]
				idsUpper = lowVoltageTest[i+1][3]

				vgsLower = lowVoltageTest[i][1]
				vgsUpper = lowVoltageTest[i+1][1]

				break


		plt.plot(lowVoltageTest[:, 1], lowVoltageTest[:, 3])
		plt.show()

		transconductance = (idsUpper - idsLower) / (vgsUpper - vgsLower)
		self.mobility = transconductance / (self.widthLength*self.cap*lowVoltageTest[0][0])

		if self.verbose: print u"\t[\u2713] Extracted mobility of " + str(self.mobility) + " cm^2 / volt-sec."

		return

	def extractThreshold(self):

		# Extract the relevant tests
		validTests = []
		minVds = 10

		for test in self.idvgTests:
			if test[0][0] <= minVds + 0.01:
				minVds = test[0][0]
				validTests.append(test)

		# Use the relevant tests to figure out the threshold
		for test in validTests:
			pass


		return


	def triodeEquation(self):


		return

	def runExtraction(self):

		self.extractCap()
		self.extractMobility()
		self.extractThreshold()

		return



extractor = Extractor('W25000L25000')