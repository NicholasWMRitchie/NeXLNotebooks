# -*- coding: utf-8 -*-
# Description: A template script for performing multiple different RCA analyses
# Author: Nicholas W. M. Ritchie
# Modified: 2-Mar-2017

analyst = "Ritchie, Ortiz-Montalvo, Vicenzi" # jl.System.getProperty("user.name") # Replace with your analyst = "your name"! 

beamEnergy = 20.0 # keV (Checked, not set)
ruleFile = "default.zrr" # or you can specify one in 'basicRCA', 'highZRCA' etc.
realTime = 0.4 # seconds per particle

# Define the elements and the associated standards in the vector file
stds = { "C": "C std.msa", "Al": "Al std.msa", "Cu":"Cu std.msa", "Ni":"Ni std.msa", 
	"Na": "NaF std.msa", "O": "MgO std.msa", "Si": "Si std.msa", # "Cl": "NaCl std.msa",
	"Fe": "Fe std.msa", "Ca": "CaF2 std.msa", "Cr": "Cr std.msa", "Ni":"Ni std.msa", 
	"Cu": "Cu std.msa", "Ti": "Ti std.msa", "Mg": "Mg std.msa", "Mn":"Mn std.msa", 
	"S" : "FeS2 std.msa", "Zn" : "Zn std.msa", "Ba":"BaF2 std.msa", "Mo": "Mo std.msa",
	"Ge": "Ge std.msa", "Bi" : "Bi std.msa"
}

#Set these to point to the directory contain your standard spectra and rules
stdsPath = "%s\\standards\\Combined\\25 keV" % rootPath
rulePath = "%s\\standards\\Rules" % (rootPath, )

# Default APA settings
def basicRCA(project, sample, tiling, stageZ):
	fov = 0.256 # mm
	maxFields = 0 # or 0 to analyse all fields
	maxPartPerField = 1000 #
	maxParticles = 10000
	seed = None # or an integer value (randomizes field selection for maxFields)
	measureThreshold = (32, 255) # lower, upper
	measureDwell = 4 # microseconds
	searchThreshold = (128, 255) # lower, upper
	searchDwell = 1 # microseconds

	vecs = buildVectors(stds,strip=["C","O"],path=stdsPath)
	rules = graf.BasicRuleSet("%s\\%s" % (rulePath, ruleFile))

	# Create the RCA object, initialize and execute the analysis
	rca=buildRCA(project, sample, vecs, rules, realTime=realTime, analyst=analyst)
	rca.setStageZ(stageZ)
	rca.setFieldOfView(fov, overlap = 1.0)
	rca.setMeasureThreshold(low=measureThreshold[0], high = measureThreshold[1], dwell = int(measureDwell*1000), measureStep = 8)
	rca.setSearchThreshold(low=searchThreshold[0], high = searchThreshold[1], dwell = int(searchDwell*1000), maxPartPerField=maxPartPerField, maxPart=maxParticles)
	try:
		rca.perform(randomizedTiling(tiling, seed=seed))
		#rca.perform(tiling)
	except jl.Throwable, ex:
		print "Error analyzing %s" % sample
		print str(ex)
	finally:
		rca.postSummary(jl.System.out)
	return rca.getZeppelin()

# APA settings biased towards high Z particles	
def highZRCA(project, sample, tiling, stageZ):
	fov = 0.256 # mm
	maxFields = 0 # or 0 to analyse all fields
	maxPartPerField = 1000 #
	maxParticles = 100000
	seed = None # or an integer value (randomizes field selection for maxFields)
	measureThreshold = (128, 255) # lower, upper
	measureDwell = 4 # microseconds
	searchThreshold = (160, 255) # lower, upper
	searchDwell = 1 # microseconds

	vecs = buildVectors(stds,strip=["C","O"],path=stdsPath)
	rules = graf.BasicRuleSet("%s\\%s" % (rulePath, ruleFile))

	# Create the RCA object, initialize and execute the analysis
	rca=buildRCA(project, sample, vecs, rules, realTime=realTime, analyst=analyst)
	rca.setStageZ(stageZ)
	rca.configEDS(vecs, rules, realTime=realTime, mode=POINT_MODE)
	rca.setFieldOfView(fov, overlap = 1.0)
	rca.setMeasureThreshold(low=measureThreshold[0], high = measureThreshold[1], dwell = int(measureDwell*1000), measureStep = 8)
	rca.setSearchThreshold(low=searchThreshold[0], high = searchThreshold[1], dwell = int(searchDwell*1000), maxPartPerField=maxPartPerField, maxPart=maxParticles)
	try:
		rca.perform(randomizedTiling(tiling, seed=seed))
		#rca.perform(tiling)
	except jl.Throwable, ex:
		print "Error analyzing %s" % sample
		print str(ex)
	finally:
		rca.postSummary(jl.System.out)
	return rca.getZeppelin()
	
# APA settings biased towards high Z particles	
def basicRCALowMag(project, sample, tiling, stageZ):
	fov = 0.512 # mm
	maxFields = 0 # or 0 to analyse all fields
	maxPartPerField = 1000 #
	maxParticles = 100000
	seed = None # or an integer value (randomizes field selection for maxFields)
	measureThreshold = (32, 255) # lower, upper
	measureDwell = 4 # microseconds
	searchThreshold = (128, 255) # lower, upper
	searchDwell = 1 # microseconds

	vecs = buildVectors(stds,strip=["C","O"],path=stdsPath)
	rules = graf.BasicRuleSet("%s\\%s" % (rulePath, ruleFile))

	# Create the RCA object, initialize and execute the analysis
	rca=buildRCA(project, sample, vecs, rules, realTime=realTime, analyst=analyst)
	rca.setStageZ(stageZ)
	rca.setFieldOfView(fov, overlap = 1.0)
	rca.setMeasureThreshold(low=measureThreshold[0], high = measureThreshold[1], dwell = int(measureDwell*1000), measureStep = 8)
	rca.setSearchThreshold(low=searchThreshold[0], high = searchThreshold[1], dwell = int(searchDwell*1000), maxPartPerField=maxPartPerField, maxPart=maxParticles)
	try:
		# rca.perform(randomizedTiling(tiling, seed=seed))
		rca.perform(tiling)
	except jl.Throwable, ex:
		print "Error analyzing %s" % sample
		print str(ex)
	finally:
		rca.postSummary(jl.System.out)
	return rca.getZeppelin()

# A search setting for the smallest particles (A = 0.05 to 10 sq micrometers)
def basicHighMagRCA(project, sample, tiling, stageZ):
	fov = 0.100 # mm
	maxFields = 0 # or 0 to analyse all fields
	maxPartPerField = 1000 #
	maxParticles = 10000
	seed = None # or an integer value (randomizes field selection for maxFields)
	measureThreshold = (32, 255) # lower, upper
	measureDwell = 4 # microseconds
	searchThreshold = (128, 255) # lower, upper
	searchDwell = 1 # microseconds

	vecs = buildVectors(stds,strip=["C","O"],path=stdsPath)
	rules = graf.BasicRuleSet("%s\\%s" % (rulePath, ruleFile))

	# Create the RCA object, initialize and execute the analysis
	rca=buildRCA(project, sample, vecs, rules, realTime=realTime, analyst=analyst)
	rca.setStageZ(stageZ)
	rca.setFieldOfView(fov, overlap = 1.0)
	rca.setMeasureThreshold(low=measureThreshold[0], high = measureThreshold[1], dwell = int(measureDwell*1000), measureStep = 8)
	rca.setSearchThreshold(low=searchThreshold[0], high = searchThreshold[1], dwell = int(searchDwell*1000), maxPartPerField=maxPartPerField, maxPart=maxParticles)
	rca.setMorphologyCriterion(semtr.RcaTranslator.AreaCriterion(0.05, 1.0))
	try:
		rca.perform(randomizedTiling(tiling, seed=seed))
		#rca.perform(tiling)
	except jl.Throwable, ex:
		print "Error analyzing %s" % sample
		print str(ex)
	finally:
		rca.postSummary(jl.System.out)
	return rca.getZeppelin()

	
	
# Begin: CONFIGURATION SECTION
# This section contains most of the items that should be changed from analysis set to analysis set.
# 1.  Create pts1, pts2, ..., ptsN outlining the N samples
# 2.  Modify 'analyses' to describe the various different analyses 
#     a. You can perform multiple analyses per sample
#     b. Choices of tilings inlude 'circularTiling', 'rectangularTiling' and 'boundaryTiling'
#     c. The broad definition of the analysis is defined by 'highZRCA', 'basicRCA', 'basicRCALowMag' (You can create your own or modified these.)

s2 = parseCoords("[{X:1.732,Y:14.410,Z:25.010,Rotate:-0.00,Tilt:-0.00}, {X:-3.991,Y:16.012,Z:25.010,Rotate:-0.00,Tilt:-0.00}, {X:-5.239,Y:9.838,Z:25.010,Rotate:-0.00,Tilt:-0.00}, {X:0.299,Y:8.510,Z:25.010,Rotate:-0.00,Tilt:-0.00}]")

s3 = parseCoords("[{X:-8.198,Y:8.420,Z:24.993,Rotate:-0.00,Tilt:-0.00}, {X:-14.445,Y:8.734,Z:24.993,Rotate:-0.00,Tilt:-0.00}, {X:-14.687,Y:1.932,Z:24.993,Rotate:-0.00,Tilt:-0.00}, {X:-8.697,Y:2.214,Z:24.993,Rotate:-0.00,Tilt:-0.00}]")

s4 = parseCoords("[{X:-8.667,Y:-3.207,Z:24.993,Rotate:-0.00,Tilt:-0.00}, {X:-14.367,Y:-3.827,Z:24.993,Rotate:-0.00,Tilt:-0.00}, {X:-12.781,Y:-11.652,Z:24.993,Rotate:-0.00,Tilt:-0.00}, {X:-7.500,Y:-10.509,Z:24.993,Rotate:-0.00,Tilt:-0.00}]")

s5 = parseCoords("[{X:5.100,Y:-15.840,Z:24.972,Rotate:-0.00,Tilt:-0.00}, {X:-0.616,Y:-18.030,Z:24.972,Rotate:-0.00,Tilt:-0.00}, {X:-3.015,Y:-11.614,Z:24.972,Rotate:-0.00,Tilt:-0.00}, {X:3.466,Y:-9.553,Z:24.972,Rotate:-0.00,Tilt:-0.00}]")

s6 = parseCoords("[{X:14.667,Y:-5.196,Z:24.993,Rotate:-0.00,Tilt:-0.00}, {X:10.537,Y:-3.385,Z:24.993,Rotate:-0.00,Tilt:-0.00}, {X:8.058,Y:-8.722,Z:24.993,Rotate:-0.00,Tilt:-0.00}, {X:11.817,Y:-10.501,Z:24.993,Rotate:-0.00,Tilt:-0.00}]")

s7 = parseCoords("[{X:-2.785,Y:-0.878,Z:24.993,Rotate:-0.00,Tilt:-0.00}, {X:-1.256,Y:-3.647,Z:24.993,Rotate:-0.00,Tilt:-0.00}, {X:3.156,Y:0.233,Z:24.993,Rotate:-0.00,Tilt:-0.00}, {X:1.417,Y:3.011,Z:24.993,Rotate:-0.00,Tilt:-0.00}]")


#tilt = parseCoords("[{X:-1.530,Y:-6.955,Z:19.064,Rotate:-0.00,Tilt:-0.00}, {X:2.383,Y:-6.235,Z:20.231,Rotate:-0.00,Tilt:-0.00}, {X:3.897,Y:-1.011,Z:21.849,Rotate:-0.00,Tilt:-0.00}, {X:0.152,Y:2.567,Z:22.299,Rotate:-0.00,Tilt:-0.00}, {X:-6.250,Y:-4.028,Z:19.203,Rotate:-0.00,Tilt:-0.00}, {X:-4.048,Y:-6.568,Z:18.579,Rotate:-0.00,Tilt:-0.00}, {X:-0.796,Y:-2.100,Z:20.781,Rotate:-0.00,Tilt:-0.00}]")

# WARNING: Working distance has been replaced by stageZ.  Now a stage motion instead of a change in focal distance.
analyses = (
   # ( "project", "sample", tiling, stageZ, rcaFunc ),
   ( "Castle Dust", "CD2 7-31", boundaryTiling(s2), 25.051, basicRCA ),
   ( "Castle Dust", "CD2 7-30", boundaryTiling(s3), 25.008, basicRCA ),
   ( "Castle Dust", "CD2 8-17", boundaryTiling(s4), 25.022, basicRCA ),
   ( "Castle Dust", "CD2 8-06", boundaryTiling(s5), 25.014, basicRCA ),
   ( "Castle Dust", "CD2 8-02", boundaryTiling(s6), 25.114, basicRCA ),
   ( "Castle Dust", "CD2 8-01", boundaryTiling(s7), 24.724, basicRCA ),
)

# End: CONFIGURATION SECTION

# Begin: RUN SECTION 
# Avoid modification
run = True
if jl.Math.abs(_ts.hvGetVoltage()-beamEnergy*1.0e3)>100.0:
	print "ERROR: Instrument not configured for %0.1f keV" % beamEnergy
	run = False
if _ts.hvGetBeam()<>1:
	print "ERROR: Electron beam not on!"
	run = False

# Check Z motions to reduce risk of collision

if run:
	results=[]
	try:
		for project, sample, tiling, stageZ, rcaFunc in analyses:
			try:
				result = rcaFunc(project, sample, tiling, stageZ)
				results.append(result)
			except jl.Throwable, ex:
				print "Error analyzing %s - %s" % (project, sample)
				print str(ex)
			if terminated:
				all = (jop.showConfirmDialog(MainFrame, "Terminate all analyses?", "multiRCA", jop.YES_NO_OPTION) == jop.YES_OPTION)
				if all:
					break
				else:
					terminated=False
	finally:
		if not terminated:
			turnOff()
else:
	print "Correct the problems and re-run."
# END: RUN SECTION
