all: generation digitizer clusterizer validation_digitizer validation_clusterizer

validation: validation_digitizer validation_clusterizer

generation:
	cmsRun Generation/python/Generation.py

digitizer:
	cmsRun Digitizer/python/Digitizer.py

clusterizer:
	cmsRun Clusterizer/python/Clusterizer.py

validation_digitizer:
	cmsRun Digitizer/test/RunDigitizerValidation.py

validation_clusterizer:
	cmsRun Clusterizer/test/RunClusterizerValidation.py
