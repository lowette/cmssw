all:
	cmsRun Generation/python/Generation.py _n $(n)
	cmsRun Digitizer/python/Digitizer.py 
	cmsRun Clusterizer/python/Clusterizer.py
	cmsRun Digitizer/test/RunDigiValidation.py 
	cmsRun Clusterizer/test/RunClusterValidation.py

validation: validation_digitizer validation_clusterizer

generation:
	cmsRun Generation/python/Generation.py _n $(n) _output $(output)

digitizer:
	cmsRun Digitizer/python/Digitizer.py _input $(input) _output $(output)

clusterizer:
	cmsRun Clusterizer/python/Clusterizer.py _input $(input) _output $(output)

validation_digitizer:
	cmsRun Digitizer/test/RunDigiValidation.py _input $(input) _output $(output)

validation_clusterizer:
	cmsRun Clusterizer/test/RunClusterValidation.py _input $(input) _output $(output)

help:
	@echo "+--------------------------------------------------------+"
	@echo "|                        RunSteps                        |"
	@echo "+--------------------------------------------------------+"
	@echo "| List of available commands:                            |"
	@echo "+--------------------------------------------------------+"
	@echo "| make generation n=<int> output=<path>                  |"
	@echo "| make digitizer input=<path> output=<path>              |"
	@echo "| make clusterizer input=<path> output=<path>            |"
	@echo "| make validation_digitizer input=<path> output=<path>   |"
	@echo "| make validation_clusterizer input=<path> output=<path> |"
	@echo "| make all n=<int>                                       |"
	@echo "| make validation                                        |"
	@echo "+--------------------------------------------------------+" 
