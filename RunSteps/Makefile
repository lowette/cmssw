all:
	cmsRun Generation/python/Generation.py _n $(n)
	cmsRun Digitizer/python/Digitizer.py
	cmsRun Clusterizer/python/Clusterizer.py
	cmsRun Digitizer/test/RunDigiValidation.py
	cmsRun L1Trigger/test/RunL1TriggerValidation.py
	cmsRun ValidaThor/python/RunClusterValidation.py
	cmsRun ValidaThor/python/RunRecHitValidation.py

generation:
	cmsRun Generation/python/Generation.py _n $(n) _output $(output)

digitizer:
	cmsRun Digitizer/python/Digitizer.py _input $(input) _output $(output)

old_digitizer:
	cmsRun Digitizer/python/OldDigitizer.py _input $(input) _output $(output)

clusterizer:
	cmsRun Clusterizer/python/Clusterizer.py _input $(input) _output $(output)

validation_digitizer:
	cmsRun Digitizer/test/RunDigiValidation.py _input $(input) _output $(output)

validation_l1trigger:
	cmsRun L1Trigger/test/RunL1TriggerValidation.py _input $(input) _output $(output)

validation_clusterizer:
	cmsRun ValidaThor/python/RunClusterValidation.py _input $(input) _output $(output)

validation_rechits:
	cmsRun ValidaThor/python/RunRecHitValidation.py _input $(input) _output $(output)

validation: validation_digitizer validation_clusterizer validation_l1trigger validation_rechits

configure:
	mkdir -p Output
	scram b -j 16

help:
	@echo "+--------------------------------------------------------+"
	@echo "|                        RunSteps                        |"
	@echo "+--------------------------------------------------------+"
	@echo "| List of available commands:                            |"
	@echo "+--------------------------------------------------------+"
	@echo "| make all n=<int>                                       |"
	@echo "| make generation n=<int> output=<path>                  |"
	@echo "| make digitizer input=<path> output=<path>              |"
	@echo "| make old_digitizer input=<path> output=<path>          |"
	@echo "| make clusterizer input=<path> output=<path>            |"
	@echo "| make validation_digitizer input=<path> output=<path>   |"
	@echo "| make validation_l1trigger input=<path> output=<path>   |"
	@echo "| make validation_clusterizer input=<path> output=<path> |"
	@echo "| make validation_rechits input=<path> output=<path>     |"
	@echo "| make validation                                        |"
	@echo "| make configure                                         |"
	@echo "| make help                                              |"
	@echo "+--------------------------------------------------------+"
