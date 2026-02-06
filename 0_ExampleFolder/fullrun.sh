#!/bin/bash

cd 1_MAPPING
bash mapping1.sh | tee MAPPING1.log
bash mapping2.sh | tee MAPPING2.log
cd ..
cd 2_TRACKING
bash tracking1.sh | tee TRACKING1.log
bash tracking2.sh | tee TRACKING2.log
bash tracking3.sh | tee TRACKING3.log
cd ..
cd 3_GROWTHSIMULATION
bash growthSim1.sh | GROWTHSIM1.log
bash growthSim2.sh | GROWTHSIM1.log
