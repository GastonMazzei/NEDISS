# Create directories
mkdir -p graphic/program-output
mkdir -p graphic/program-output-backup
mkdir -p graphic/timeseries
mkdir -p graphic/timeseries/oscillation-files
mkdir -p graphic/timeseries/timeseries-sweep
mkdir -p graphic/image

# Remove old files
rm -f graphic/program-output/*
rm -f graphic/program-output-backup/*
rm -f graphic/timeseries/oscillation-files/*

# Load variables
source simulation.conf.sh

# Simulate
mpirun -oversubscribe  -x SAMPLING_FREQ=$SAMPLING_FREQ -x EQNUMBER=$EQNUMBER -x SOLVER=$SOLVER -x SOLVERDEPTH=$SOLVERDEPTH -x K1=$K1 -x K2=$K2 -x K3=$K3 -x K4=$K4  -x TOPOLOGY=$TOPOLOGY -x J=$J -x WMIN=$WMIN -x WMAX=$WMAX -x kneigh=$kneigh -x proba=$proba -x NRUNS=$NRUNS -x NNODES=$NNODES -x TEST=-1 -x SEED=$SEED -x OMP_THREAD_LIMIT=$OMP_THREAD_LIMIT  -x  OMP_NESTED=true -np $NP  cmake-build-debug/cppprojct

# Do a backup just in case
echo "Doing backup"
cp graphic/program-output/* graphic/program-output-backup/.

# Extract timeseries
echo "Extracting timeseries"
python3 python3/extract_timeseries_from_graphviz.py

# More Python preprocessing
echo "Plotting distance over time"
python3 python3/kuramoto_distance_over_time.py

echo "The script has ended! :-) the configuration and result are appended to graphic/timeseries/total-results.csv and an illustration of the fit is available at graphic/image/maxdistovertime. The data  is available at graphic/timeseries/result.pkl. The image and the data will be overwritten by subsequent calls to this program, but the results will be accumulated in the .csv. Calls to other programs can wipe the .csv, as for example does the period_analysis.sh before starting."
