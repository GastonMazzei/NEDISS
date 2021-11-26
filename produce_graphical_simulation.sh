

source simulation.conf.sh

mpirun -oversubscribe  -x SAMPLING_FREQ=$SAMPLING_FREQ -x EQNUMBER=$EQNUMBER -x SOLVER=$SOLVER -x SOLVERDEPTH=$SOLVERDEPTH -x K1=$K1 -x K2=$K2 -x K3=$K3 -x K4=$K4  -x TOPOLOGY=$TOPOLOGY -x kneigh=$kneigh -x proba=$proba -x NRUNS=$NRUNS -x NNODES=$NNODES -x TEST=-1 -x SEED=$SEED -x OMP_THREAD_LIMIT=$OMP_THREAD_LIMIT  -x  OMP_NESTED=true -np $NP  cmake-build-debug/cppprojct

# give colors to data
python3 python3/make_undirected.py

# render data
SAMPLES=$((NRUNS / SAMPLING_FREQ))
for ((i=0;i<SAMPLES;++i)); do
	dot -Kcirco -Tpng "graphic/processed-program-output/graphviz.$i.dot" -o "graphic/rendered-program-output/graphviz.$i.png";
done

python3 python3/extract_timeseries_from_graphviz.py
rm graphic/program-output-backup/*
mv graphic/program-output/* graphic/program-output-backup/.
rm graphic/video/result.avi
rm graphic/video/result-oscillating.avi
python3 python3/correct_names.py
python3 python3/plot_oscillating.py
rm graphic/processed-program-output/*
rm graphic/rendered-program-output/*
rm graphic/timeseries/oscillation-files/*


