# Create directories
mkdir graphic/program-output
mkdir graphic/processed-program-output
mkdir graphic/rendered-program-output
mkdir graphic/program-output-backup
mkdir graphic/timeseries
mkdir graphic/timeseries/oscillation-files
mkdir graphic/video

# Remove old files
rm graphic/program-output-backup/*
rm graphic/video/result.avi
rm graphic/video/result-oscillating.avi
rm graphic/processed-program-output/*
rm graphic/rendered-program-output/*
rm graphic/timeseries/oscillation-files/*

# Load variables
source simulation.conf.sh

# Simulate
mpirun -oversubscribe  -x SAMPLING_FREQ=$SAMPLING_FREQ -x EQNUMBER=$EQNUMBER -x SOLVER=$SOLVER -x SOLVERDEPTH=$SOLVERDEPTH -x K1=$K1 -x K2=$K2 -x K3=$K3 -x K4=$K4  -x TOPOLOGY=$TOPOLOGY -x J=$J -x WMIN=$WMIN -x WMAX=$WMAX -x kneigh=$kneigh -x proba=$proba -x NRUNS=$NRUNS -x NNODES=$NNODES -x TEST=-1 -x SEED=$SEED -x OMP_THREAD_LIMIT=$OMP_THREAD_LIMIT  -x  OMP_NESTED=true -np $NP  cmake-build-debug/cppprojct

# Do a backup just in case
cp graphic/program-output/* graphic/program-output-backup/.

# Color the data according to the value :-)
python3 python3/make_undirected.py

# Render data 
SAMPLES=$((NRUNS / SAMPLING_FREQ))
MAX_N_CIRCOENGINE=550
CIRCO_CAN_WORK=true
TIMETOL=3 #in seconds
for ((i=0;i<SAMPLES;++i)); do
        echo "Starting lap $i"
        if [ $TOPOLOGY == 0 ] && [ $CIRCO_CAN_WORK ]; then
                if [ $i -eq 0 ]; then echo "About to try Circo"; fi
                timeout $TIMETOL dot -Kcirco -Tpng "graphic/processed-program-output/graphviz.$i.dot" -o "graphic/rendered-program-output/graphviz.$i.png";
                EXIT_STATUS=$?
                if [ $EXIT_STATUS -eq 124 ]; then
                        CIRCO_CAN_WORK=false;
                        echo "The circo engine won't work... defaulting to patchwork :-)"
                        dot -Kpatchwork -Tpng "graphic/processed-program-output/graphviz.$i.dot" -o "graphic/rendered-program-output/graphviz.$i.png";
                fi
        elif [ $TOPOLOGY == 1 ] && [ $MAX_N_CIRCOENGINE -gt $NNODES  ] && [ $CIRCO_CAN_WORK ]; then
                if [ $i -eq 0 ]; then echo "About to try Circo"; fi
                timeout $TIMETOL dot -Kcirco -Tpng "graphic/processed-program-output/graphviz.$i.dot" -o "graphic/rendered-program-output/graphviz.$i.png";
                EXIT_STATUS=$?
                if [ $EXIT_STATUS -eq 124 ]; then
                        CIRCO_CAN_WORK=false;
                        echo "The circo engine won't work... defaulting to patchwork :-)"
                        dot -Kpatchwork -Tpng "graphic/processed-program-output/graphviz.$i.dot" -o "graphic/rendered-program-output/graphviz.$i.png";
                fi
        elif [ $TOPOLOGY == 3 ] && [ $MAX_N_CIRCOENGINE -gt $NNODES ] && [ $CIRCO_CAN_WORK ]; then
                if [ $i -eq 0 ]; then echo "About to try Circo"; fi
                timeout $TIMETOL dot -Kcirco -Tpng "graphic/processed-program-output/graphviz.$i.dot" -o "graphic/rendered-program-output/graphviz.$i.png";
                EXIT_STATUS=$?
                if [ $EXIT_STATUS -eq 124 ]; then
                        CIRCO_CAN_WORK=false;
                        echo "The circo engine will take too long... defaulting to patchwork :-)"
                        dot -Kpatchwork -Tpng "graphic/processed-program-output/graphviz.$i.dot" -o "graphic/rendered-program-output/graphviz.$i.png";
                fi
        else
                dot -Kpatchwork -Tpng "graphic/processed-program-output/graphviz.$i.dot" -o "graphic/rendered-program-output/graphviz.$i.png";
        fi
done


# Extract timeseries
python3 python3/extract_timeseries_from_graphviz.py


# More Python preprocessing
python3 python3/correct_names.py
python3 python3/plot_oscillating.py
python3 python3/kuramoto_distance_over_time.py


# Final section
./make_youtube_video.sh

