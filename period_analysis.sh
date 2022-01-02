# Create directories
mkdir -p graphic/image
mkdir -p graphic/program-output
mkdir -p graphic/program-output-backup
mkdir -p graphic/timeseries
mkdir -p graphic/timeseries/oscillation-files
mkdir -p graphic/timeseries/timeseries-sweep

# Remove old files
rm -r graphic/timeseries/timeseries-sweep/*
rm -f graphic/program-output/*
rm -f graphic/program-output-backup/*
rm -f graphic/timeseries/oscillation-files/*

# Modify the period_analysis_configuration.py file if required
echo "Trying to use the configuration found on period_analysis_configuration.py :-)"

# Run the capture
# OBS: there are several criteria
# either (1) discard samples where the exponential decay section was not detected
#        (2) default to the entire dataset for an estimate
#        (3) ask the user to mark it manually with a simple GUI
# The current status is (1), but it can be modified at python3/kuramoto_distance_over_time.py
python3 python3/capture_period_range.py

# Plot the results
python3 python3/plot_captured_periods.py

# End :-)
rm graphic/image/wholeseries.png
rm tmp.log
echo "The script has ended! :-). The final plot is available at graphic/image/period_analysis_result.png. Subsequent applications of this program will overwrite that image. In the following section, user has the possibility of removing outliers to produce a better plot."

# Ask user if the script should continue
echo "\n\n"
read -p "Continue? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1

# Show the images and store the outliers to remove
bash ./filtering_engine.sh

# Recompute
python3 python3/plot_captured_periods.py

# End :-)
rm graphic/image/wholeseries.png
echo "The script has ended! :-). The final plot is available at graphic/image/period_analysis_result.png. Subsequent applications of this program will overwrite that image."
