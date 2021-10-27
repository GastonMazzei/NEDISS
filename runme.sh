echo "arg is $1"
(cd cmake-build-debug ;  ./cppprojct $1 )
python3 plot.py
