
# Variables defined for the analysis of periods as a function of one of the following vars
SAMPLES=1 
TOPOLOGY=0
TOPOLOGY_NAME_DECODER = {
        0:"Ring",
        1:"Clique",
        2:"Erdos Renyi",
        3:"Small World",
        }


# Only one of the following can be a list of length>1
# as the result is converted into a 1-D plot and we use some nested looping like
# for N in NODES
#   for P in PROBA
#      ....
NNODES = [54 + i * 110 for i in range(10)]
PROBA = [0.5]
KNEIGH = [3]
J = [30]
