#!/bin/bash
# Compute HGG for networks of global indices
#
# To run this script "HyperbolicGraphGenerator-1.0.3" folder is needed, which can be found: https://github.com/mwarning/Hyperbolic-Graph-Generator 
## Terminal comand
# ./Compute_HGG.sh 1000 3
# 
# Input 1: number_of_nodes 
# Input 2: targeted_ave_degree
#----------------------------------------------------------------------------------------------------------------------


no=$1
k=$2
g=2
n=1
seed=1

if [ ! -d "../Data/Model_Networks/HGG/"$no"_$k" ]; then
    mkdir -p "../Data/Model_Networks/HGG/"$no"_$k"
fi
while [ $n -le $seed ]; do
	echo $n $no $k
	HyperbolicGraphGenerator-1.0.3/tools/hyperbolic_graph_generator -n $no -k $k -g $g -f "../Data/Model_networks/HGG/"$no"_$k/HGG_seed_"$n
	n=$((n+1))
done

echo 'Done'
