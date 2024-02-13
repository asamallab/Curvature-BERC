#!/bin/bash
# Compute Forman, Augmented Forman and Ollivier curvature for undirected and unweighted networks of global indices 
## g++ FormanUndirected.cpp -o FormanCurvature
## g++ FormanTriangleUndirected.cpp -o FormanTriangleUndirected
#
# Terminal command 
# ./Compute_Curvatures.sh BA 1000_3
# 
# Input 1: ModelName 
# Input 2: numberOfNodes_Parameter
# -----------------------------------------------------------
# Written by: Madhumita Mondal
#------------------------------------------------------------

n=1

if [ ! -d "../Data/Model_Networks/$1/$2/Forman" ]; then
    mkdir -p "../Data/Model_Networks/$1/$2/Forman"
fi
if [ ! -d "../Data/Model_Networks/$1/$2/AugForman" ]; then
    mkdir -p "../Data/Model_Networks/$1/$2/AugForman"
fi
if [ ! -d "../Data/Model_Networks/$1/$2/Ollivier" ]; then
    mkdir -p "../Data/Model_Networks/$1/$2/Ollivier"
fi

seed=1

while [ $n -le $seed ]; do
	echo $n
	./FormanCurvature 0 "../Data/Model_Networks/$1/$2/$1_node_seed_$n.txt" 0 "../Data/Model_Networks/$1/$2/$1_edge_seed_$n.txt" "../Data/Model_Networks/$1/$2/Forman/$1_edge_seed_$n.txt" "../Data/Model_Networks/$1/$2/Forman/$1_node_seed_$n.txt"

	./FormanTriangleUndirected 0 "../Data/Model_Networks/$1/$2/$1_node_seed_$n.txt" 0 "../Data/Model_Networks/$1/$2/$1_edge_seed_$n.txt" "../Data/Model_Networks/$1/$2/AugForman/$1_edge_seed_$n.txt" "../Data/Model_Networks/$1/$2/AugForman/$1_node_seed_$n.txt"
    
	python3 ./OR_NotNorm.py 0 "../Data/Model_Networks/$1/$2/$1_edge_seed_$n.txt" "../Data/Model_Networks/$1/$2/Ollivier/$1_edge_seed_$n.txt" "../Data/Model_Networks/$1/$2/Ollivier/$1_node_seed_$n.txt"

	n=$((n+1))
done 

echo "Done"