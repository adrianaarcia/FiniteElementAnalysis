# FiniteElementAnalysis
Finite Element Analysis in MATLAB using the Direct Stiffness Method.

## Overview

`femData.m` applies the Direct Stiffness Method (DSM) to calculate deformation (i.e., the displacements of each of the nodes) of an object represented 
by a two-dimensional finite element discretization. I've included comments that outline each of the steps involved in the DSM process. 

The program plots the original shape along with the deformed shape. I included a magnification factor of 10 to make the deformed shape more evident. 
The script prints out the node with the largest displacement and the magnitude of that displacement (actual displacement, not the magnified displacement). 
It also prints out the magnitude of the maximum stress and the element where that occurs.  

The following supporting files were authored by Marshall Long. 

`femData.m` loads information on the 2D mesh from data files:
* `mesh.txt` contains the x-y coordinates in an N x 2 array; N is the number of nodes. 
* `elements.txt` will contain the element data in a M x 2 array; M is the number of elements. 
   The ith row of the elements array contains the node numbers of the ith element. 

`femData.m` also loads information on boundary conditions from files as well: 
* `fixed.txt` contains a 1D array of the nodes that are fixed (both x and y are fixed). 
* `force.txt` contains a K x 3 matrix; K is the number of nodes that have external forces 
  applied. Each row of the force array contains the node number where the force is applied, 
  the x force, and the y force. For example force = [6 0 -5000; 8 0 -5000; 14 0 -5000]; 
  would mean that nodes 6, 8 and 14, experience forces of fx = 0 and fy = -5000.
