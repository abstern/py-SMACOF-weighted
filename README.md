This is a Fortran implementation of the weighted SMACOF algorithm.
It has a Python interface which is intended to be used with f2py3.

This software is provided as-is and probably needs modification to fit in any existing workflow.

# What and why?
SMACOF is an algorithm meant to embed a decorated graph in some Euclidean space.
The embedding is constructed such that the following stress function is minimized:

![equation](https://latex.codecogs.com/svg.latex?%5Csigma_G%28X%29%20%3D%20%5Csum_%7Be%20%5Cin%20E%28G%29%7D%20w_e%20%28l_e%28X%29%20-%20l_e%28G%29%29%5E2)

Here G is a graph, E(G) are its edges, which are all decorated with a length l_e(G) and a weight w_e,
and the intent is to let the length l_e(X) of the embedded edge correspond as closely as possible to the 'exact' length l_e(G).

Many existing open implementations of SMACOF do not allow for the weights to differ from 1 or even for the lengths to differ from the unit graph distance.

One possible use of the weighted algorithm is to try to find an 'isometric' embedding of a graph on a non-Euclidean manifold into some Euclidean space. 'Isometric' is meant here in the sense of Riemannian geometry, that is, locally.

For example, consider a set of points on the sphere labeled with their great-circle distance.
A locally isometric embedding into \R^3 clearly exists, but the Euclidean distance between points disagrees more strongly with the great-circle distance as the points considered are farther apart.
Weighing the edges by their length (e.g. w_e = \exp{-l_e(G)}) allows for the nonlocal distances to stop distorting the embedding and, thus, for the algorithm to recover the sphere in the limit of sufficiently dense point sets.


# Compilation
Compilation:
```
f2py3 -c smacof.f90 -m smacof_numerical
```

Then run `python3 -i smacof.py`.

# References
Details of the algorithm can be found in Modern Multidimensional Scaling: Theory and Applications (Borg, Groenen), chapter 8.
