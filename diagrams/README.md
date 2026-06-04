# CC Diagrams

This library contains ways to generate coupled-cluster diagrams to arbitrary order. It also comes with a tester program which will print out the equations for the diagrams, using the assumption that the indices
on the amplitudes are in strictly increasing order. Permutation operators are ignored, since every index will need to be permuted in some way due to the ordering restriction.
Internally, the diagrams are represented in a similar way as is used by Kállay and Surján[^c1].

[^c1]: Kállay, M.; Surján, P. R. Higher Excitations in Coupled-Cluster Theory. The Journal of Chemical Physics 2001, 115 (7), 2945–2954. [https://doi.org/10.1063/1.1383290]().


