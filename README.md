# iterative5point
Fast five point epipolar geometry hypothesis generation for RANSAC using Gauss Newton on se(3)

This is an extremely rough sketch of the ideas in the paper 

* "An Iterative 5-pt Algorithm for Fast and Robust Essential Matrix Estimation (with Vincent Lui)" by Lui and Drummond.

* "Improved RANSAC performance using simple, iterative minimal-set solvers" by Rosten, Reitmayr, and Drummond.

These papers show subtle care: a couple of equations seem to break with the convention established in the rest of the work, but upon consideration they do this to maximize exploitable sparsity, allowing a careful implementer to save a couple cycles. 
