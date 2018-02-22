# iterative5point
Fast five point epipolar geometry hypothesis generation for RANSAC using Gauss Newton on the Lie algebra se(3)

This particular implementation is not fast, but the technique is. It is a quite rough sketch of the ideas in the papers

* "An Iterative 5-pt Algorithm for Fast and Robust Essential Matrix Estimation (with Vincent Lui)" by Lui and Drummond.

* "Improved RANSAC performance using simple, iterative minimal-set solvers" by Rosten, Reitmayr, and Drummond.

These papers show subtle care: a couple of equations seem to break with the convention established in the rest of the work, but upon consideration they do this to maximize exploitable sparsity, allowing a careful implementer to save a couple cycles. 

Anyone looking to accelerate visual odometry or just avoid some needless SVD's is encouraged to review the papers and use this as a first rough guide.

For an introduction to Lie algebras, please see Drummond's notes

* http://twd20g.blogspot.com/p/notes-on-lie-groups.html
