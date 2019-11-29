## Description

This repository explores the theory and application of statistical machine learning, focusing on statistics, Bayesian and graphical model. It follows the structure of the classical machine learning textbooks ([Machine Learning: a probabilistic approach](https://www.amazon.com/Machine-Learning-Probabilistic-Perspective-Computation/dp/0262018020)) by Kevin Murphy.

### Topics covered in the repository
1. [Mixture models and the EM algorithm](https://github.com/jlin-vt/SML/blob/master/docs/ch11.pdf)
2. [Sparse linear models](https://github.com/jlin-vt/SML/blob/master/docs/ch13.pdf)
3. [Kernels](https://github.com/jlin-vt/SML/blob/master/docs/ch14.pdf)
4. [Gaussian Process](https://github.com/jlin-vt/SML/blob/master/docs/ch15.pdf)
5. [Adaptive basis function models](https://github.com/jlin-vt/SML/blob/master/docs/ch16.pdf)
6. [Graphical models](https://github.com/jlin-vt/SML/blob/master/docs/ch17.pdf)
7. [Variational Inference](https://github.com/jlin-vt/SML/blob/master/docs/ch21.pdf)
8. [Monte Carlo inference](https://github.com/jlin-vt/SML/blob/master/docs/ch23.pdf)
9. [Markov chain Monte Carlo (MCMC) inference](https://github.com/jlin-vt/SML/blob/master/docs/ch24.pdf)
10. [Clustering](https://github.com/jlin-vt/SML/blob/master/docs/ch25.pdf)

## Software
This toolkit contains many demos of different methods applied to many different kinds of data sets. The demos are listed [here](https://github.com/jlin-vt/SML/tree/master/R/vignettes). The vast majority of the code is written in `R`. In a few cases, I provide wrappers to implementations written in `Julia`, for speed reasons. Both programs are math based and have their own advantages: `Julia` (like `Matlab`) is the one for designing algorithms (e.g. matrix operations), while `R` is great for data analysis and statistics.

### Dependencies
If you choose `R`, you can download it from the [CRAN](https://cran.r-project.org/).  [R Studio](https://www.rstudio.com/) is an excellent graphical interface. Also, you should install a few packages:

- `mvtnorm` package computes multivariate normal and t probabilities, quantiles, random deviates and densities.

- `ggplot2`: a system for 'declaratively' creating graphics, based on "The Grammar of Graphics".

For `Julia` users, I recommend installing the latest versions of [Julia](https://julialang.org/) and [Juno](http://junolab.org/). Alternatively, you can run `Julia` on [Atom](https://atom.io/), which is another powerful editor.
Some demos may depend on the packages below:

- `Distributions`: package for probability distributions and associated functions.

##  Resources

### Print Textbooks and Online References
- [Machine Learning: a Probabilistic Perspective](http://www.cs.ubc.ca/~murphyk/MLbook/index.html)
- [Pattern Recognition and Machine Learning](https://www.amazon.com/Pattern-Recognition-Learning-Information-Statistics/dp/0387310738/ref=pd_lpo_sbs_14_img_1?_encoding=UTF8&psc=1&refRID=T2B4ZKDZR78F2J42E0FR)
- [The Elements of Statistical Learning](https://web.stanford.edu/~hastie/ElemStatLearn/)
- [All of Statistics](http://www.stat.cmu.edu/~larry/all-of-statistics/index.html)
- [CS142: Machine Learning](http://cs.brown.edu/courses/cs142/index.html)

### Programming Languages
- [Learning Julia](http://julialang.org/learning/)
- [Advanced R](http://adv-r.had.co.nz/)

## Disclaimer
The above material is purely for self-study purposes and will not be published in any form or by any means, including photocopying, recording, or other electronic or mechanical methods, without the prior written permission of original authors.
