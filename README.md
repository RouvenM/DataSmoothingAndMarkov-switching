# On-the-combination-of-data-smoothing-and-Markov-switching-models
This repository contains R code to provide evidence that smoothing creates autocorrelation which has to be kept in mind when applying, e.g., Markov-switching models.

In particular, we found a paper by Sandri et al. (2020) in which data has been smoothed and afterwards fit by a Markov-switching model. In this regard, we demonstrate that such a Markov-switching model can be fit even to i.i.d. data.
To show this, we replicate the pre-processing steps with the data used and with simulated i.i.d. data. We thereby demonstrate that the approach as suggested by the authors does in fact (erroneously) identify state-switching dynamics even in i.i.d. data.

References:
M. Sandri, P. Zuccolotto, and M. Manisera. Markov switching
modelling of shooting performance variability and teammate
interactions in basketball. Journal of the Royal Statistical
Society: Series C, 69(5):1337–1356, 2020.
