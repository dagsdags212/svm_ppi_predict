# SVM for Predicting PPIs

Goals:
1. Implement the **conjoint triad** descriptor
2. Implement the kernel function
3. Train the model using SVM

## Classification of Amino Acids

|Category|Amino Acids|
|---|---------|
|1|A, G, V|
|2|I, L, F, P|
|3|Y, M, T, S|
|4|H, N, Q, W|
|5|R, K|
|6|D, E|
|7|C|


## Conjoint Triad Method

- maps a protein sequence into a binary space (V, F) where V is the number of features and F is the frequency of each triad type
- using 7 features, the length of V is 7x7x7=343 to represent triads
- length of F is equal to the length of V
- values of vector F are normalize to be length-invariant:

$$ d_i = \frac{(f_i - min\{f_1, f_2, ..., f_{343}\})}{max\{f_1, f_2, ..., f_{343}\}} $$


## References

- Paper by [Shen _et al._ (2007)](https://www.pnas.org/doi/10.1073/pnas.0607879104)
