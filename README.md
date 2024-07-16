# Normal Ordering Repository

This repository includes a single main function `order_prob` that can be used to calculate the probability that the sampled values of a set of independent normally-distributed random variables with provided means and variances are found in a particular order. The algorithm used here is a special case of the procedure described in Tetsuhisa Miwa, A. J. Hayter, and Satoshi Kuriki's article, "The Evaluation of General Non-Centred Orthant Probabilities", published in the *Journal of the Royal Statistical Society*, Volume 65, Number 1, pages 223-234 in 2003.

The file also includes a function `sample` for finding the same probability by sampling from the random variables many times, which can be used to confirm the accuracy of the main function.

## Requirements

`order_prob` uses only imports from the standard library `math` module; `sample` uses the standard library `random` as well. The code should work for all versions of Python 3, and I believe many versions of Python 2 as well.

## Usage

To use the function, simply import it into your Python script. You must pass it two arrays of equal length, the first of which contains the means of the random variables (which may be any numerical value), and the second of which contains the variances of the random variables (which must be positive). The function calculates the probability that the random variables are sampled in index order: that the sample from the variable with mean and variance given by the first elements of the arrays is less than the sample from the variable with mean and variance given by the second elements of the arrays, and so on.

```python
from orsch import order_prob

print(str(order_prob([-1,1,0,2],[1,3,2,2])))
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.