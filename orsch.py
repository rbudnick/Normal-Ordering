from math import sqrt,exp,pi
from random import gauss

def order_prob(means, vars):
    """Probability that normal distributions are in index order.
    
    Arguments:
    means -- the means of each normal random variable
    vars -- the variances of each normal random variable
    """
    # Find difference variable means
    # (Pr[X>Y] = P[X-Y>0])
    mus = [j - i for i, j in zip(means[: -1], means[1 :])]
    m = len(mus)

    # Find Cholensky decomposition of covariance matrix
    diags = [] # Main diagonal entries of the matrix
    off_diags = [0] # Off-diagonal entries, with extra b_{1,0}=0
    for i in range(m):
        diags.append(sqrt((vars[i]+vars[i+1])-off_diags[-1]*off_diags[-1]))
        if i < m-1:
            off_diags.append(-vars[i+1]/diags[-1])

    # grid is the set of points on which function values are calculated
    # This is the main way to tune accuracy of the result
    grid = list(map(lambda x: (x-5000)/500., range(10001)))
    # You may want to increase density of the points near zero, as below
    # grid = list(map(lambda x: (x-5000)*abs(x-5000)/2500000., range(10001)))

    normal = list(map(lambda x: 1/sqrt(2*pi)*exp(-x*x/2), grid))

    # Calculate recursive integrals
    f = [1]*len(grid)
    for i in range(m):
        g = [0]*len(grid) # New function being constructed

        # offset and coefficient used to find lower bound of integral piece
        offset = -mus[m-1-i]/diags[m-1-i]
        coefficient = -off_diags[m-1-i]/diags[m-1-i]

        lower, upper = grid[-1], grid[-1] # Bounds of integral piece
        j = len(grid)-2 # Index of grid block that we've summed to
        f_top = f[-1]*normal[-1] # Function value at top of integral piece
        partial = 0 # Portion of trapezoid included in previous piece
        for z in range(len(grid)-1,-1,-1): # Integrate from above
            upper = lower
            lower = offset + coefficient * grid[z]
            # Skip if we're not in range yet
            if lower > grid[-1]:
                lower = grid[-1]
                g[z-1] = g[z]
                continue
            # Add all of the whole trapezoids
            while lower < grid[j] and j >= 0:
                f_bottom = f[j] * normal[j]
                g[z] += (f_top + f_bottom)/2 * (grid[j+1]-grid[j])
                f_top = f_bottom
                j -= 1
            if j >= 0:
                # Add the bottom partial-trapezoid
                prop = (grid[j+1]-lower)/(grid[j+1]-grid[j])
                f_bottom = prop*(f[j] * normal[j]) + (1-prop)*f_top
                partial = (f_top + f_bottom)/2 * (grid[j+1]-lower)
                g[z] += partial
            else:
                # Or, if bottomed-out
                for k in range(z,0,-1):
                    g[k-1] = g[k]
                break
            if z > 0:
                g[z-1] = g[z] - partial
        f = g
    
    return g[0]

def sample(means, vars, res=1000000):
    """Sample probability that normal distributions are in index order.
    
    Arguments:
    means -- the means of each normal random variable
    vars -- the variances of each normal random variable
    res -- the number of samples to draw
    """
    counter = 0
    for i in range(res):
        values = []
        for j in range(len(means)):
            values.append(gauss(means[j],sqrt(vars[j])))
        if sorted(values) == values:
            counter += 1
    return (counter/res)

def main():
    my_means = [0,2,1,3]
    my_vars = [1,3,2,2]
    print("I have " + str(len(my_means)) +
          " independent normally-distributed random variables")
    for i in range(len(my_means)):
        print("N_" + str(i) + " has mean " + str(my_means[i]) +
              " and variance " + str(my_vars[i]))
    print("What is the probability that N_0", end="")
    for i in range(1,len(my_means)):
        print(" < N_" + str(i), end="")
    print("?")
    print("Computing integrals following Miwa et al. (2003) gives " +
          str(order_prob(my_means,my_vars)))
    my_res = 100000
    print("Sampling the variables " + str(my_res) +
          " times gives " + str(sample(my_means,my_vars,my_res)))

if __name__ == "__main__":
    main()