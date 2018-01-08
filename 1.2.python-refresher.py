
# coding: utf-8

# # Quick introduction to Python and Jupyter notebooks
# 
# Jupyter notebooks consist of a sequence of *cells*.  Cells can contain code (for example Python), or text.  A selected cell can be executed by pressing Shift-Enter.
# 
# ## Powers of numbers
# 
# The following cell prints the first 7 square numbers.  Execute it!

# In[ ]:

for i in range(7):
    print(i**2)


# #### Your turn!
# 
# The next cell is a copy of the previous one.  Modify it to print, say, the first 9 *cube* numbers.

# In[ ]:

for i in range(7):
    print(i**2)


# ## Prime numbers
# 
# Let’s now look at a function that checks whether a number is prime.

# In[ ]:

def is_prime_naive(n):
    for i in range(2, n):       # Note: range(n, m) -> n, n+1, ..., m-1
        if n % i == 0:          # % is the modulo operation.
           return False
    return True    


# Let’s try it out:

# In[ ]:

is_prime_naive(2**19 - 1), is_prime_naive(1000000)


# Note that the notebook simply prints the value of the last expression in a cell, we did not have to use print.
# 
# The above routine is, of course, absurdly inefficient.  Let’s try out the next Mersenne prime.

# In[ ]:

is_prime_naive(2**31 - 1)


# Since the above calculation is taking a very long time, let’s interrupt the computation by using the pull-down menu “Kernel”.
# 
# To speed up the routine, the first thing we can do is avoid to divide by even numbers:

# In[ ]:

def is_prime_slow(n):
    if n % 2 == 0:
        return False
    for i in range(3, n, 2):
        if n % i == 0:
           return False
    return True


# To better understand the built-in `range` function, place the cursor on it in the above cell and press “shift + Tab” first once, and then a second time.  This works for most functions.

# #### Your turn!
# 
# The above speedup is not significant.  To improve the runtime cost from $O(n)$ to $O(\sqrt{n})$, observe the following:
# 
# * The square root of a number `x` is given by `math.sqrt(x)`.
# * This function is part of the `math` package that needs to be imported with the `import math` statement, preferably before the function definition.
# * ‘range‘ only accepts integer arguments.  To truncate the real number `x` towards zero, use `int(x)`.
# 
# Copy the cell of `is_prime_slow` by using the “Edit” menu.  Paste it below this cell, and rename the function to `is_prime`.  Now use the above hints to make it significantly faster.

# Check it with the huge prime number below.  It should evaluate instantly.

# In[ ]:

is_prime(2**31 - 1)


# Did you bear in mind that the result of `range` does not include the upper boundary?  Check it:

# In[ ]:

is_prime(9)


# ## Linear algebra
# 
# The `numpy` package offers a wealth of functionality for multi-dimensional numerical arrays.  It gives Python MATLAB-like capabilities -- only that Python is a much nicer programming language.
# 
# Let's create a square matrix and fill both of its diagonals next to the main one with `-1`.  This is a simple model for an elastic string.
# 
# (Note that, for the sake of demonstration, we construct a *dense* matrix, even though most of its entries are zero.  Sparse linear algebra in Python is provided by the `scipy` package, for example.)

# In[ ]:

import numpy as np

n = 100
M = np.zeros((n, n))            # Create an matrix filled with zeros.
for i in range(n - 1):
    M[i, i + 1] = M[i + 1, i] = -1

print(M)


# The above code is not in the spirit of NumPy, since one of NumPy's main goals is to avoid loops by *vectorizing* code.  Vectorization can bring numerical Python up to C speed by moving loops from slow Python to fast machine code.
# 
# The below cell re-creates the same matrix `M` as the above cell.  Note that `i` is no longer a scalar, but a vector of 99 integers.

# In[ ]:

n = 100

i = np.arange(n - 1)
M = np.zeros((n, n))
M[i, i + 1] = M[i + 1, i] = -1

print(M)


# Now let's use NumPy to calculate the eigenvalues and eigenvectors of the Hermitian matrix `M`:

# In[ ]:

vals, vecs = np.linalg.eigh(M)


# The `matplotlib` package is the standard way to plot data in Python.  Let’s use it to plot the *third* eigenvector.  Note how we use NumPy’s “fancy indexing” to extract the the third (= index 2) column of the ‘vecs‘ matrix.
# 
# (The trailing semicolon has no meaning in Python.  It makes the notebook suppress the output of the preceding command.)

# In[ ]:

get_ipython().magic('run matplotlib_setup.ipy')
from matplotlib import pyplot

pyplot.plot(vecs[:, 2]);


# #### Your turn!
# As an exercise, write a loop that computes the product of all the eigenvalues of `M`.  Verify that it is equal to the determinant given by `np.linalg.det`.

# Here is a faster way to calculate the product of the entries of a NumPy array:

# In[ ]:

vals.prod()


# Note that the thing before the dot is no longer a package or module but a general object -- a NumPy array in this case.  `prod` is not a function, but a *method*, i.e. an operation of the object `vals`.
