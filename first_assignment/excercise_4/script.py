import numpy as np
import matplotlib.pyplot as plt
import os

def rel_path(filename):
    """
    Given a filename, returns the absolute path to the file relative to the caller's location.

    Args:
        filename (str): The name of the file.

    Returns:
        str: The absolute path to the file relative to the caller's location.
    """
    caller_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(caller_dir, filename)

pi_recursion = np.loadtxt(rel_path('pi_recursion_data.txt'),delimiter=';')
pi_kahan = np.loadtxt(rel_path('pi_Kahan_data.txt'),delimiter=';')

plt.figure(figsize=(10, 5))
plt.scatter(pi_recursion[:,0],abs( pi_recursion[:,1] - np.pi), label='Pi recursion', color='blue',marker='x')
plt.scatter(pi_kahan[:,0], abs(pi_kahan[:,1] - np.pi), label='Pi Kahan', color='red', marker = 'x')
plt.yscale('log')
plt.grid()
plt.title('Pi recursion vs Pi Kahan')
plt.xlabel('Number of iterations')
plt.ylabel(r'$|A_n - \pi|$')
plt.legend()
plt.show()