import numpy as np
import matplotlib.pyplot as plt

pi_recursion = np.loadtxt(r'/home/sasch/Studium/comp-astro/first_assignment/excercise_4/pi_recursion_data.txt',delimiter=';')
pi_kahan = np.loadtxt(r'/home/sasch/Studium/comp-astro/first_assignment/excercise_4/pi_Kahan_data.txt',delimiter=';')

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