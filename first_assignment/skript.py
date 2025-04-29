import numpy as np
import matplotlib.pyplot as plt

single_precision_forward_data = np.array([1.64393484592438, 1.64472532272339, 1.64472532272339, 1.64472532272339])
single_precision_backwards_data = np.array([1.64393448829651, 1.64493298530579, 1.6449339389801, 1.64493405818939])
double_precision_forward_data = np.array([1.64393456668156, 1.64493306684877, 1.64493396684726, 1.64493405783458])
double_precision_backwards_data = np.array([1.64393456668156, 1.64493306684873, 1.64493396684823, 1.64493405684823])

s_inf = np.pi**2/6
print(f's_inf = {s_inf}')
k = np.array([1e3, 1e6, 1e7, 1e8])

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
ax1.scatter(k, (np.abs(single_precision_forward_data - s_inf))/s_inf, label='Single Precision Forward',marker = 'x')
ax2.scatter(k, (np.abs(single_precision_backwards_data - s_inf))/s_inf, label='Single Precision Backward', marker = 'x')
ax3.scatter(k, (np.abs(double_precision_forward_data - s_inf))/s_inf, label='Double Precision Forward', marker = 'x')
ax4.scatter(k, (np.abs(double_precision_backwards_data - s_inf))/s_inf, label='Double Precision Backward', marker = 'x')
ax1.set_title('Single precision forward')
ax2.set_title('Single precision backward')
ax3.set_title('Double precision forward')
ax4.set_title('Double precision backward')
ax1.set_xlabel('k')
ax2.set_xlabel('k')
ax3.set_xlabel('k')
ax4.set_xlabel('k')
ax1.set_ylabel('Relative error log scaled')
ax2.set_ylabel('Relative error log scaled')
ax3.set_ylabel('Relative error log scaled')
ax4.set_ylabel('Relative error log scaled')
ax1.set_yscale('log')
ax2.set_yscale('log')
ax3.set_yscale('log')
ax4.set_yscale('log')
plt.suptitle('Relative error of forward and backward harmonic sum')

plt.show()