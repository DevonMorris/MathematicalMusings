import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

phi = norm.cdf
z = np.linspace(-5, 5, num=1000)
d = np.linspace(0, 4, num=9)

for d_ in d:
    P_FA = 1 - phi(z)
    P_D = 1 - phi(z - d_)

    plt.plot(P_FA,P_D, label="$d = {}$".format(d_))

plt.title("ROC Curves")
plt.xlabel("$P_{FA}$")
plt.xlabel("$P_D$")
plt.legend(loc="best")
plt.show()

d = 0
confident = False

while not confident:
    P_FA = 1 - phi(z)
    P_D = 1 - phi(z - d)
    t1 = P_FA < 0.01
    t2 = P_D > 0.99
    p = t1*t2
    confident = np.max(p)
    d += 0.01

print(d)
