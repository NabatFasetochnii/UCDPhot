import numpy as np

X = np.arange(0, 10, 0.01)
X1 = 3 * X
new_x = np.zeros(len(X))
index = np.where(X < 0.05)
new_x[index] = np.power(X[index]-X1[index], 2)
print(new_x)
