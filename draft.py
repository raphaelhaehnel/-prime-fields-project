import numpy as np
import matplotlib.pyplot as plt

# Define the mean and covariance matrix
mu = np.array([0, 0])
cov = np.array([[1, 0], [0, 1]])

# Create a meshgrid of x and y values
x, y = np.meshgrid(np.linspace(-3, 3, 100), np.linspace(-3, 3, 100))

# Create a 2D Gaussian distribution
pos = np.empty(x.shape + (2,))
pos[:, :, 0] = x
pos[:, :, 1] = y
z = np.exp(-0.5 * np.sum((pos - mu) @ np.linalg.inv(cov) * (pos - mu),
           axis=2)) / (2 * np.pi * np.sqrt(np.linalg.det(cov)))

# Plot the distribution
plt.contourf(x, y, z)
plt.colorbar()
plt.show()
