import matplotlib.pyplot as plt
import numpy as np

data = np.genfromtxt("../resources/error.log", delimiter=",")
time = data[:, 0]
iter = data[:, 1]
error = data[:, 2]

ref = error[0] / (iter[1] ** 2)
second_order = ref * (iter ** 2)

plt.figure(figsize=(12, 6))
#plt.plot(time, error, "r-", label="Error over time")
plt.plot(iter, error, "g-", label="Error over iterations")
plt.plot(iter, second_order, "b-", label="Second order")

plt.yscale("log")
plt.xscale("log")

plt.xlabel("Time steps/Iterations")
plt.ylabel("Error")
plt.legend()
plt.grid(True)

plt.savefig("../resources/error_plot.png")
