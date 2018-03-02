import pickle
import matplotlib.pyplot as plt

eps = {}
eps["0.0001"] = pickle.load(open("eps0.0001.p"))
eps["0.001"] = pickle.load(open("eps0.001.p"))
eps["0.01"] = pickle.load(open("eps0.01.p"))
plt.plot(eps['0.01'],'g-', label="eps = 0.01")
plt.plot(eps['0.001'],'b-', label="eps = 0.001")
plt.plot(eps['0.0001'],'r-', label="eps = 0.0001")



plt.xlabel("$t$",fontsize=24)
plt.ylabel("$L_2$",fontsize=24)
plt.legend()
plt.savefig("L2error_vs_eps.png", bbox_inches="tight")

