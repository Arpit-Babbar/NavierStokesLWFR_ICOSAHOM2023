#timestep  time            dt               Viscous_drag_corr   Viscous_lift_corr   Pressure_drag    Pressure_lift
from mpl import *
import numpy as np
import matplotlib.pyplot as plt
import os

data = np.loadtxt("analysis.dat")

t = data[:,1]
drag_f = data[:,3]
lift_f = data[:,4]
drag_p = data[:,5]
lift_p = data[:,6]

drag = drag_p - drag_f
lift = lift_p - lift_f

plt.figure()
plt.plot(t, drag)
plt.xlabel("$t$")
plt.ylabel("$c_d$")
plt.ylim((3.075,3.135))
plt.xlim((5.055,9.1))
plt.grid()
plt.savefig("cd.pdf")
os.system("bash pdfbb cd.pdf")

plt.figure()
plt.plot(t, lift)
plt.xlabel("$t$")
plt.ylabel("$c_l$")
plt.ylim((-1,1))
plt.xlim((5.055,9.1))
plt.grid()
plt.savefig("cl.pdf")
os.system("bash pdfbb cl.pdf")
plt.show()
