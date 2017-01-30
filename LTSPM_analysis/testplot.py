import numpy as np
import matplotlib.pyplot as plt

fig1, axes1 = plt.subplots()
axes1.plot(range(10))
#figure one

fig2, axes2 = plt.subplots()
axes2.plot(range(3))
#figure two

plt.get_current_fig_manager().window.setGeometry(600,400,1000,800)

plt.show()