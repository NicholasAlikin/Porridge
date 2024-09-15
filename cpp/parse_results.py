import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    filename = r'E:\projects\Porridge\cpp\math\build\res.txt'

    res = np.loadtxt(filename, comments='F')

    fig, ax = plt.subplots()

    k = 1.; m = 1.; d = 0.1; f0 = 1.
    A = lambda w: (k-m*w**2)*f0/( (k-m*w**2)**2 + (w*d)**2)
    B = lambda w: w*d*f0/( (k-m*w**2)**2 + (w*d)**2)
    C = lambda w: np.sqrt(A(w)**2 + B(w)**2)

    w = np.linspace(0,3,100)
    ax.plot(w, C(w), '-b', lw=2)
    ax.plot(res[:,1], res[:,0], '.r', markersize=6)

    plt.show()