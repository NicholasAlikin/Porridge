import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    filename = r'E:\projects\Porridge\cpp\build\res.txt'

    res = np.loadtxt(filename, comments='#').T
    if 0:
        fig, ax = plt.subplots()

        k = 1.; m = 1.; d = 0.1; f0 = 1.
        A = lambda w: (k-m*w**2)*f0/( (k-m*w**2)**2 + (w*d)**2)
        B = lambda w: w*d*f0/( (k-m*w**2)**2 + (w*d)**2)
        C = lambda w: np.sqrt(A(w)**2 + B(w)**2)

        w = np.linspace(0,3,100)
        ax.plot(w, C(w), '-b', lw=2)
        w = res[1]
        A = res[0]
        if res.shape[0]>2:
            Y = res[2:]
        ax.plot(w, A, '.r', markersize=6)
          

  
    if 1:
        ndof = res.shape[0]-1
        fig, ax = plt.subplots(ncols=ndof)
        w = res[-1]
        
        for k in range(ndof):

            ax[k].plot(w,res[k],'.')
            ax[k].legend(str(k+1))
        # ax[1].plot(w,res[1],'.')
        # ax[2].plot(w,res[2],'.')
        # ax[3].plot(w,res[3],'.')
        # ax[4].plot(w,res[4],'.')
        
        # ax[0].legend('ux')
        # ax[1].legend('uy')
        # ax[2].legend('uy')
        # ax[3].legend('ty')
        # ax[4].legend('tz')
    
    plt.show()
