import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

# def opticalvignet_correctarr(params=[-1.46961740e-19,  6.24888387e-16, -5.33987703e-13, -1.78593045e-10, -2.51044719e-06, 
#                                              6.66156473e-10, -3.55316143e-06,  5.06805490e-02, -8.54433272e+01, -3.26824001e+05]):

# def quartic(x, a, b, c, d, e):
#     y = a*x**4 + b*x**3 + c*x**2 + d*x + e
#     return y

# def fit_3d(Ax, ax, bx, cx, dx, ex, ay, by, cy, dy, ey):
#     x, y = Ax
#     z = quartic(x, ax, bx, cx, dx, ex)*quartic(y, ay, by, cy, dy, ey)
#     return z

def cos4_3d(Ax, N, A, x_0, y_0):
    x, y = Ax
    z = N*(np.cos(A * np.sqrt((x - x_0)**2 + (y - y_0)**2)))**4
    return z

def opticalvignet_correctarr(params=[9.98828800e-01, -2.52880385e-04,  9.21335131e+02,  9.17392750e+02]):
    fit_x = np.linspace(0, 2048, 2048)
    fit_y = np.linspace(0, 2048, 2048)
    X, Y = np.meshgrid(fit_x, fit_y)
    fit_z = cos4_3d((X, Y), *params)
    return fit_z

def quartic(x, a, b, c, d, e):
    y = a*x**4 + b*x**3 + c*x**2 + d*x + e
    return y

def fit_3d(Ax, ax, bx, cx, dx, ex, ay, by, cy, dy, ey):
    x, y = Ax
    z = quartic(x, ax, bx, cx, dx, ex)*quartic(y, ay, by, cy, dy, ey)
    return z

def vignet_corr(Subts):
    corrarr = opticalvignet_correctarr()
    CorrectedSubts = [Subts[i]/(corrarr/np.amax(corrarr)) for i in range(len(Subts))]
    return CorrectedSubts

if __name__=="__main__":
    correctarray = opticalvignet_correctarr()
    plt.imshow(correctarray)
    plt.show()