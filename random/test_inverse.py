"""Script for testing if we can approximate inverse biexponential function"""
from __future__ import division

import numpy as np
from numpy.polynomial.polynomial import polyval
from matplotlib import pyplot as plt


def biexp(x, a=0.5, b=1, c=0.5, d=1, f=0):
    """Biexponential transform. Default to sinh"""
    return a * np.exp(b * x) - c * np.exp(- d * x) + f


def comp_series(a, b, c, d, N=7):
    ret = []
    term1, term2, denom = a, c, 1
    print("a={}, b={}, c={}, d={}".format(a, b, c, d))
    for n in range(1, N+1):
        denom *= n
        term1 *= b
        term2 *= -d
        ret.append((term1 - term2)/denom)
        print(n, term1, term2, denom)
    return ret


def comp_series_reversion(an):
    a1, a2, a3, a4, a5, a6, a7 = an[:7]
    A1 = a1**(-1)
    A2 = a1**(-3) * (-a2)
    A3 = a1**(-5) * (2*(a2**2) - a1*a3)
    A4 = a1**(-7) * (5*a1*a2*a3 - (a1**2)*a4 - 5*(a2**3))
    A5 = a1**(-9) * (
        6*(a1**2)*a2*a4 + 3*(a1**2)*(a3**2) +
        14*(a2**4) - (a1**3)*a5 - 21*a1*(a2**2)*a3
    )
    A6 = a1**(-11) * (
        7*(a1**3)*a2*a5 + 7*(a1**3)*a3*a4 + 84*a1*(a2**3)*a3 - (a1**4)*a6 -
        28*(a1**2)*a2*(a3**2) - 42*(a2**5) - 28*(a1**2)*(a2**2)*a4
    )
    A7 = a1**(-13) * (
        8*(a1**4)*a2*a6 + 8*(a1**4)*a3*a5 + 4*(a1**4)*(a4**2) +
        120*(a1**2)*(a2**3)*a4 + 180*(a1**2)*(a2**2)*(a3**2) + 132*(a2**6) -
        (a1**5)*a7 - 36*(a1**3)*(a2**2)*a5 - 72*(a1**3)*a2*a3*a4 -
        12*a1**3*(a3**3) - 330*a1*(a2**4)*a3
    )
    An = [A1, A2, A3, A4, A5, A6, A7]
    return An


def inv_biexp(y, a=0.5, b=1, c=0.5, d=1, f=0):
    an = comp_series(a, b, c, d)
    An = comp_series_reversion(an)
    x = polyval(y, [a-c] + An[:3])
    return x


def test(a=0.5, b=1, c=0.5, d=1, f=0):
    x0 = np.linspace(-1.5, 1.5, num=100)
    y0 = biexp(x0, a, b, c, d)
    y1 = np.linspace(min(y0), max(y0), num=100)
    x1 = inv_biexp(y1, a, b, c, d)

    fig = plt.figure()
    # ax = fig.add_subplot(2, 1, 1)
    plt.scatter(x0, y0, color='r')
    # ax = fig.add_subplot(2, 1, 2)
    plt.scatter(x1, y1, color='b')
    return np.sum(np.power(x0 - inv_biexp(y0, a, b, c, d), 2))


def check_sinh():
    print(test())


def check_biexp():
    print(test(a=1, b=1, c=1, d=1, f=0))


def main():
    check_sinh()
    check_biexp()
    plt.show()

if __name__ == '__main__':
    main()
