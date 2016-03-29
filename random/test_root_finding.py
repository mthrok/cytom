from __future__ import print_function

from time import time

import numpy as np
import scipy.optimize


def f(p, w):
    return -w + 2*p*np.log(p)/(p+1)


def fprime(p, w):
    return 2 * (1/(p+1) + np.log(p)/(p+1)**2)


def fprime2(p, w):
    return 2 * (-1/(p+1)**2 + 1/p/(p+1)**2 - 2*np.log(p)/(p+1)**3)


def test():
    r = 262144
    d = 4.5
    range_ = 4096
    cutoff = -111

    c_newton1, c_newton2, c_bisect, c_brentq, c_brenth = 0, 0, 0, 0, 0
    t_newton1, t_newton2, t_bisect, t_brentq, t_brenth = 0, 0, 0, 0, 0
    d_newton1, d_newton2, d_bisect, d_brentq, d_brenth = 0, 0, 0, 0, 0
    for w in np.linspace(0, d, 101)[1:-1]:
        print('w:{}'.format(w))
        try:
            print('... {:>15}'.format('Newton (1st):'), end=' ')
            t0 = time()
            p = scipy.optimize.newton(f, w + d, fprime, args=(w,))
            t_newton1 += time() - t0
            diff = f(p, w)
            if abs(diff) < 1e-10:
                c_newton1 += 1
                d_newton1 += diff
                print('SUCCESS: {}'.format(p))
            else:
                raise ValueError((p, diff))
        except RuntimeError:
            print('Root finding failed.')
        except ValueError as e:
            print('Incorrect value: {}'.format(e))

        try:
            print('... {:>15}'.format('Newton (2nd):'), end=' ')
            t0 = time()
            p = scipy.optimize.newton(
                f, w + d, fprime, args=(w,), fprime2=fprime2)
            t_newton2 += time() - t0
            diff = f(p, w)
            if abs(diff) < 1e-10:
                c_newton2 += 1
                d_newton2 += diff
                print('SUCCESS: {}'.format(p))
            else:
                raise ValueError((p, diff))
        except RuntimeError:
            print('Root finding failed.')
        except ValueError as e:
            print('Incorrect value: {}'.format(e))

        try:
            print('... {:>15}'.format('Bisect:'), end=' ')
            t0 = time()
            p = scipy.optimize.bisect(f, 1e-10, 2 * (w + d), args=(w,))
            t_bisect = time() - t0
            diff = f(p, w)
            if abs(diff) < 1e-10:
                c_bisect += 1
                d_bisect += diff
                print('SUCCESS: {}'.format(p))
            else:
                raise ValueError((p, diff))
        except RuntimeError:
            print('Root finding failed.')
        except ValueError as e:
            print('Incorrect value: {}'.format(e))

        try:
            print('... {:>15}'.format('BrentQ:'), end=' ')
            t0 = time()
            p = scipy.optimize.brentq(f, 1e-10, 2 * (w + d), args=(w,))
            t_brentq = time() - t0
            diff = f(p, w)
            if abs(diff) < 1e-10:
                c_brentq += 1
                d_brentq += diff
                print('SUCCESS: {}'.format(p))
            else:
                raise ValueError((p, diff))
        except RuntimeError:
            print('Root finding failed.')
        except ValueError as e:
            print('Incorrect value: {}'.format(e))

        try:
            print('... {:>15}'.format('BrentH:'), end=' ')
            t0 = time()
            p = scipy.optimize.brenth(f, 1e-10, 2 * (w + d), args=(w,))
            t_brenth = time() - t0
            diff = f(p, w)
            if abs(diff) < 1e-10:
                c_brenth += 1
                d_brenth += diff
                print('SUCCESS: {}'.format(p))
            else:
                raise ValueError((p, diff))
        except RuntimeError:
            print('Root finding failed.')
        except ValueError as e:
            print('Incorrect value: {}'.format(e))

    print()
    print('Time:')
    print('... Newton (1st): {}'.format(t_newton1/c_newton1))
    print('... Newton (2nd): {}'.format(t_newton2/c_newton2))
    print('...       Bisect: {}'.format(t_bisect/c_bisect))
    print('...       BrentQ: {}'.format(t_brentq/c_brentq))
    print('...       BrentH: {}'.format(t_brenth/c_brenth))
    print()
    print('Diff:')
    print('... Newton (1st): {}'.format(d_newton1/c_newton1))
    print('... Newton (2nd): {}'.format(d_newton2/c_newton2))
    print('...       Bisect: {}'.format(d_bisect/c_bisect))
    print('...       BrentQ: {}'.format(d_brentq/c_brentq))
    print('...       BrentH: {}'.format(d_brenth/c_brenth))

test()
