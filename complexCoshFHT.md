
# Numerical experiments on the inverse cosh-weighted Hilbert transform

This notebook is to study the computational characteristics of the Bertola-Katsevich-Tovbis inversion formulas for cosh-weighted Hilbert transform that solves the image reconstruciton problem in several practical imaging applications with considering the attenuation compenation.

There are four pairs of known analytical functions and their cosh-weighted Hilbert transforms. These four pairs are used to study the numerical behavior of the inversion formulas on error propagation and noise response. In particular, it is necessary to find out the range of $\mu$ for accurate and stable inversions with considering noise.

The implementation is based the Chebyshev series representation of finite Hilbert transform on the unit interval. The sampling points at the Chebyshev nodes for using the sine and cosine transforms.

Two main references are:

https://arxiv.org/abs/2002.12723

Bertola M, Katsevich A, and Tovbis A, "Inversion formulae for the cosh-weighted Hilbert transform", Proceeding of the American mathematical society, vol. 141, pp. 2703-2718, 2013.

## Setup Python libraries


```python
import numpy as np
import matplotlib.pyplot as plt
```

## Generate data


```python
sample_dim = 1000
step = np.pi/sample_dim

phi_step = np.arange(step/2.0, np.pi, step)
hilbert_mat = hilbert_cgl(sample_dim)
```


```python
mu = (1.0 + 1.0*1j)*np.pi
test_func = 3 # 0~4 for five test functions, otherwise goes to 0 for Chebyshev polynomials
poly_order = 32

f_0, F_0 = generate_test_func(phi_step, test_func, mu=mu, poly_order = poly_order)

sigma = 0.1
np.random.seed(11)
noise = np.random.normal(0.0, sigma, sample_dim)
F_1 = F_0 + noise*(1 + 1j)
mse_F = 0.0

if sigma > 0.0:
    der_F = der(F_0, F_1)
    print("source data-error-ratio: sigma = {0}, DER = {1:.2f}".format(sigma, der_F))
```

## Perform the numerical inversion and show result


```python
f_1 = inv_chfht(F_1, hilbert_mat, phi_step, mu)
show_result_complex(phi_step, f_0, f_1, F_0, F_1, mu)
```


```python
show_result_real(phi_step, f_0, f_1, F_0, F_1, mu/np.pi, True)
```

# Test functions, Hilbert transform matrix and visualization


```python
# five test functions
def test_sin(t, mu):
    f = np.sin(mu*np.sin(t))
    F = np.sinh(mu*np.cos(t))
    return f, F

def test_cos(t, mu):
    _cos = np.cos(t)
    _sin = np.sin(t)
    f = np.cos(mu*_sin)*_sin
    F = (_cos*np.cosh(mu*_cos)-0.5*mu*np.sinh(mu*_cos))
    return f, F

def test_ech1(t, mu):
    f = np.sin(t-mu*np.sin(t))
    _cos = np.cos(t)
    F = 0.5*np.exp(-mu*_cos)*(2*_cos+mu)
    return f, F
    
def test_ech2(t, mu):
    f = np.sin(2*t-mu*np.sin(t))
    _cos = np.cos(t)
    F = 0.5*np.exp(-mu*_cos)*(np.square(2*_cos)+2*mu*_cos+0.5*mu*mu-2.0)
    return f, F
    
def test_cheb(t, n):
    f = np.sin(n*t)
    F = np.cos(n*t)
    return f, F

# utility method to generate test functions
def generate_test_func(phi_step, test_func, mu=0.0, poly_order = 32):
    if test_func == 0:
        return test_cheb(phi_step, poly_order)
    elif test_func == 1:
        return test_cos(phi_step, mu)
    elif test_func == 2:
        return test_sin(phi_step, mu)
    elif test_func == 3:
        return test_ech1(phi_step, mu)
    elif test_func == 4:
        return test_ech2(phi_step, mu)
    else:
        return test_cheb(phi_step, poly_order)

# calculate data-error-ratio
def der(d0, d1):
    err = d1 - d0
    dat_norm = np.multiply(d0, np.conj(d0)) 
    dat_norm = np.mean(np.abs(dat_norm))
    err_norm = np.multiply(err, np.conj(err))
    err_norm = np.mean(np.abs(err_norm))
    der = np.log10(np.sqrt(dat_norm/err_norm))
    return der
    
# Chebyshev matrix for Hilbert transform
def hilbert_cgl(dim):
    grid = np.arange(0, 4*dim, 1)*(0.5*np.pi/dim)
    _cos = np.cos(grid)*np.sqrt(2.0/dim)
    _sin = np.sin(grid)*np.sqrt(2.0/dim)

    C3 = np.empty((dim, dim))
    S1 = np.empty((dim, dim))
    dim4 = 4*dim
    mm = 0
    while mm < dim:
        nn = 0
        while nn < dim:
            C3[mm][nn] = _cos[((2*mm+1)*nn)%dim4];
            S1[mm][nn] = _sin[((2*mm+1)*nn)%dim4];
            nn += 1
        mm += 1
    
    ht_cgl = np.matmul(S1, C3.transpose())
    _scale = np.pi/np.sin(np.arange(0.5, dim, 1)*(np.pi/dim))
    ht_cgl = np.multiply(_scale, ht_cgl.transpose()).transpose()
    
    return ht_cgl

# Bertola-Katsevich-Tovbis inversion formula by sine/cosine transform
def inv_chfht(F_1, ht_mat, phi_step, mu):
    sin_step = np.sin(phi_step)
    cos_w = np.cos(mu*sin_step)
    sin_w = np.sin(mu*sin_step)

    fht_cos = np.multiply(cos_w, F_1)
    fht_cos = np.matmul(ht_mat, fht_cos)
    fht_cos = np.multiply(cos_w, fht_cos)
    fht_cos = np.multiply(sin_step, fht_cos)

    fht_sin = np.multiply(sin_w, F_1)
    fht_sin = np.multiply(sin_step, fht_sin)
    fht_sin = np.matmul(ht_mat, fht_sin)
    fht_sin = np.multiply(sin_w, fht_sin)

    f_1 = (fht_cos + fht_sin)/np.pi
    return f_1

# visualize the inversion results
def show_result_complex(phi_step, f_0, f_1, F_0, F_1, mu):
    mse_F = F_1 - F_0
    mse_F = np.multiply(mse_F, np.conj(mse_F))
    mse_F = np.mean(np.abs(mse_F))   
    der_F = der(F_0, F_1)
    der_f = der(f_0, f_1)

    plt.figure(figsize=(12, 6))

    ax1 = plt.subplot(221)
    if mse_F > 0.0:
        plt.plot(phi_step, F_1.real, 'r', label='real ' + r'$\~F_{\mu}$')
        plt.title(r'$\mu={0:.1f}$'.format(mu) + ', DER = {0:.2f}'.format(der_F))
    else:
        plt.title(r'$\mu={0:.1f}$'.format(mu))

    plt.plot(phi_step, F_0.real, 'b', label='real ' + r'$F_{\mu}$')

    leg = plt.legend(loc='best', ncol=1)
    leg.get_frame().set_alpha(0.5)
    plt.setp(ax1.get_xticklabels(), visible=False)

    ax3 = plt.subplot(223)
    if mse_F > 0.0:
        plt.plot(phi_step, F_1.imag, 'r', label='imag ' + r'$\~F_{\mu}$')
    plt.plot(phi_step, F_0.imag, 'b', label='imag ' + r'$F_{\mu}$')

    leg = plt.legend(loc='best', ncol=1)
    leg.get_frame().set_alpha(0.5)

    ax2 = plt.subplot(222)
    if mse_F > 0.0:
        plt.plot(phi_step, f_1.real, 'r', label='real ' + r'$H_{\mu}^{-1}\~F_{\mu}$')
    else:
        plt.plot(phi_step, f_1.real, 'r', label='real ' + r'$H_{\mu}^{-1}F_{\mu}$')

    plt.plot(phi_step, f_0.real, 'b', label='real f')
    plt.title('inverse for ' + r'$\mu={0:.1f}$'.format(mu) + ', DER = {0:.2f}'.format(der_f))
    leg = plt.legend(loc='best', ncol=1)
    leg.get_frame().set_alpha(0.5)
    plt.setp(ax2.get_xticklabels(), visible=False)

    plt.subplot(224)
    if mse_F > 0.0:
        plt.plot(phi_step, f_1.imag, 'r', label='imag ' + r'$H_{\mu}^{-1}\~F_{\mu}$')
    else:
        plt.plot(phi_step, f_1.imag, 'r', label='imag ' + r'$H_{\mu}^{-1}F_{\mu}$')

    plt.plot(phi_step, f_0.imag, 'b', label='imag f')

    leg = plt.legend(loc='best', ncol=1)
    leg.get_frame().set_alpha(0.5)

    plt.show()
    
def show_result_real(phi_step, f_0, f_1, F_0, F_1, mu, show_real):
    mse_F = F_1 - F_0
    mse_F = np.multiply(mse_F, np.conj(mse_F))
    mse_F = np.mean(np.abs(mse_F))
    if show_real:
        f_0 = f_0.real
        f_1 = f_1.real
        F_0 = F_0.real
        F_1 = F_1.real
    else:
        f_0 = f_0.imag
        f_1 = f_1.imag
        F_0 = F_0.imag
        F_1 = F_1.imag
    der_F = der(F_0, F_1)
    der_f = der(f_0, f_1)

    plt.figure(figsize=(12, 3.0))

    plt.subplot(121)
    if mse_F > 0.0:
        plt.plot(phi_step, F_1, 'r', label=r'$\~F_{\mu}$')
    plt.plot(phi_step, F_0, 'b', label=r'$F_{\mu}$')
    if mse_F > 0.0:
        plt.title(r'$\mu={0}\pi$'.format(mu) + ': DER = {0:.2f}'.format(der_F))
    else:
        plt.title(r'$\mu={0}\pi$'.format(mu))
    leg = plt.legend(loc='best', ncol=2)
    leg.get_frame().set_alpha(0.5)

    plt.subplot(122)
    if mse_F > 0.0:
        plt.plot(phi_step, f_1, 'r', label=r'$H_{\mu}^{-1}\~F_{\mu}$')
    else:
        plt.plot(phi_step, f_1, 'r', label=r'$H_{\mu}^{-1}F_{\mu}$')
    plt.plot(phi_step, f_0, 'b', label='f')
    
    plt.title('inverse for ' + r'$\mu={0}\pi$'.format(mu) + ': DER = {0:.2f}'.format(der_f))
    leg = plt.legend(loc='best', ncol=2)
    leg.get_frame().set_alpha(0.5)

    plt.show()    
```


```python

```
