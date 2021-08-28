"""
Débora Vilela Bueno
Gabriela Rocha
George Harrison Rocha
Larissa Rodrigues de Almeida

Grupo 2 - Cálculo Numérico
Professor: André Pierro de Camargo
"""

from math import comb, pi, sin, cos, log
import numpy as np


c_m_q = [-1] * 31  # coeficientes minimos quadrados

# EXERCÍCIO 1

# Tarefa 1

# 1

k = 10
x = 0.7688


def P(j):
    """
    Função que retorna o j-ésimo polinômio de Legendre
    """
    return lambda x: sum(comb(j, i) * comb(j + i, i) * ((x - 1) / 2) ** i for i in range(0, j + 1))


for j in range(0, k):
    print(f'{P(j)(x):19.16f}')

# Tarefa 2

# 1
print()

m = 10 ** 4


def trapezios(a, b, m, f):
    t = np.linspace(a, b, m + 1)  # Cria um vetor com m + 1 números igualmente espaçados entre a e b.
    area = (b - a) / (2 * m) * (f(t[0]) + sum(2 * f(t[i]) for i in range(1, m)) + f(t[m]))
    return area


def r(x):
    return -pi ** 2 * (sin(pi * x) + cos(pi * x))


def alpha_j(j, m, regra_integracao):
    """
    Função que retorna o j_ésimo coeficiente mínimo quadrado.
    Se o j-ésimo coeficiente mínimo quadrado ainda não foi calculado, então calcula e guarda o resultado em c_m_q[j].
    """
    if c_m_q[j] == -1:
        p_r = lambda x: P(j - 1)(x) * r(x)
        c_m_q[j] = (2 * j - 1) / 2 * regra_integracao(-1, 1, m, p_r)
    return c_m_q[j]


k = 9
for j in range(1, k + 1):
    print(f'{alpha_j(j, m, trapezios):21.14e}')


# 2
print()


def Fj(j):
    """
    Função que calcula a integral dupla de fj, onde sua primitiva foi calculada usando conceitos de Cálculo.
    """
    j -= 1
    return lambda x: sum(comb(j, i) * comb(j + i, i) * ((x - 1) ** (i + 2) / ((i + 1) * (i + 2) * (2 ** i))) for i in range(0, j + 1))


def F(k, m, regra_integracao):
    """
    Função que calcula a integral dupla de f*
    """
    return lambda x: sum(alpha_j(j, m, regra_integracao) * Fj(j)(x) for j in range(1, k + 1))


alpha = beta = -1


def G(k, alpha, beta, m, regra_integracao):
    return lambda x: F(k, m, regra_integracao)(x) + (alpha - F(k, m, regra_integracao)(-1)) * (1 - x) / 2 + (beta - F(k, m, regra_integracao)(1)) * (x + 1) / 2


xs = [-1, -0.7, 0, 0.3, 1]
k = 7
Gk = G(k, alpha, beta, m, trapezios)
for x in xs:
    print(f'{Gk(x):17.14f}')


# 3
print()


def y(x):
    return sin(pi * x) + cos(pi * x)


def Emk(y, G, k, alpha, beta, m, regra_integracao):
    xs = np.linspace(-1, 1, 10118)  # Cria um vetor com 10118 números igualmente espaçados entre -1 e 1.
    return max(abs(y(x) - G(k, alpha, beta, m, regra_integracao)(x)) for x in xs)


print('  k', f'  {"Emk":^20}')
for k in range(2, 31):
    print(f'{k:3}: {Emk(y, G, k, alpha, beta, m, trapezios):21.14e}')

# TAREFA 3
m = 10 ** 5

# 1
print()
c_m_q = [-1] * 31  # "Zerando" os coeficientes minimos quadrados para recalculá-los com m = 10 ** 5

print('  k', f'  {"Emk":^20}')
for k in range(2, 31):
    print(f'{k:3}: {Emk(y, G, k, alpha, beta, m, trapezios):21.14e}')

# 2
print()
c_m_q = [-1] * 31  # "Zerando" os coeficientes minimos quadrados para recalculá-los com m = 10 ** 5 e usando simpson como regra de integração


def simpson(a, b, m, f):
    ts = np.linspace(a, b, m + 1)  # Cria um vetor com m + 1 números igualmente espaçados entre a e b.
    meio = sum(4 * f(ts[i]) if i % 2 != 0 else 2 * f(ts[i]) for i in range(1, m))
    area = ((b - a) / (3 * m)) * (f(ts[0]) + meio + f(ts[m]))
    return area


print('  k', f'  {"Emk":^20}')
for k in range(2, 31):
    print(f'{k:3}: {Emk(y, G, k, alpha, beta, m, simpson):21.14e}')

# 3 - GRÁFICO

# Tarefa 4


# 1
print()


def f_i(i):
    return lambda x: x ** (i - 1)


def fi_fj(i, j):
    fi = f_i(i)
    fj = f_i(j)

    return sum(fi(xk) * fj(xk) for xk in range(2, 17))


yks = list(map(log, [2.11503492462615e+00, 4.65926463996488e-01, 1.29124439761996e-01, 2.69084508442367e-02, 5.33968374774285e-03, 9.08993610154563e-04, 1.43101463317530e-04, 2.02304504324236e-05,2.64879763656189e-06,3.18669965637675e-07,3.57451284127563e-08,3.73634190253824e-09,3.71193520365409e-10,7.75705055744424e-11, 6.76469991134354e-11]))


def fi_y(i):
    fi = f_i(i)

    return sum(fi(xk) * yk for xk, yk in zip(range(2, 17), yks))


matriz = [[fi_fj(i, j) for j in range(1, 4)] for i in range(1, 4)]

print(*matriz, sep='\n')

vetor = [fi_y(i) for i in range(1, 4)]

print('\n',*vetor)


# 2
print()


def resolucao_sistema_gauss(A, y):
    n = len(y)

    # 1ª parte: Escalonamento (ou eliminação)

    for j in range(n - 1):
        ell = j

        while A[ell][j] == 0:
            ell += 1

            if ell > n:
                return None  # sistema singular

        for k in range(j, n):
            A[j][k], A[ell][k] = A[ell][k], A[j][k]

        y[j], y[ell] = y[ell], y[j]

        for i in range(j + 1, n):
            m = -A[i][j] / A[j][j]

            for k in range(j, n):
                A[i][k] += m*A[j][k]

            y[i] += m * y[j]

    # 2ª parte: substituição pra trás

    x = [0] * n
    for k in range(n - 1, -1, -1):
        x[k] = y[k]

        if k < n:
            for j in range(k + 1, n):
                x[k] -= A[k][j] * x[j]

        x[k] /= A[k][k]

    return x


print(resolucao_sistema_gauss(matriz, vetor))


# 3 - GRÁFICO


# EXERCÍCIO 2

ti = [0.000, 0.125, 0.250, 0.375, 0.500, 0.625, 0.750, 0.875, 1.000]
xi = [0.7416, 0.2685, 0.3333, 0.3982, -0.0749, -0.3089, 0.3333, 0.9756, 0.7416]
yi = [0.7416, 0.9756, 0.3333, -0.3089, -0.0749, 0.3982, 0.3333, 0.2685, 0.7416]
zi = [-0.4832, -0.2441, 0.3334, 0.9107, 1.1498, 0.9107, 0.3334, -0.2441, -0.4832]


def pol_interpolador(t, x, y):

    n = len(x) - 1
    M = [[0.0] * (n + 1) for _ in range(n + 1)]

    for i in range(n + 1):
        M[i][0] = y[i]

    for j in range(n):
        for i in range(n - j):
            M[i][j + 1] = (M[i + 1][j] - M[i][j]) / (x[i + j + 1] - x[i])
    #### Calculando o polinômio
    resp = 0
    fator = 1
    for i in range(n + 1):
        resp += M[0][i] * fator
        fator *= (t - x[i])

    return resp


px = []
py = []
pz = []

# 1
print()
print('t_j  ', *map('{:^18}'.format, ['x_j', 'y_j', 'z_j']))
for j in range(0, 51):
    tj = j / 50
    px.append(pol_interpolador(tj, ti, xi))
    py.append(pol_interpolador(tj, ti, yi))
    pz.append(pol_interpolador(tj, ti, zi))
    print(f'{tj:.2f} {px[-1]:18.15f} {py[-1] :18.15f} {pz[-1]:18.15f}')


# 2
print()
print(' j  tj  px(tj)+py(tj)+pz(tj)')
for j, e in enumerate(zip(px, py, pz)):
    tj = j / 50
    print(f'{j:2} {tj:.2f}  {sum(e):.16f}')
