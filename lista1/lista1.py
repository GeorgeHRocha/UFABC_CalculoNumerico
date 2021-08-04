"""
Débora Vilela Bueno
Gabriela Rocha
George Harrison Rocha
Larissa Rodrigues de Almeida

Grupo 2 - Cálculo Numérico
Professor: André Pierro de Camargo
"""

from math import pi, sin


class P_H_Fixos:

    def __init__(self, p, h):
        self.p = p
        self.h = h

    def fi(self, x, q):
        return 2 + self.h ** 2 * q - ((1 - ((self.h ** 2) / 4) * self.p ** 2)) / x

    def sequencia(self, a0, n, q):
        seq = [a0]
        for _ in range(n - 1):
            a0 = self.fi(a0, q)
            seq.append(a0)

        return seq


# Classe usada no Exercício 1
class P_H_Q_Fixos(P_H_Fixos):

    def __init__(self, p, h, q):
        super().__init__(p, h)
        self.q = q

    def fi(self, x, q=None):
        return super().fi(x, self.q)

    def f(self, x):
        return self.fi(x) - x

    def menor_k(self, e):
        k = 0
        while (3/4) ** k >= e:
            k += 1

        return k

    def sequencia(self, a0, n, q=None):
        return super().sequencia(a0, n, self.q)


# Classe usada no Exercício 2
class P_H_N_Fixos(P_H_Fixos):

    def __init__(self, p, h, n):
        super().__init__(p, h)
        self.n = n

    def min_diag(self, q):
        return min(self.sequencia(2 + self.h ** 2 * q, self.n - 1, q))

    @staticmethod
    def sec_met(f, x_0, x_1, n_max, prec):
        iter = []
        i = 1
        x_k = x_0
        x_k_next = x_1

        while i < n_max and f(x_k_next - prec) * f(x_k_next + prec) > 0:
            temp = x_k - f(x_k) * (x_k_next - x_k) / (f(x_k_next) - f(x_k))
            x_k = x_k_next
            x_k_next = temp

            linha = [i, x_k, x_k_next, f(x_k_next), f(x_k_next - prec) * f(x_k_next + prec)]
            iter.append(linha)
            i += 1

        # Mostrando os resultados em uma tabela
        print('\nTODAS as iterações geradas pelo método das secantes')
        print(f'{"NumIter":8} {"iter":^23} {"new iter":^23} {"f(new iter)":^23} {"controle":^23}')
        for l in iter:
            print(f'{l[0]:7}', *map('{:23.16e}'.format, l[1:]))

        return iter


# Classe usada no Exercício 3
class DiferencasFinitas:

    def __init__(self, p, q, r, a, b, alpha, beta):
        self.p = p
        self.q = q
        self.r = r

        self.a = a
        self.b = b

        self.alpha = alpha
        self.beta = beta

    def ah(self, n):
        h = self.h(n)
        A = [[0] * (n - 1) for _ in range(n - 1)]

        for i in range(n - 1):
            A[i][i] = 2 + h ** 2 * self.q(self.x(i + 1, h))

        for i in range(n - 2):
            A[i + 1][i] = -1 - (h / 2) * self.p(self.x(i + 2, h))
            A[i][i + 1] = -1 + (h / 2) * self.p(self.x(i + 1, h))

        return A

    def x(self, i, h):
        return self.a + i * h

    @staticmethod
    def gauss_seidel(A, y, n_max, x0=None):
        m = len(y)
        if m == 1:
            return y
        if x0 is None:
            x0 = [0] * m
            for i in range(m):
                x0[i] = y[i] / A[i][i]

        iter = []
        for k in range(n_max):
            for i in range(m):
                x0[i] = y[i]
                for j in range(m):
                    if i != j:
                        x0[i] -= A[i][j] * x0[j]
                x0[i] /= A[i][i]

            iter.append(x0.copy())

        # Mostrando os resultados em uma tabela
        print('\n', ' ' * 4, *map('{:^18}'.format, [f'Iteracao {i + 1}' for i in range(n_max)]))
        for e, k in zip(zip(*iter), range(m)):
            print(f'[{k + 1},]', *map('{:.16f}'.format, e))

        return x0

    def aplicar_gauss_seidel(self, n, n_max, x0):
        return self.gauss_seidel(self.ah(n), self.vh(n), n_max, x0)

    @staticmethod
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

    def h(self, n):
        return (self.b - self.a) / n

    def vh(self, n):
        h = self.h(n)

        vh = [-self.r(self.x(1, h))
              + (1 / (h ** 2) + 1 / (2 * h) * self.p(self.x(1, h))) * self.alpha]

        for i in range(2, n - 1):
            vh.append(-self.r(self.x(i, h)))

        vh.append(-self.r(self.x(n - 1, h))
                  + (1 / (h ** 2) - 1 / (2 * h) * self.p(self.x(n - 1, h))) * self.beta)

        return list(map(lambda x: x * h ** 2, vh))

    def projh(self, n, y):
        h = self.h(n)
        projh = []

        for i in range(1, n):
            projh.append(y(self.x(i, h)))

        return projh

    def eh(self, n, y):
        projh = self.projh(n, y)
        yh = self.resolucao_sistema_gauss(self.ah(n), self.vh(n))

        return max(map(lambda x: abs(x[0] - x[1]), zip(projh, yh)))


####################


# EXERCÍCIO 1
p = 10
h = 1 / 10
q = 0

exercicio1 = P_H_Q_Fixos(p, h, q)

# 1 - Teórico

# 2 - Teórico

# 3 - a
print('f(1):', exercicio1.f(1))
print('f(2):', exercicio1.f(2))

# 3 - b
e = 10 ** -8
k = exercicio1.menor_k(e)
print('\nMenor valor de k que satisfaz a inequacao dada:', k)

# 3 - c
print('\nCalculando os primeiros k termos da sequência a0, a1, a2 ...:')
termos = exercicio1.sequencia(2, k)
for i, termo in enumerate(termos):
    print(f'{("α" + str(i)):>3}: {termo:.16f}')


# 3 - d - Teórico

# EXERCÍCIO 2

a = 0
b = 1
p = 0
n = 20
h = (b - a) / n

exercicio2 = P_H_N_Fixos(p, h, n)

# 1
print('\nminDiag(q) para os valores q = -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0:')
print(f'{"q":^4} {"minDiag(q)":^18}')
for q in range(-10, 1):
    print(f'{q:3}: {exercicio2.min_diag(q):18.15f}')

# 2 - Gráfico

# 3
prec = 10 ** -8
q0 = -10
q1 = -9.9
n_max = 1000

exercicio2.sec_met(exercicio2.min_diag, q0, q1, n_max, prec)


# EXERCÍCIO 3

p = lambda x: 0
q = lambda x: pi ** 2
r = lambda x: -2 * pi ** 2 * sin(pi * x)
a = 0
b = 1
alpha = beta = 0
equacao2 = DiferencasFinitas(p, q, r, a, b, alpha, beta)


# 1
n = 10
n_max = 7
equacao2.aplicar_gauss_seidel(n, n_max, [0] * (n - 1))

# 2
y = lambda x: sin(pi * x)

print(f'\n{"n":^5} {"En":^22}')
for n in range(100, 1001, 100):
    print(f'{n:4}: {equacao2.eh(n, y):.16e}')
