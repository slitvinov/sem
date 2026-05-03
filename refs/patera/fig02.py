import math


def cheb(N):
    x = [math.cos(math.pi * j / N) for j in range(N + 1)]
    c = [(2.0 if (j == 0 or j == N) else 1.0) * (-1) ** j for j in range(N + 1)]
    D = [[0.0] * (N + 1) for _ in range(N + 1)]
    for i in range(N + 1):
        for j in range(N + 1):
            if i != j:
                D[i][j] = (c[i] / c[j]) / (x[i] - x[j])
        D[i][i] = -sum(D[i][k] for k in range(N + 1) if k != i)
    return D, x


def matmul(A, B):
    n, p, m = len(A), len(B), len(B[0])
    return [[sum(A[i][k] * B[k][j] for k in range(p))
             for j in range(m)] for i in range(n)]


def solve_dense(A, b):
    n = len(b)
    A = [row[:] for row in A]
    b = b[:]
    for k in range(n):
        p = k
        for i in range(k + 1, n):
            if abs(A[i][k]) > abs(A[p][k]):
                p = i
        if p != k:
            A[k], A[p] = A[p], A[k]
            b[k], b[p] = b[p], b[k]
        for i in range(k + 1, n):
            f = A[i][k] / A[k][k]
            for j in range(k, n):
                A[i][j] -= f * A[k][j]
            b[i] -= f * b[k]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = b[i]
        for j in range(i + 1, n):
            s -= A[i][j] * x[j]
        x[i] = s / A[i][i]
    return x


def f(x):
    return math.cos(math.pi * x + math.pi / 4)


def exact(x):
    return -(math.cos(math.pi * x + math.pi / 4) + math.sqrt(2) / 2) / math.pi ** 2


def chebyshev_polys(N):
    T = [[0.0] * (N + 1) for _ in range(N + 1)]
    T[0][0] = 1.0
    if N >= 1:
        T[1][1] = 1.0
    for n in range(2, N + 1):
        for k in range(N + 1):
            a = 2.0 * T[n - 1][k - 1] if k > 0 else 0.0
            T[n][k] = a - T[n - 2][k]
    return T


def poly_deriv(p):
    return [k * p[k] for k in range(1, len(p))]


def poly_mult(p, q):
    r = [0.0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            r[i + j] += a * b
    return r


def poly_int_pm1(p):
    s = 0.0
    for k, c in enumerate(p):
        if k % 2 == 0:
            s += c * 2 / (k + 1)
    return s


def a_modal(N):
    T = chebyshev_polys(N)
    Tp = [poly_deriv(T[n]) for n in range(N + 1)]
    a = [[0.0] * (N + 1) for _ in range(N + 1)]
    for n in range(N + 1):
        for m in range(n, N + 1):
            v = poly_int_pm1(poly_mult(Tp[n], Tp[m])) if Tp[n] and Tp[m] else 0.0
            a[n][m] = v
            a[m][n] = v
    return a


def b_modal(N):
    b = [[0.0] * (N + 1) for _ in range(N + 1)]
    for n in range(N + 1):
        for m in range(N + 1):
            if (n + m) % 2 == 0:
                b[n][m] = 1.0 / (1 - (n + m) ** 2) + 1.0 / (1 - (n - m) ** 2)
    return b


def patera_AB(N, L):
    a = a_modal(N)
    b = b_modal(N)
    cb = [2.0 if (k == 0 or k == N) else 1.0 for k in range(N + 1)]
    T = [[math.cos(n * j * math.pi / N) for j in range(N + 1)] for n in range(N + 1)]
    A = [[0.0] * (N + 1) for _ in range(N + 1)]
    B = [[0.0] * (N + 1) for _ in range(N + 1)]
    pre = 4.0 / (N * N)
    for j in range(N + 1):
        for k in range(N + 1):
            sa = 0.0
            sb = 0.0
            for n in range(N + 1):
                for m in range(N + 1):
                    fac = T[n][j] * T[m][k] / (cb[n] * cb[m])
                    sa += fac * a[n][m]
                    sb += fac * b[n][m]
            A[j][k] = (2.0 / L) * pre * sa / (cb[j] * cb[k])
            B[j][k] = (L / 2.0) * pre * sb / (cb[j] * cb[k])
    return A, B


def spectral_1elem(N):
    A, B = patera_AB(N, 2.0)
    _, xb = cheb(N)
    fvals = [f(xb[j]) for j in range(N + 1)]
    Fg = [-sum(B[i][j] * fvals[j] for j in range(N + 1)) for i in range(N + 1)]
    K = [[A[i][j] for j in range(1, N)] for i in range(1, N)]
    rhs = [Fg[i] for i in range(1, N)]
    u_int = solve_dense(K, rhs)
    err = 0.0
    for i in range(N - 1):
        err = max(err, abs(u_int[i] - exact(xb[i + 1])))
    return err


def spectral_2elem(N):
    A, B = patera_AB(N, 1.0)
    _, xb = cheb(N)
    n_global = 2 * N + 1
    Ag = [[0.0] * n_global for _ in range(n_global)]
    Bg = [[0.0] * n_global for _ in range(n_global)]
    x_global = [0.0] * n_global
    for j in range(N + 1):
        x_global[N - j] = (xb[j] - 1) / 2
        x_global[2 * N - j] = (xb[j] + 1) / 2
    for j in range(N + 1):
        for k in range(N + 1):
            Ag[N - j][N - k] += A[j][k]
            Bg[N - j][N - k] += B[j][k]
    for j in range(N + 1):
        for k in range(N + 1):
            Ag[2 * N - j][2 * N - k] += A[j][k]
            Bg[2 * N - j][2 * N - k] += B[j][k]
    fvals = [f(x_global[i]) for i in range(n_global)]
    Fg = [-sum(Bg[i][j] * fvals[j] for j in range(n_global)) for i in range(n_global)]
    Asub = [[Ag[i][j] for j in range(1, n_global - 1)] for i in range(1, n_global - 1)]
    rhs = [Fg[i] for i in range(1, n_global - 1)]
    u_int = solve_dense(Asub, rhs)
    err = 0.0
    for i in range(n_global - 2):
        err = max(err, abs(u_int[i] - exact(x_global[i + 1])))
    return err


def quad_fem(M_elem):
    n = 2 * M_elem + 1
    h = 2.0 / M_elem
    K_ref = [[7/6, -4/3, 1/6],
             [-4/3, 8/3, -4/3],
             [1/6, -4/3, 7/6]]
    M_ref = [[4/15, 2/15, -1/15],
             [2/15, 16/15, 2/15],
             [-1/15, 2/15, 4/15]]
    x_global = [-1.0 + h * i / 2 for i in range(n)]
    f_global = [f(x) for x in x_global]
    Kg = [[0.0] * n for _ in range(n)]
    Mg = [[0.0] * n for _ in range(n)]
    for e in range(M_elem):
        i0 = 2 * e
        for i in range(3):
            for j in range(3):
                Kg[i0 + i][i0 + j] += (2 / h) * K_ref[i][j]
                Mg[i0 + i][i0 + j] += (h / 2) * M_ref[i][j]
    Fg = [-sum(Mg[i][j] * f_global[j] for j in range(n)) for i in range(n)]
    A = [[Kg[i][j] for j in range(1, n - 1)] for i in range(1, n - 1)]
    b = [Fg[i] for i in range(1, n - 1)]
    u_int = solve_dense(A, b)
    err = 0.0
    for i in range(n - 2):
        err = max(err, abs(u_int[i] - exact(x_global[i + 1])))
    return err


if __name__ == "__main__":
    print("# 1-element spectral")
    for N in [4, 5, 6, 7, 8, 10, 12, 15, 20]:
        Nt = N + 1
        e = spectral_1elem(N)
        print(Nt, math.log10(e) if e > 0 else -16)
    print("# 2-element spectral element")
    for N in [2, 3, 4, 5, 6, 7, 8, 10, 12]:
        Nt = 2 * N + 1
        e = spectral_2elem(N)
        print(Nt, math.log10(e) if e > 0 else -16)
    print("# quadratic FEM")
    for M_elem in [2, 3, 4, 5, 6, 8, 10, 12, 15, 20]:
        Nt = 2 * M_elem + 1
        e = quad_fem(M_elem)
        print(Nt, math.log10(e) if e > 0 else -16)
