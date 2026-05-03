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


def spectral_1elem(N):
    D, x = cheb(N)
    D2 = matmul(D, D)
    A = [[D2[i][j] for j in range(1, N)] for i in range(1, N)]
    b = [f(x[i]) for i in range(1, N)]
    u_int = solve_dense(A, b)
    err = 0.0
    for i in range(N - 1):
        err = max(err, abs(u_int[i] - exact(x[i + 1])))
    return err


def spectral_2elem(N):
    D, xi = cheb(N)
    D2 = matmul(D, D)
    n = 2 * N + 1
    A = [[0.0] * n for _ in range(n)]
    b = [0.0] * n

    A[0][0] = 1.0
    for j in range(1, N):
        r = N - j
        for jj in range(N + 1):
            A[r][N - jj] += 4 * D2[j][jj]
        b[r] = f((xi[j] - 1) / 2)

    for k in range(N + 1):
        A[N][N - k] += 2 * D[0][k]
        A[N][2 * N - k] -= 2 * D[N][k]

    for j in range(1, N):
        r = 2 * N - j
        for jj in range(N + 1):
            A[r][2 * N - jj] += 4 * D2[j][jj]
        b[r] = f((xi[j] + 1) / 2)

    A[2 * N][2 * N] = 1.0

    u = solve_dense(A, b)
    err = 0.0
    for j in range(N + 1):
        err = max(err, abs(u[N - j] - exact((xi[j] - 1) / 2)))
        err = max(err, abs(u[2 * N - j] - exact((xi[j] + 1) / 2)))
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
