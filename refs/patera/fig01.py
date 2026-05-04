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


def matvec(M, v):
    return [sum(M[i][j] * v[j] for j in range(len(v))) for i in range(len(M))]


def exact(x, t):
    return math.sin(math.pi * (t - x - 1)) if x < t - 1 else 0.0


def solve(N, nsteps):
    T = 5
    D, xb = cheb(N)
    dt = T / nsteps
    A = [[2.0 * dt * D[i][j] for j in range(N + 1)] for i in range(N + 1)]
    omega = math.pi * dt
    x1 = [(xb[j] - 1.0) / 2.0 for j in range(N + 1)]
    x2 = [(xb[j] + 1.0) / 2.0 for j in range(N + 1)]

    u1 = [0.0] * (N + 1)
    u2 = [0.0] * (N + 1)

    dx1_prev = None
    dx2_prev = None
    for step in range(1, nsteps + 1):
        dx1 = matvec(A, u1)
        dx2 = matvec(A, u2)
        dx2[N] = dx1[0]

        if dx1_prev is None:
            u1 = [u1[i] - dx1[i] for i in range(N + 1)]
            u2 = [u2[i] - dx2[i] for i in range(N + 1)]
        else:
            u1 = [u1[i] - 1.5 * dx1[i] + 0.5 * dx1_prev[i]
                  for i in range(N + 1)]
            u2 = [u2[i] - 1.5 * dx2[i] + 0.5 * dx2_prev[i]
                  for i in range(N + 1)]

        u1[N] = math.sin(omega * step)

        dx1_prev, dx2_prev = dx1, dx2

    err = 0.0
    for i in range(N + 1):
        err = max(err, abs(u1[i] - exact(x1[i], T)))
    for i in range(N):
        err = max(err, abs(u2[i] - exact(x2[i], T)))
    return err


ref = { }
with open("fig01.dat") as f:
    for line in f:
        n, x = line.split()
        ref[int(n)] = float(x)

for n in [7, 9, 11, 13, 15, 17, 19]:
    m = n
    ep = math.log10(solve((n - 1) // 2, m))
    while True:
        m *= 2
        e = math.log10(solve((n - 1) // 2, m))
        if abs(e - ep) < 0.01:
            break
        ep = e
    print("%2d%8.2f%8.2f%8d" % (n, e, ref[n], m))
