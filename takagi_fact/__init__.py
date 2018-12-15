from mpmath import mp

EPS = 1e-10


def approx_eq(a, b):
    return EPS > a - b > -EPS


def symmetric_svd(a):
    # Here we implicitly assume a is symmetric.
    A = mp.matrix(a)

    reA = A.apply(mp.re)
    imA = A.apply(mp.im)
    n = len(a)

    Bmat = [[0 for _ in range(2 * n)] for _ in range(2 * n)]
    for i in range(n):
        for j in range(n):
            Bmat[i][j] = reA[i, j]
            Bmat[i + n][j] = imA[i, j]
            Bmat[i][j + n] = imA[i, j]
            Bmat[i + n][j + n] = -reA[i, j]

    B = mp.matrix(Bmat)
    # Q.T * mp.diag(ev) * Q == B
    ev, Q = mp.eigsy(B)

    Qmat = [[Q[j, i] for j in range(2 * n)] for i in range(2 * n)]
    Umat = [[0 for _ in range(n)] for _ in range(n)]
    vs = [Qmat[i] for i, v in enumerate(ev) if v > 0]

    assert(len(vs) == n)

    for i in range(n):
        for j in range(n):
            Umat[i][j] = vs[j][i] - 1j * vs[j][i + n]

    Q = mp.matrix(Umat)

    # Accordingly, Q.T * A * Q is a diagonal matrix
    sing_mat = Q.T * A * Q
    sing_vs = [sing_mat[i, i] for i in range(len(sing_mat))]

    return sing_vs, Q


def takagi_fact(a):
    A = mp.matrix(a)
    sing_vs, Q = symmetric_svd(A)
    phase_mat = mp.diag([mp.exp(-1j * mp.arg(sing_v) / 2.0)
                         for sing_v in sing_vs])

    vs = [mp.fabs(sing_v) for sing_v in sing_vs]
    Qp = Q * phase_mat

    return vs, Qp


def set_dps(dps):
    mp.dps = dps
