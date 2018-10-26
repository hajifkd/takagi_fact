from mpmath import mp


def symmetric_svd(a):
    # Here we implicitly assume a is symmetric.
    A = mp.matrix(a)
    B = A.H * A

    # Q.H * mp.diag(ev) * Q == B
    # Note that we CANNOT use eigsy for complex symmetric matrix
    ev, Q = mp.eig(B)

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
