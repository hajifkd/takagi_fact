import unittest
import takagi_fact
import itertools

import mpmath


class TakagiTest(unittest.TestCase):

    def test_symmtric_svd(self):
        N_REPEAT = 100
        N_DIM_UPPER = 2
        for _ in range(N_REPEAT):
            for n in range(2, N_DIM_UPPER + 1):
                ar = mpmath.randmatrix(n)
                ai = mpmath.randmatrix(n)
                a = ar + 1j * ai
                A = a + a.T

                vs, Q = takagi_fact.symmetric_svd(A)

                diag = mpmath.diag(vs)
                qtaq = Q.T * A * Q

                for i, j in itertools.product(range(n), repeat=2):
                    self.assertAlmostEqual(diag[i, j], qtaq[i, j])

    def test_takagi_fact(self):
        N_REPEAT = 100
        N_DIM_UPPER = 2
        for _ in range(N_REPEAT):
            for n in range(2, N_DIM_UPPER + 1):
                ar = mpmath.randmatrix(n)
                ai = mpmath.randmatrix(n)
                a = ar + 1j * ai
                A = a + a.T

                vs, Q = takagi_fact.takagi_fact(A)

                diag = mpmath.diag(vs)
                qtaq = Q.T * A * Q

                for i, j in itertools.product(range(n), repeat=2):
                    self.assertAlmostEqual(mpmath.im(diag[i, j]), 0)
                    self.assertAlmostEqual(diag[i, j], qtaq[i, j])
