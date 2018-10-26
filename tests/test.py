import unittest
import takagi_fact
import itertools

import mpmath


class TakagiTest(unittest.TestCase):

    def test_for_spec_a(self):
        n = 2
        A = mpmath.matrix([[0j, 1 + 0j], [1 + 0j, 0j]])
        vs, Q = takagi_fact.symmetric_svd(A)

        diag = mpmath.diag(vs)
        qtaq = Q.T * A * Q

        unit = mpmath.diag([1 for _ in range(n)])
        qdaq = Q.H * Q

        for i, j in itertools.product(range(n), repeat=2):
            self.assertAlmostEqual(diag[i, j], qtaq[i, j])
            self.assertAlmostEqual(unit[i, j], qdaq[i, j])

    def test_symmtric_svd(self):
        N_REPEAT = 2
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

                unit = mpmath.diag([1 for _ in range(n)])
                qdaq = Q.H * Q

                for i, j in itertools.product(range(n), repeat=2):
                    self.assertAlmostEqual(diag[i, j], qtaq[i, j])
                    self.assertAlmostEqual(unit[i, j], qdaq[i, j])

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

                unit = mpmath.diag([1 for _ in range(n)])
                qdaq = Q.H * Q

                for i, j in itertools.product(range(n), repeat=2):
                    self.assertAlmostEqual(mpmath.im(diag[i, j]), 0)
                    self.assertAlmostEqual(diag[i, j], qtaq[i, j])
                    self.assertAlmostEqual(unit[i, j], qdaq[i, j])
