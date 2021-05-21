import unittest
import numpy as np
from math import inf as infinity

from source.tbn import Tbn
from source.formulations.hilbert_basis import Formulation
from source.monomer import Monomer
from source.domain import Domain


class TestGraverBasis(unittest.TestCase):
    def setUp(self):
        pass

    def test_project_tbn_to_column_matrix(self):
        with self.subTest("slack version"):
            test_tbn = Tbn.from_string("a* b* \n a b \n a* \n b*")
            actual_matrix = Formulation._project_tbn_to_column_matrix(test_tbn)
            expected_matrix = np.array([
                [ 1, 1],  # a b
                [-1,-1],  # a* b*
                [-1, 0],  # a*
                [ 0,-1],  # b*
                [ 1, 0],  # a slack
                [ 0, 1],  # b slack
            ], np.int64).T
            self.assertTrue(np.array_equal(expected_matrix, actual_matrix))

        with self.subTest("no slack version"):
            test_tbn = Tbn.from_string("a* b* \n a b")  # no slack this time
            actual_matrix = Formulation._project_tbn_to_column_matrix(test_tbn)
            expected_matrix = np.array([
                [-1,-1],  # a b
                [ 1, 1],  # a* b*
            ], np.int64).T
            self.assertTrue(np.array_equal(expected_matrix, actual_matrix))

    def test_get_hilbert_basis_from_matrix(self):
        A = np.array([
            [1, -2, 1],
        ], np.int64)
        expected_hilbert_basis = np.array([
            [2, 1, 0],
            [1, 1, 1],
            [0, 1, 2],
        ], np.int64).T
        hilbert_basis = Formulation._get_hilbert_basis_from_matrix(A, quiet=True)
        self.assertEqual(expected_hilbert_basis.shape, hilbert_basis.shape)  # dimensions should be equal
        for expected_vector in expected_hilbert_basis.T:
            found = False
            for potential_match_vector in hilbert_basis.T:
                if np.array_equal(expected_vector, potential_match_vector):
                    found = True
            self.assertTrue(found)  # each element should have a match in the basis set
