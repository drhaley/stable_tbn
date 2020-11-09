import unittest
import numpy as np
from math import inf as infinity

from source.tbn import Tbn
from source.formulations.graver_basis import Formulation
from source.monomer import Monomer
from source.domain import Domain


class TestGraverBasis(unittest.TestCase):
    def setUp(self):
        pass

    def test_project_tbn_to_column_matrix(self):
        test_tbn = Tbn.from_string("a* b* \n a b \n a* \n b*")
        actual_matrix = Formulation._project_tbn_to_column_matrix(test_tbn)
        transposed_of_expected_matrix = np.array([
            [-1,-1],
            [ 1, 1],
            [-1, 0],
            [ 0,-1],
        ], np.int64)
        expected_matrix = np.transpose(transposed_of_expected_matrix)
        self.assertTrue(np.array_equal(expected_matrix, actual_matrix))

    def test_convert_matrix_to_hermite_normal_form(self):
        # test cases 1 and 2 taken from
        #  De Loera, Jesus A., Raymond Hemmecke, and Matthias Koeppe, eds.
        #  Algebraic and geometric ideas in the theory of discrete optimization.
        #  Society for Industrial and Applied Mathematics, 2012.
        with self.subTest("case 1"):
            test_matrix = np.array([
                [1, 1, 1, 1],
                [1, 5,10,25],
            ], np.int64)
            expected_H = np.array([
                [1, 0],
                [0, 1],
            ], np.int64)
            expected_C = np.array([
                [ 0, 1,-5, 5],
                [ 2,-2, 9,-6],
                [-1, 1,-4, 0],
                [ 0, 0, 0, 1],
            ], np.int64)
            actual_H, actual_C = Formulation._convert_matrix_to_hermite_normal_form(test_matrix)
            self.assertTrue(np.array_equal(expected_H, actual_H))
            self.assertTrue(np.array_equal(expected_C, actual_C))
        with self.subTest("case 2"):
            test_matrix = np.array([
                [2, 6, 1],
                [4, 7, 7],
            ], np.int64)
            expected_H = np.array([
                [1, 0],
                [2, 5],
            ], np.int64)
            expected_C = np.array([
                [ 4, 3, 7],
                [-1,-1,-2],
                [-1, 0,-2],
            ], np.int64)
            actual_H, actual_C = Formulation._convert_matrix_to_hermite_normal_form(test_matrix)
            self.assertTrue(np.array_equal(expected_H, actual_H))
            self.assertTrue(np.array_equal(expected_C, actual_C))

    def test_get_kernel_basis_from_hermite_normal_form(self):
        # test cases 1 and 2 taken from
        #  De Loera, Jesus A., Raymond Hemmecke, and Matthias Koeppe, eds.
        #  Algebraic and geometric ideas in the theory of discrete optimization.
        #  Society for Industrial and Applied Mathematics, 2012.
        with self.subTest("case 1"):
            H = np.array([
                [1, 0],
                [0, 1],
            ], np.int64)
            C = np.array([
                [ 0, 1,-5, 5],
                [ 2,-2, 9,-6],
                [-1, 1,-4, 0],
                [ 0, 0, 0, 1],
            ], np.int64)
            expected_kernel_basis = np.array([
                [-5, 5],
                [ 9,-6],
                [-4, 0],
                [ 0, 1],
            ], np.int64)
            actual_kernel_basis = Formulation._get_kernel_basis_from_hermite_normal_form(H, C)
            self.assertTrue(np.array_equal(expected_kernel_basis, actual_kernel_basis))
        with self.subTest("case 2"):
            H = np.array([
                [1, 0],
                [2, 5],
            ], np.int64)
            C = np.array([
                [ 4, 3, 7],
                [-1,-1,-2],
                [-1, 0,-2],
            ], np.int64)
            expected_kernel_basis = np.array([
                [7],
                [-2],
                [-2],
            ], np.int64)
            actual_kernel_basis = Formulation._get_kernel_basis_from_hermite_normal_form(H, C)
            self.assertTrue(np.array_equal(expected_kernel_basis, actual_kernel_basis))

    def test_get_graver_basis_from_kernel_basis(self):
        self.assertTrue(False)

    def test_integration(self):
        # this case from en.wikipedia.org/wiki/Graver_basis
        A = np.array([
            [1, 2, 1],
        ], np.int64)
        transpose_of_expected_graver_basis = np.array([
            [ 2,-1, 0],
            [ 0,-1, 2],
            [ 1, 0,-1],
            [ 1,-1, 1],
            [-2, 1, 0],
            [ 0, 1,-2],
            [-1, 0, 1],
            [-1, 1,-1],
        ], np.int64)
        expected_graver_basis = np.transpose(transpose_of_expected_graver_basis)
        H, C = Formulation._convert_matrix_to_hermite_normal_form(A)
        kernel_basis = Formulation._get_kernel_basis_from_hermite_normal_form(H, C)
        graver_basis = Formulation._get_graver_basis_from_kernel_basis(kernel_basis)
        for expected_vector in expected_graver_basis.T:
            found = False
            for potential_match_vector in graver_basis.T:
                if np.array_equal(expected_vector, potential_match_vector):
                    found = True
            self.assertTrue(found)
