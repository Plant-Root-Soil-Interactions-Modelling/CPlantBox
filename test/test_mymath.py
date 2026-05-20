"""
Tests for the math primitives in mymath.h (Vector2i, Vector2d, Vector3d, Matrix3d)
exposed via the plantbox Python binding.

Each test computes the expected result independently with numpy / pure Python
so the assertions are not circular.
"""

import math
import unittest

import numpy as np

import plantbox as pb


def v3(v):
    """Convert plantbox.Vector3d to a numpy array."""
    return np.array([v.x, v.y, v.z])


def m3(m):
    """Convert plantbox.Matrix3d to a 3x3 numpy array (row-major)."""
    return np.array(
        [
            [m.r0.x, m.r0.y, m.r0.z],
            [m.r1.x, m.r1.y, m.r1.z],
            [m.r2.x, m.r2.y, m.r2.z],
        ]
    )


EPS = 1e-12  # tolerance for floating-point comparisons


class TestVector2i(unittest.TestCase):

    def test_default_constructor(self):
        v = pb.Vector2i()
        self.assertEqual(v.x, 0)
        self.assertEqual(v.y, 0)

    def test_constructor(self):
        v = pb.Vector2i(3, -7)
        self.assertEqual(v.x, 3)
        self.assertEqual(v.y, -7)

    def test_copy_constructor(self):
        v = pb.Vector2i(5, 6)
        c = pb.Vector2i(v)
        self.assertEqual(c.x, 5)
        self.assertEqual(c.y, 6)

    def test_list_constructor(self):
        v = pb.Vector2i([2, 9])
        self.assertEqual(v.x, 2)
        self.assertEqual(v.y, 9)

    def test_readwrite_fields(self):
        v = pb.Vector2i(1, 2)
        v.x = 10
        v.y = 20
        self.assertEqual(v.x, 10)
        self.assertEqual(v.y, 20)

    def test_buffer_protocol(self):
        v = pb.Vector2i(3, 4)
        a = np.array(v, copy=False)
        self.assertEqual(a[0], 3)
        self.assertEqual(a[1], 4)

    def test_str(self):
        v = pb.Vector2i(1, 2)
        self.assertIn("1", str(v))
        self.assertIn("2", str(v))


class TestVector2d(unittest.TestCase):

    def test_default_constructor(self):
        v = pb.Vector2d()
        self.assertAlmostEqual(v.x, 0.0)
        self.assertAlmostEqual(v.y, 0.0)

    def test_constructor(self):
        v = pb.Vector2d(1.5, -2.5)
        self.assertAlmostEqual(v.x, 1.5)
        self.assertAlmostEqual(v.y, -2.5)

    def test_copy_constructor(self):
        v = pb.Vector2d(3.0, 4.0)
        c = pb.Vector2d(v)
        self.assertAlmostEqual(c.x, 3.0)
        self.assertAlmostEqual(c.y, 4.0)

    def test_list_constructor(self):
        v = pb.Vector2d([7.0, 8.0])
        self.assertAlmostEqual(v.x, 7.0)
        self.assertAlmostEqual(v.y, 8.0)

    def test_buffer_protocol(self):
        v = pb.Vector2d(1.0, 2.0)
        a = np.array(v, copy=False)
        self.assertAlmostEqual(a[0], 1.0)
        self.assertAlmostEqual(a[1], 2.0)

    def test_str(self):
        v = pb.Vector2d(1.0, 2.0)
        self.assertIn("1", str(v))
        self.assertIn("2", str(v))


class TestVector3d(unittest.TestCase):

    # ---- constructors --------------------------------------------------------

    def test_default_constructor(self):
        v = pb.Vector3d()
        self.assertEqual(v.x, 0.0)
        self.assertEqual(v.y, 0.0)
        self.assertEqual(v.z, 0.0)

    def test_xyz_constructor(self):
        v = pb.Vector3d(1.0, 2.0, 3.0)
        self.assertAlmostEqual(v.x, 1.0)
        self.assertAlmostEqual(v.y, 2.0)
        self.assertAlmostEqual(v.z, 3.0)

    def test_copy_constructor(self):
        v = pb.Vector3d(4.0, 5.0, 6.0)
        c = pb.Vector3d(v)
        self.assertAlmostEqual(c.x, v.x)
        self.assertAlmostEqual(c.y, v.y)
        self.assertAlmostEqual(c.z, v.z)

    def test_list_constructor(self):
        v = pb.Vector3d([7.0, 8.0, 9.0])
        self.assertAlmostEqual(v.x, 7.0)
        self.assertAlmostEqual(v.y, 8.0)
        self.assertAlmostEqual(v.z, 9.0)

    def test_readwrite_fields(self):
        v = pb.Vector3d(1.0, 2.0, 3.0)
        v.x = 10.0
        v.y = 20.0
        v.z = 30.0
        self.assertAlmostEqual(v.x, 10.0)
        self.assertAlmostEqual(v.y, 20.0)
        self.assertAlmostEqual(v.z, 30.0)

    # ---- buffer protocol (numpy zero-copy) -----------------------------------

    def test_buffer_protocol_values(self):
        v = pb.Vector3d(1.0, 2.0, 3.0)
        a = np.array(v, copy=False)
        np.testing.assert_array_almost_equal(a, [1.0, 2.0, 3.0])

    def test_buffer_protocol_zero_copy(self):
        """Modifying the numpy view should update the Vector3d in place."""
        v = pb.Vector3d(1.0, 2.0, 3.0)
        a = np.array(v, copy=False)
        a[0] = 99.0
        self.assertAlmostEqual(v.x, 99.0)

    # ---- length --------------------------------------------------------------

    def test_length_known(self):
        v = pb.Vector3d(3.0, 4.0, 0.0)
        self.assertAlmostEqual(v.length(), 5.0)

    def test_length_unit(self):
        v = pb.Vector3d(0.0, 0.0, 1.0)
        self.assertAlmostEqual(v.length(), 1.0)

    def test_length_general(self):
        x, y, z = 1.2, 3.4, 5.6
        v = pb.Vector3d(x, y, z)
        expected = math.sqrt(x * x + y * y + z * z)
        self.assertAlmostEqual(v.length(), expected)

    # ---- normalize -----------------------------------------------------------

    def test_normalize_in_place(self):
        v = pb.Vector3d(3.0, 4.0, 0.0)
        v.normalize()
        self.assertAlmostEqual(v.x, 0.6)
        self.assertAlmostEqual(v.y, 0.8)
        self.assertAlmostEqual(v.z, 0.0)
        self.assertAlmostEqual(v.length(), 1.0)

    def test_normalize_arbitrary(self):
        x, y, z = 2.0, -1.0, 3.0
        v = pb.Vector3d(x, y, z)
        v.normalize()
        self.assertAlmostEqual(v.length(), 1.0, places=12)

    # ---- dot product (times with Vector3d) -----------------------------------

    def test_dot_product_orthogonal(self):
        v1 = pb.Vector3d(1.0, 0.0, 0.0)
        v2 = pb.Vector3d(0.0, 1.0, 0.0)
        self.assertAlmostEqual(v1.times(v2), 0.0)

    def test_dot_product_parallel(self):
        v = pb.Vector3d(1.0, 2.0, 3.0)
        self.assertAlmostEqual(v.times(v), v.x**2 + v.y**2 + v.z**2)

    def test_dot_product_general(self):
        v1 = pb.Vector3d(1.0, 2.0, 3.0)
        v2 = pb.Vector3d(4.0, -5.0, 6.0)
        expected = 1 * 4 + 2 * (-5) + 3 * 6  # = 12
        self.assertAlmostEqual(v1.times(v2), expected)

    # ---- scalar multiply (times with double) ---------------------------------

    def test_scalar_multiply(self):
        v = pb.Vector3d(1.0, 2.0, 3.0)
        result = v.times(2.0)
        np.testing.assert_array_almost_equal(v3(result), [2.0, 4.0, 6.0])

    def test_scalar_multiply_zero(self):
        v = pb.Vector3d(1.0, 2.0, 3.0)
        result = v.times(0.0)
        np.testing.assert_array_almost_equal(v3(result), [0.0, 0.0, 0.0])

    def test_scalar_multiply_negative(self):
        v = pb.Vector3d(1.0, 2.0, 3.0)
        result = v.times(-1.0)
        np.testing.assert_array_almost_equal(v3(result), [-1.0, -2.0, -3.0])

    # ---- plus / minus --------------------------------------------------------

    def test_plus(self):
        v1 = pb.Vector3d(1.0, 2.0, 3.0)
        v2 = pb.Vector3d(4.0, 5.0, 6.0)
        result = v1.plus(v2)
        np.testing.assert_array_almost_equal(v3(result), [5.0, 7.0, 9.0])

    def test_minus(self):
        v1 = pb.Vector3d(4.0, 5.0, 6.0)
        v2 = pb.Vector3d(1.0, 2.0, 3.0)
        result = v1.minus(v2)
        np.testing.assert_array_almost_equal(v3(result), [3.0, 3.0, 3.0])

    def test_minus_self(self):
        v = pb.Vector3d(1.0, 2.0, 3.0)
        result = v.minus(v)
        np.testing.assert_array_almost_equal(v3(result), [0.0, 0.0, 0.0])

    # ---- cross product -------------------------------------------------------

    def test_cross_basis(self):
        """e_x × e_y = e_z"""
        ex = pb.Vector3d(1.0, 0.0, 0.0)
        ey = pb.Vector3d(0.0, 1.0, 0.0)
        result = ex.cross(ey)
        np.testing.assert_array_almost_equal(v3(result), [0.0, 0.0, 1.0])

    def test_cross_anticommutative(self):
        v1 = pb.Vector3d(1.0, 2.0, 3.0)
        v2 = pb.Vector3d(4.0, 5.0, 6.0)
        r1 = v3(v1.cross(v2))
        r2 = v3(v2.cross(v1))
        np.testing.assert_array_almost_equal(r1, -r2)

    def test_cross_orthogonal_to_inputs(self):
        v1 = pb.Vector3d(1.0, 2.0, 3.0)
        v2 = pb.Vector3d(4.0, 5.0, 6.0)
        c = v3(v1.cross(v2))
        self.assertAlmostEqual(float(np.dot(c, v3(v1))), 0.0, places=12)
        self.assertAlmostEqual(float(np.dot(c, v3(v2))), 0.0, places=12)

    def test_cross_known_values(self):
        v1 = pb.Vector3d(1.0, 2.0, 3.0)
        v2 = pb.Vector3d(4.0, 5.0, 6.0)
        # (2*6-3*5, 3*4-1*6, 1*5-2*4) = (-3, 6, -3)
        np.testing.assert_array_almost_equal(v3(v1.cross(v2)), [-3.0, 6.0, -3.0])

    # ---- rotAB (static) ------------------------------------------------------

    def test_rotAB_is_unit(self):
        """rotAB(a, b) should return a unit vector (first column of Rx(b)*Rz(a))."""
        for a, b in [(0.0, 0.0), (0.5, 0.3), (math.pi, math.pi / 4)]:
            result = pb.Vector3d.rotAB(a, b)
            self.assertAlmostEqual(result.length(), 1.0, places=12, msg=f"rotAB({a},{b}) not unit")

    def test_rotAB_zero_angles(self):
        """rotAB(0,0) should be (1,0,0)."""
        result = pb.Vector3d.rotAB(0.0, 0.0)
        np.testing.assert_array_almost_equal(v3(result), [1.0, 0.0, 0.0])

    def test_rotAB_matches_matrix_rotAB(self):
        """Vector3d.rotAB(a,b) == first column of Matrix3d.rotAB(a,b)."""
        a, b = 0.7, 0.4
        vec_result = v3(pb.Vector3d.rotAB(a, b))
        mat_result = m3(pb.Matrix3d.rotAB(a, b))[:, 0]
        np.testing.assert_array_almost_equal(vec_result, mat_result)

    # ---- str -----------------------------------------------------------------

    def test_str(self):
        v = pb.Vector3d(1.0, 2.0, 3.0)
        s = str(v)
        self.assertIn("1", s)
        self.assertIn("2", s)
        self.assertIn("3", s)


class TestMatrix3d(unittest.TestCase):

    # ---- constructors --------------------------------------------------------

    def test_default_identity(self):
        m = pb.Matrix3d()
        np.testing.assert_array_almost_equal(m3(m), np.eye(3))

    def test_nine_double_constructor(self):
        m = pb.Matrix3d(1, 2, 3, 4, 5, 6, 7, 8, 9)
        expected = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=float)
        np.testing.assert_array_almost_equal(m3(m), expected)

    def test_column_vector_constructor(self):
        """Matrix3d(c1, c2, c3) constructs from columns."""
        c0 = pb.Vector3d(1, 0, 0)
        c1 = pb.Vector3d(0, 1, 0)
        c2 = pb.Vector3d(0, 0, 1)
        m = pb.Matrix3d(c0, c1, c2)
        np.testing.assert_array_almost_equal(m3(m), np.eye(3))

    def test_copy_constructor(self):
        m = pb.Matrix3d(1, 2, 3, 4, 5, 6, 7, 8, 9)
        c = pb.Matrix3d(m)
        np.testing.assert_array_almost_equal(m3(c), m3(m))

    # ---- row / column accessors ----------------------------------------------

    def test_column(self):
        m = pb.Matrix3d(1, 2, 3, 4, 5, 6, 7, 8, 9)
        np.testing.assert_array_almost_equal(v3(m.column(0)), [1, 4, 7])
        np.testing.assert_array_almost_equal(v3(m.column(1)), [2, 5, 8])
        np.testing.assert_array_almost_equal(v3(m.column(2)), [3, 6, 9])

    def test_row(self):
        m = pb.Matrix3d(1, 2, 3, 4, 5, 6, 7, 8, 9)
        np.testing.assert_array_almost_equal(v3(m.row(0)), [1, 2, 3])
        np.testing.assert_array_almost_equal(v3(m.row(1)), [4, 5, 6])
        np.testing.assert_array_almost_equal(v3(m.row(2)), [7, 8, 9])

    # ---- rotation matrices ---------------------------------------------------

    def _check_rotation(self, M, a, axis):
        """Rotation matrix must be orthogonal and have det = +1."""
        Mn = m3(M)
        np.testing.assert_array_almost_equal(Mn @ Mn.T, np.eye(3), err_msg=f"rot{axis}({a}) not orthogonal")
        self.assertAlmostEqual(np.linalg.det(Mn), 1.0, places=12, msg=f"rot{axis}({a}) det != 1")

    def test_rotX_orthogonal(self):
        for a in [0, 0.3, math.pi / 2, math.pi]:
            self._check_rotation(pb.Matrix3d.rotX(a), a, "X")

    def test_rotY_orthogonal(self):
        for a in [0, 0.3, math.pi / 2, math.pi]:
            self._check_rotation(pb.Matrix3d.rotY(a), a, "Y")

    def test_rotZ_orthogonal(self):
        for a in [0, 0.3, math.pi / 2, math.pi]:
            self._check_rotation(pb.Matrix3d.rotZ(a), a, "Z")

    def test_rotX_known(self):
        """rotX(pi/2) should map e_y to e_z."""
        Rx = m3(pb.Matrix3d.rotX(math.pi / 2))
        result = Rx @ np.array([0.0, 1.0, 0.0])
        np.testing.assert_array_almost_equal(result, [0.0, 0.0, 1.0])

    def test_rotY_known(self):
        """rotY(pi/2) should map e_z to e_x."""
        Ry = m3(pb.Matrix3d.rotY(math.pi / 2))
        result = Ry @ np.array([0.0, 0.0, 1.0])
        np.testing.assert_array_almost_equal(result, [1.0, 0.0, 0.0])

    def test_rotZ_known(self):
        """rotZ(pi/2) should map e_x to e_y."""
        Rz = m3(pb.Matrix3d.rotZ(math.pi / 2))
        result = Rz @ np.array([1.0, 0.0, 0.0])
        np.testing.assert_array_almost_equal(result, [0.0, 1.0, 0.0])

    def test_rotX_zero(self):
        """rotX(0) == identity."""
        np.testing.assert_array_almost_equal(m3(pb.Matrix3d.rotX(0.0)), np.eye(3))

    def test_rotAB_is_Rx_b_times_Rz_a(self):
        """rotAB(a,b) == rotX(b) * rotZ(a) via numpy."""
        a, b = 0.7, 0.4
        expected = m3(pb.Matrix3d.rotX(b)) @ m3(pb.Matrix3d.rotZ(a))
        result = m3(pb.Matrix3d.rotAB(a, b))
        np.testing.assert_array_almost_equal(result, expected)

    def test_rotAB_orthogonal(self):
        Rab = m3(pb.Matrix3d.rotAB(0.5, 0.3))
        np.testing.assert_array_almost_equal(Rab @ Rab.T, np.eye(3))
        self.assertAlmostEqual(np.linalg.det(Rab), 1.0, places=12)

    # ---- det -----------------------------------------------------------------

    def test_det_identity(self):
        self.assertAlmostEqual(pb.Matrix3d().det(), 1.0)

    def test_det_known(self):
        # [[1,2,3],[4,5,6],[7,8,10]] -> det = -3
        m = pb.Matrix3d(1, 2, 3, 4, 5, 6, 7, 8, 10)
        expected = np.linalg.det(m3(m))
        self.assertAlmostEqual(m.det(), expected, places=10)

    def test_det_singular(self):
        # Rows 1 and 2 are identical -> det = 0
        m = pb.Matrix3d(1, 2, 3, 1, 2, 3, 7, 8, 9)
        self.assertAlmostEqual(m.det(), 0.0, places=10)

    # ---- inverse -------------------------------------------------------------

    def test_inverse_identity(self):
        inv = pb.Matrix3d().inverse()
        np.testing.assert_array_almost_equal(m3(inv), np.eye(3))

    def test_inverse_roundtrip(self):
        """M * M^{-1} == I"""
        m = pb.Matrix3d(1, 2, 3, 4, 5, 6, 7, 8, 10)
        Mn = m3(m)
        inv_np = np.linalg.inv(Mn)
        inv_pb = m3(m.inverse())
        np.testing.assert_array_almost_equal(inv_pb, inv_np, decimal=10)

    def test_inverse_product_is_identity(self):
        m = pb.Matrix3d(1, 2, 3, 0, 1, 4, 5, 6, 0)
        Mn = m3(m)
        inv = m3(m.inverse())
        np.testing.assert_array_almost_equal(Mn @ inv, np.eye(3), decimal=10)

    # ---- times (matrix * matrix) --------------------------------------------

    def test_times_matrix_identity(self):
        m = pb.Matrix3d(1, 2, 3, 4, 5, 6, 7, 8, 9)
        orig = m3(m).copy()
        m.times(pb.Matrix3d())  # multiply by identity (in-place)
        np.testing.assert_array_almost_equal(m3(m), orig)

    def test_times_matrix_general(self):
        A = pb.Matrix3d(1, 2, 3, 4, 5, 6, 7, 8, 10)
        B = pb.Matrix3d(9, 8, 7, 6, 5, 4, 3, 2, 1)
        An, Bn = m3(A).copy(), m3(B).copy()
        A.times(B)
        np.testing.assert_array_almost_equal(m3(A), An @ Bn)

    def test_times_matrix_non_commutative(self):
        # A = rotX(0.3), B = rotZ(0.7) — generic rotations never commute
        A = pb.Matrix3d.rotX(0.3)
        B = pb.Matrix3d.rotZ(0.7)
        An, Bn = m3(A).copy(), m3(B).copy()
        A.times(B)
        result_AB = m3(A).copy()
        expected_AB = An @ Bn
        np.testing.assert_array_almost_equal(result_AB, expected_AB)
        # verify AB != BA
        result_BA = Bn @ An
        self.assertFalse(np.allclose(result_AB, result_BA), "rotX and rotZ should not commute for these angles")

    # ---- times (matrix * vector) --------------------------------------------

    def test_times_vector_identity(self):
        v = pb.Vector3d(1.0, 2.0, 3.0)
        result = pb.Matrix3d().times(v)
        np.testing.assert_array_almost_equal(v3(result), v3(v))

    def test_times_vector_general(self):
        M = pb.Matrix3d(1, 2, 3, 4, 5, 6, 7, 8, 9)
        v = pb.Vector3d(1.0, 0.0, 0.0)
        result = M.times(v)
        expected = m3(M) @ np.array([1.0, 0.0, 0.0])
        np.testing.assert_array_almost_equal(v3(result), expected)

    def test_times_vector_rotation(self):
        """Rotation by 90° around Z should map (1,0,0) to (0,1,0)."""
        Rz = pb.Matrix3d.rotZ(math.pi / 2)
        ex = pb.Vector3d(1.0, 0.0, 0.0)
        result = Rz.times(ex)
        np.testing.assert_array_almost_equal(v3(result), [0.0, 1.0, 0.0])

    # ---- ons (orthonormal system) --------------------------------------------

    def test_ons_first_column_is_normalized_input(self):
        """After ons(v), the first column should be the normalized v."""
        vec = pb.Vector3d(1.0, 2.0, 3.0)
        v_norm = v3(vec) / np.linalg.norm(v3(vec))
        M = pb.Matrix3d.ons(vec)
        np.testing.assert_array_almost_equal(v3(M.column(0)), v_norm)

    def test_ons_is_orthonormal(self):
        """Columns of ons result should be mutually orthogonal unit vectors."""
        for xyz in [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 2, 3), (-1, 0.5, 2), (0.1, -0.2, 0.9)]:
            vec = pb.Vector3d(*xyz)
            M = pb.Matrix3d.ons(vec)
            Mn = m3(M)
            # columns should be unit vectors
            for i in range(3):
                self.assertAlmostEqual(np.linalg.norm(Mn[:, i]), 1.0, places=12, msg=f"column {i} not unit for input {xyz}")
            # columns should be mutually orthogonal
            self.assertAlmostEqual(float(np.dot(Mn[:, 0], Mn[:, 1])), 0.0, places=12, msg=f"col0·col1 != 0 for {xyz}")
            self.assertAlmostEqual(float(np.dot(Mn[:, 0], Mn[:, 2])), 0.0, places=12, msg=f"col0·col2 != 0 for {xyz}")
            self.assertAlmostEqual(float(np.dot(Mn[:, 1], Mn[:, 2])), 0.0, places=12, msg=f"col1·col2 != 0 for {xyz}")

    def test_ons_determinant_is_plus_one(self):
        """ons should produce a proper rotation (not a reflection)."""
        for xyz in [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 2, 3)]:
            vec = pb.Vector3d(*xyz)
            M = pb.Matrix3d.ons(vec)
            det = np.linalg.det(m3(M))
            self.assertAlmostEqual(det, 1.0, places=12, msg=f"ons det != 1 for {xyz}")

    # ---- buffer protocol -----------------------------------------------------

    def test_buffer_protocol_values(self):
        m = pb.Matrix3d(1, 2, 3, 4, 5, 6, 7, 8, 9)
        a = np.array(m, copy=False)
        expected = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=float)
        np.testing.assert_array_almost_equal(a, expected)

    # ---- str -----------------------------------------------------------------

    def test_str(self):
        m = pb.Matrix3d()
        s = str(m)
        self.assertIn("1", s)
        self.assertIn("0", s)


class TestTurtle3D(unittest.TestCase):

    # ---- construction --------------------------------------------------------

    def test_default_construction(self):
        t = pb.Turtle3D()
        np.testing.assert_array_almost_equal(v3(t.getPosition()), [0.0, 0.0, 0.0])
        np.testing.assert_array_almost_equal(v3(t.heading()), [1.0, 0.0, 0.0])
        np.testing.assert_array_almost_equal(v3(t.left()), [0.0, 1.0, 0.0])
        np.testing.assert_array_almost_equal(v3(t.up()), [0.0, 0.0, 1.0])

    def test_custom_construction(self):
        pos = pb.Vector3d(1.0, 2.0, 3.0)
        frame = pb.Matrix3d.rotZ(math.pi / 4)
        t = pb.Turtle3D(pos, frame)
        np.testing.assert_array_almost_equal(v3(t.getPosition()), [1.0, 2.0, 3.0])
        np.testing.assert_array_almost_equal(m3(t.getFrame()), m3(frame))

    # ---- forward -------------------------------------------------------------

    def test_forward_default_heading(self):
        """Default heading is +x; forward(d) moves by d along +x."""
        t = pb.Turtle3D()
        t.forward(5.0)
        np.testing.assert_array_almost_equal(v3(t.getPosition()), [5.0, 0.0, 0.0])

    def test_forward_along_rotated_heading(self):
        """After turnLeft(pi/2) heading is +y; forward(d) moves along +y."""
        t = pb.Turtle3D()
        t.turnLeft(math.pi / 2)
        t.forward(3.0)
        np.testing.assert_array_almost_equal(v3(t.getPosition()), [0.0, 3.0, 0.0], decimal=12)

    def test_forward_does_not_change_frame(self):
        t = pb.Turtle3D()
        orig = m3(t.getFrame()).copy()
        t.forward(10.0)
        np.testing.assert_array_almost_equal(m3(t.getFrame()), orig)

    # ---- yaw (turnLeft / turnRight) ------------------------------------------

    def test_turnLeft_90(self):
        """turnLeft(pi/2) rotates heading from +x to +y."""
        t = pb.Turtle3D()
        t.turnLeft(math.pi / 2)
        np.testing.assert_array_almost_equal(v3(t.heading()), [0.0, 1.0, 0.0], decimal=12)

    def test_turnRight_90(self):
        """turnRight(pi/2) rotates heading from +x to -y."""
        t = pb.Turtle3D()
        t.turnRight(math.pi / 2)
        np.testing.assert_array_almost_equal(v3(t.heading()), [0.0, -1.0, 0.0], decimal=12)

    def test_turnLeft_turnRight_inverse(self):
        t = pb.Turtle3D()
        orig = m3(t.getFrame()).copy()
        t.turnLeft(0.7)
        t.turnRight(0.7)
        np.testing.assert_array_almost_equal(m3(t.getFrame()), orig, decimal=12)

    # ---- pitch (pitchUp / pitchDown) -----------------------------------------

    def test_pitchUp_90(self):
        """pitchUp(pi/2): heading rotates around L=(0,1,0); result = (0,0,-1)."""
        t = pb.Turtle3D()
        t.pitchUp(math.pi / 2)
        np.testing.assert_array_almost_equal(v3(t.heading()), [0.0, 0.0, -1.0], decimal=12)

    def test_pitchDown_90(self):
        """pitchDown(pi/2): heading rotates to (0,0,1)."""
        t = pb.Turtle3D()
        t.pitchDown(math.pi / 2)
        np.testing.assert_array_almost_equal(v3(t.heading()), [0.0, 0.0, 1.0], decimal=12)

    def test_pitchUp_pitchDown_inverse(self):
        t = pb.Turtle3D()
        orig = m3(t.getFrame()).copy()
        t.pitchUp(0.5)
        t.pitchDown(0.5)
        np.testing.assert_array_almost_equal(m3(t.getFrame()), orig, decimal=12)

    # ---- roll (rollLeft / rollRight) -----------------------------------------

    def test_rollLeft_90(self):
        """rollLeft(pi/2): left moves from +y to +z."""
        t = pb.Turtle3D()
        t.rollLeft(math.pi / 2)
        np.testing.assert_array_almost_equal(v3(t.left()), [0.0, 0.0, 1.0], decimal=12)

    def test_rollRight_90(self):
        """rollRight(pi/2): left moves from +y to -z."""
        t = pb.Turtle3D()
        t.rollRight(math.pi / 2)
        np.testing.assert_array_almost_equal(v3(t.left()), [0.0, 0.0, -1.0], decimal=12)

    def test_rollLeft_rollRight_inverse(self):
        t = pb.Turtle3D()
        orig = m3(t.getFrame()).copy()
        t.rollLeft(1.1)
        t.rollRight(1.1)
        np.testing.assert_array_almost_equal(m3(t.getFrame()), orig, decimal=12)

    # ---- frame orthonormality ------------------------------------------------

    def test_frame_stays_orthonormal_after_many_rotations(self):
        t = pb.Turtle3D()
        for _ in range(10):
            t.turnLeft(0.3)
            t.pitchUp(0.7)
            t.rollLeft(1.1)
            t.turnRight(0.4)
            t.pitchDown(0.2)
        cols = [v3(t.heading()), v3(t.left()), v3(t.up())]
        for i, c in enumerate(cols):
            self.assertAlmostEqual(np.linalg.norm(c), 1.0, places=10, msg=f"column {i} not unit after rotations")
        self.assertAlmostEqual(float(np.dot(cols[0], cols[1])), 0.0, places=10)
        self.assertAlmostEqual(float(np.dot(cols[0], cols[2])), 0.0, places=10)
        self.assertAlmostEqual(float(np.dot(cols[1], cols[2])), 0.0, places=10)

    # ---- setters -------------------------------------------------------------

    def test_setPosition(self):
        t = pb.Turtle3D()
        t.setPosition(pb.Vector3d(4.0, 5.0, 6.0))
        np.testing.assert_array_almost_equal(v3(t.getPosition()), [4.0, 5.0, 6.0])

    def test_setFrame(self):
        t = pb.Turtle3D()
        f = pb.Matrix3d.rotX(0.5)
        t.setFrame(f)
        np.testing.assert_array_almost_equal(m3(t.getFrame()), m3(f))

    # ---- str -----------------------------------------------------------------

    def test_str(self):
        t = pb.Turtle3D()
        s = str(t)
        self.assertIn("pos", s)


class TestTurtlePolyline(unittest.TestCase):
    """Tests for TurtlePolyline (formerly Meristem) exposed as pb.TurtlePolyline."""

    # ---- construction --------------------------------------------------------

    def test_default_construction(self):
        """Default TurtlePolyline starts with 0 nodes; anchor is at the origin."""
        m = pb.TurtlePolyline()
        self.assertEqual(m.size(), 0)
        np.testing.assert_array_almost_equal(v3(m.getAnchor()), [0.0, 0.0, 0.0])

    def test_default_initial_node_index(self):
        """initialNodeIdx starts at 0 (deque index of the initial node)."""
        m = pb.TurtlePolyline()
        self.assertEqual(m.getInitialNodeIndex(), 0)

    def test_custom_construction(self):
        pos = pb.Vector3d(1.0, 2.0, 3.0)
        m = pb.TurtlePolyline(pos, pb.Matrix3d())
        np.testing.assert_array_almost_equal(v3(m.getAnchor()), [1.0, 2.0, 3.0])
        self.assertEqual(m.size(), 0)

    def test_custom_construction_first_node_at_anchor(self):
        """getNode(0) after addNodeBack(0) returns the anchor position."""
        pos = pb.Vector3d(1.0, 2.0, 3.0)
        m = pb.TurtlePolyline(pos, pb.Matrix3d())
        m.addNodeBack(0.0)
        np.testing.assert_array_almost_equal(v3(m.getNode(0)), [1.0, 2.0, 3.0])

    # ---- addNodeBack ---------------------------------------------------------

    def test_addNodeBack_increases_size(self):
        m = pb.TurtlePolyline()
        m.addNodeBack(1.0)
        self.assertEqual(m.size(), 1)
        m.addNodeBack(1.0)
        self.assertEqual(m.size(), 2)

    def test_addNodeBack_does_not_change_initial_node_index(self):
        m = pb.TurtlePolyline()
        m.addNodeBack(1.0)
        m.addNodeBack(1.0)
        self.assertEqual(m.getInitialNodeIndex(), 0)

    def test_getNode0_straight(self):
        """Node 0 = anchor + dist*heading."""
        m = pb.TurtlePolyline()
        m.addNodeBack(5.0)
        np.testing.assert_array_almost_equal(v3(m.getNode(0)), [5.0, 0.0, 0.0], decimal=12)

    def test_getNode_with_yaw(self):
        """Node after yaw=pi/2 then forward(3) should land on +y axis."""
        m = pb.TurtlePolyline()
        m.addNodeBack(3.0, yaw=math.pi / 2)
        np.testing.assert_array_almost_equal(v3(m.getNode(0)), [0.0, 3.0, 0.0], decimal=12)

    def test_getNode_chained(self):
        """Two straight steps of length 2 give cumulative positions."""
        m = pb.TurtlePolyline()
        m.addNodeBack(2.0)
        m.addNodeBack(2.0)
        np.testing.assert_array_almost_equal(v3(m.getNode(0)), [2.0, 0.0, 0.0], decimal=12)
        np.testing.assert_array_almost_equal(v3(m.getNode(1)), [4.0, 0.0, 0.0], decimal=12)

    # ---- addNodeFront --------------------------------------------------------

    def test_addNodeFront_increases_size(self):
        m = pb.TurtlePolyline()
        m.addNodeFront(1.0)
        self.assertEqual(m.size(), 1)

    def test_addNodeFront_increments_initial_node_index(self):
        """Each addNodeFront() prepends a node, shifting initialNodeIdx by 1."""
        m = pb.TurtlePolyline()
        self.assertEqual(m.getInitialNodeIndex(), 0)
        m.addNodeFront(1.0)
        self.assertEqual(m.getInitialNodeIndex(), 1)
        m.addNodeFront(1.0)
        self.assertEqual(m.getInitialNodeIndex(), 2)

    def test_addNodeFront_shifts_indices(self):
        """addNodeFront prepends; old node 0 becomes node 1."""
        m = pb.TurtlePolyline()
        m.addNodeBack(5.0)  # deque: [{0,0,0,5}]; node 0=(5,0,0)
        m.addNodeFront(3.0)  # deque: [{0,0,0,3}, {0,0,0,5}]; initialNodeIdx=1
        np.testing.assert_array_almost_equal(v3(m.getNode(0)), [3.0, 0.0, 0.0], decimal=12)
        np.testing.assert_array_almost_equal(v3(m.getNode(1)), [8.0, 0.0, 0.0], decimal=12)

    def test_initial_node_index_after_addNodeFront(self):
        """initialNodeIdx is incremented by each addNodeFront call."""
        m = pb.TurtlePolyline()
        m.addNodeBack(5.0)  # initialNodeIdx stays 0
        m.addNodeFront(3.0)  # prepended; initialNodeIdx becomes 1
        self.assertEqual(m.getInitialNodeIndex(), 1)

    # ---- getLength -----------------------------------------------------------

    def test_getLength_empty(self):
        """Default polyline has no nodes; total length == 0."""
        m = pb.TurtlePolyline()
        self.assertAlmostEqual(m.getLength(), 0.0)

    def test_getLength_total(self):
        """getLength() sums all segment distances."""
        m = pb.TurtlePolyline()
        m.addNodeBack(3.0)
        m.addNodeBack(4.0)
        m.addNodeBack(5.0)
        self.assertAlmostEqual(m.getLength(), 12.0)

    def test_getLength_n_zero(self):
        """getLength(0) == 0."""
        m = pb.TurtlePolyline()
        m.addNodeBack(3.0)
        self.assertAlmostEqual(m.getLength(0), 0.0)

    def test_getLength_n_partial(self):
        """getLength(n) returns cumulative length of first n nodes."""
        m = pb.TurtlePolyline()
        # deque: [{3}, {4}, {5}]
        m.addNodeBack(3.0)
        m.addNodeBack(4.0)
        m.addNodeBack(5.0)
        self.assertAlmostEqual(m.getLength(1), 3.0)  # first segment
        self.assertAlmostEqual(m.getLength(2), 7.0)  # first two segments
        self.assertAlmostEqual(m.getLength(3), 12.0)  # all nodes

    def test_getLength_n_equals_total(self):
        """getLength(size()) == getLength()."""
        m = pb.TurtlePolyline()
        m.addNodeBack(2.0)
        m.addNodeBack(3.0)
        self.assertAlmostEqual(m.getLength(m.size()), m.getLength())

    def test_getLength_addNodeFront_counts(self):
        """getLength() includes segments added at the front."""
        m = pb.TurtlePolyline()
        m.addNodeBack(5.0)
        m.addNodeFront(3.0)
        self.assertAlmostEqual(m.getLength(), 8.0)

    # ---- getNodeIndexAtLength ------------------------------------------------

    def test_getNodeIndexAtLength_zero(self):
        """getNodeIndexAtLength(0) returns 0 (clamped)."""
        m = pb.TurtlePolyline()
        m.addNodeBack(3.0)
        m.addNodeBack(4.0)
        self.assertEqual(m.getNodeIndexAtLength(0.0), 0)

    def test_getNodeIndexAtLength_exact_boundary(self):
        """getNodeIndexAtLength at exact cumulative distance boundary."""
        m = pb.TurtlePolyline()
        # deque: [{3}, {4}]  cumulative: [3, 7]
        m.addNodeBack(3.0)
        m.addNodeBack(4.0)
        # at s=3 the running sum first reaches 3 at index 0
        # roundUp=False → return max(0, 0-1)=0; roundUp=True → return 0
        self.assertEqual(m.getNodeIndexAtLength(3.0, False), 0)
        self.assertEqual(m.getNodeIndexAtLength(3.0, True), 0)

    def test_getNodeIndexAtLength_beyond_total(self):
        """getNodeIndexAtLength(s) for s > total returns last index."""
        m = pb.TurtlePolyline()
        m.addNodeBack(3.0)
        m.addNodeBack(4.0)
        self.assertEqual(m.getNodeIndexAtLength(100.0), m.size() - 1)

    # ---- getNodeFrame --------------------------------------------------------

    def test_getNodeFrame_identity_anchor(self):
        """With identity anchor frame, frame at any straight node is also identity."""
        m = pb.TurtlePolyline()
        m.addNodeBack(5.0)
        m.addNodeBack(5.0)
        for i in range(m.size()):
            np.testing.assert_array_almost_equal(m3(m.getNodeFrame(i)), np.eye(3), decimal=12)

    def test_getNodeFrame_after_yaw(self):
        """After a yaw=pi/2 node the frame heading should be +y."""
        m = pb.TurtlePolyline()
        m.addNodeBack(1.0, yaw=math.pi / 2)
        frame = m.getNodeFrame(0)
        np.testing.assert_array_almost_equal(v3(frame.column(0)), [0.0, 1.0, 0.0], decimal=12)

    def test_getNodeFrame_stays_orthonormal(self):
        """Frame at every node must be an orthonormal matrix."""
        m = pb.TurtlePolyline()
        m.addNodeBack(2.0, yaw=0.3)
        m.addNodeBack(2.0, pitch=0.5)
        m.addNodeBack(2.0, roll=0.7)
        for i in range(m.size()):
            F = m3(m.getNodeFrame(i))
            np.testing.assert_array_almost_equal(F @ F.T, np.eye(3), decimal=10, err_msg=f"frame at node {i} not orthonormal")
            self.assertAlmostEqual(np.linalg.det(F), 1.0, places=10, msg=f"frame det != 1 at node {i}")

    def test_getNodeFrame_anchor_frame_respected(self):
        """Anchor frame rotZ(pi/2) means heading starts as +y; straight nodes keep +y."""
        m = pb.TurtlePolyline()
        m.setAnchorFrame(pb.Matrix3d.rotZ(math.pi / 2))
        m.addNodeBack(3.0)
        frame = m.getNodeFrame(0)
        np.testing.assert_array_almost_equal(v3(frame.column(0)), [0.0, 1.0, 0.0], decimal=12)

    # ---- getPolyline ---------------------------------------------------------

    def test_getPolyline_length(self):
        """getPolyline() has one entry per deque node (3 added = 3)."""
        m = pb.TurtlePolyline()
        m.addNodeBack(1.0)
        m.addNodeBack(1.0)
        m.addNodeBack(1.0)
        self.assertEqual(len(m.getPolyline()), 3)

    def test_getPolyline_straight(self):
        """Steps of dist=2 give positions [2, 4, 6] along x."""
        m = pb.TurtlePolyline()
        for _ in range(3):
            m.addNodeBack(2.0)
        pts = m.getPolyline()
        for i, p in enumerate(pts):
            np.testing.assert_array_almost_equal(v3(p), [2.0 * (i + 1), 0.0, 0.0], decimal=12)

    # ---- getNodes and getTurtleNode (raw turtle data) -----------------------

    def test_getNodes_empty_on_default_construction(self):
        """getNodes() returns an empty list for a default-constructed TurtlePolyline."""
        m = pb.TurtlePolyline()
        nodes = m.getNodes()
        self.assertEqual(len(nodes), 0)

    def test_getNodes_returns_raw_nodes(self):
        """getNodes() returns all deque entries."""
        m = pb.TurtlePolyline()
        m.addNodeBack(1.0, yaw=0.5, pitch=0.3, roll=0.1)
        nodes = m.getNodes()
        self.assertEqual(len(nodes), 1)
        self.assertAlmostEqual(nodes[0].yaw, 0.5)
        self.assertAlmostEqual(nodes[0].pitch, 0.3)
        self.assertAlmostEqual(nodes[0].roll, 0.1)
        self.assertAlmostEqual(nodes[0].dist, 1.0)

    def test_getTurtleNode_added(self):
        """getTurtleNode(i) returns the TurtleNode at deque index i."""
        m = pb.TurtlePolyline()
        m.addNodeBack(2.5, yaw=0.1, pitch=0.2, roll=0.3)
        n = m.getTurtleNode(0)
        self.assertAlmostEqual(n.dist, 2.5)
        self.assertAlmostEqual(n.yaw, 0.1)
        self.assertAlmostEqual(n.pitch, 0.2)
        self.assertAlmostEqual(n.roll, 0.3)

    def test_getTurtleNode_matches_getNodes(self):
        """getTurtleNode(i) and getNodes()[i] refer to the same data."""
        m = pb.TurtlePolyline()
        m.addNodeBack(3.0, yaw=0.7)
        nodes = m.getNodes()
        for i in range(m.size()):
            n = m.getTurtleNode(i)
            self.assertAlmostEqual(n.dist, nodes[i].dist)
            self.assertAlmostEqual(n.yaw, nodes[i].yaw)
            self.assertAlmostEqual(n.pitch, nodes[i].pitch)
            self.assertAlmostEqual(n.roll, nodes[i].roll)

    # ---- setters -------------------------------------------------------------

    def test_setAnchor(self):
        m = pb.TurtlePolyline()
        m.setAnchor(pb.Vector3d(7.0, 8.0, 9.0))
        np.testing.assert_array_almost_equal(v3(m.getAnchor()), [7.0, 8.0, 9.0])

    def test_setAnchorFrame_changes_heading(self):
        """After setAnchorFrame(rotZ(pi/2)), forward moves along +y."""
        m = pb.TurtlePolyline()
        m.setAnchorFrame(pb.Matrix3d.rotZ(math.pi / 2))
        m.addNodeBack(2.0)
        np.testing.assert_array_almost_equal(v3(m.getNode(0)), [0.0, 2.0, 0.0], decimal=12)

    def test_setAnchorFrame_roundtrip(self):
        f = pb.Matrix3d.rotX(0.8)
        m = pb.TurtlePolyline()
        m.setAnchorFrame(f)
        np.testing.assert_array_almost_equal(m3(m.getAnchorFrame()), m3(f))

    # ---- str -----------------------------------------------------------------

    def test_str(self):
        m = pb.TurtlePolyline()
        m.addNodeBack(1.0)
        s = str(m)
        self.assertIn("TurtlePolyline", s)


if __name__ == "__main__":
    unittest.main(verbosity=2)
