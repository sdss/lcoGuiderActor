"""
Test the PID class to show it works for various types of fake data.

The Kp, Ti, Td values should match those in the test data file headers. We've
picked Kp=0.6, Ti=100, Td=0 as our values for testing purposes (with Ti=0 for
the p-only cases), and the uniform dt=1.

The datafiles have the following columns:
    dt value p_corr pi_corr

where p_corr is the expected correction for only Kp, and pi_corr is with Kp and Ti.
"""

import unittest

import numpy as np

from guiderActor import PID


class TestPID(unittest.TestCase):
    def setUp(self):
        self.dt = 1
        self.Kp = 0.6
        self.Ti = 100
        self.Td = 0
        self.Imax = -1
        self.nfilt = 1
        self.names = ('dt', 'value', 'p_corr', 'pi_corr')
        self.dtype = np.dtype(zip(self.names, ['f8' for x in self.names]))
        self.linear_uniform_dt = 'data/linear_uniform_dt.txt'
        self.quadratic_uniform_dt = 'data/quadratic_uniform_dt.txt'
        self.linear_uniform_dt_noise = 'data/linear_uniform_dt_noise.txt'

    def test_p_linear_uniform_dt(self):
        """Test with only P term active, linear data, uniform dt."""
        pid = PID.PID(self.dt, self.Kp, 0, self.Td, self.Imax, self.nfilt)
        data = np.loadtxt(self.linear_uniform_dt, dtype=self.dtype)
        for line in data:
            correction = pid.update(line['value'])
            self.assertAlmostEqual(correction, line['p_corr'])

    def test_pi_linear_uniform_dt(self):
        """Test with both P and I terms active, linear data, uniform dt."""
        pid = PID.PID(self.dt, self.Kp, self.Ti, self.Td, self.Imax, self.nfilt)
        data = np.loadtxt(self.linear_uniform_dt, dtype=self.dtype)
        for line in data:
            correction = pid.update(line['value'])
            self.assertAlmostEqual(correction, line['pi_corr'])

    def test_p_quadratic_uniform_dt(self):
        """Test with only P term active, linear data, uniform dt."""
        pid = PID.PID(self.dt, self.Kp, 0, self.Td, self.Imax, self.nfilt)
        data = np.loadtxt(self.quadratic_uniform_dt, dtype=self.dtype)
        for line in data:
            correction = pid.update(line['value'])
            self.assertAlmostEqual(correction, line['p_corr'])

    def test_pi_quadratic_uniform_dt(self):
        """Test with only P term active, linear data, uniform dt."""
        pid = PID.PID(self.dt, self.Kp, self.Ti, self.Td, self.Imax, self.nfilt)
        data = np.loadtxt(self.quadratic_uniform_dt, dtype=self.dtype)
        for line in data:
            correction = pid.update(line['value'])
            self.assertAlmostEqual(correction, line['pi_corr'])

    def test_p_linear_uniform_dt_noise(self):
        """Test with only P term active, linear data, uniform dt."""
        pid = PID.PID(self.dt, self.Kp, 0, self.Td, self.Imax, self.nfilt)
        data = np.loadtxt(self.linear_uniform_dt_noise, dtype=self.dtype)
        for line in data:
            correction = pid.update(line['value'])
            self.assertAlmostEqual(correction, line['p_corr'])

    def test_pi_linear_uniform_dt_noise(self):
        """Test with both P and I terms active, linear data, uniform dt."""
        pid = PID.PID(self.dt, self.Kp, self.Ti, self.Td, self.Imax, self.nfilt)
        data = np.loadtxt(self.linear_uniform_dt_noise, dtype=self.dtype)
        for line in data:
            correction = pid.update(line['value'])
            self.assertAlmostEqual(correction, line['pi_corr'])


if __name__ == '__main__':
    unittest.main(verbosity=1)
