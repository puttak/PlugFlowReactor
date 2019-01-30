# -*- coding: utf-8 -*-
""" Provides a Python interface to PlugFlowReactor model.

Notice that this interface does not currently implement all the possible use
cases of PlugFlowReactor C++ class. This might be fixed in the future. Except
for the C shared library itself, there is no need to compile any Cython code,
meaning that this package is based simply on built-in ctypes features and thus
is submitted to its limitations.
"""

import os
import math
import time
import ctypes
import numpy
from datetime import datetime
from matplotlib import pyplot
from pandas import read_csv


class PyPlugFlowReactor(object):
    """ Interface to PFR solver.

    Parameters
    ----------

    """

    def __init__(self, mech, phase, D, T, P, X, mdot):
        self._load_library()
        self._convert_types(mech, phase, D, T, P, X, mdot)

    def _load_library(self):
        """ Load C-interface function and add mechanisms to path. """
        lib = 'libPlugFlowReactor_shared.so'
        mdir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
        os.environ['CANTERA_DATA'] = os.path.join(mdir, 'data')
        self._clib = ctypes.CDLL(os.path.join(mdir, lib))
        self._func = self._clib['C_PlugFlowReactor']

    def _convert_types(self, mech, phase, D, T, P, X, mdot):
        """ Convert types for C-compatibility.

        TODO
        ----
        - Process other types for X, such as dictionary and list.
        """
        self._mech = ctypes.create_string_buffer(mech.encode('utf-8'))
        self._phase = ctypes.create_string_buffer(phase.encode('utf-8'))
        self._X = ctypes.create_string_buffer(X.encode('utf-8'))

        self._D = ctypes.c_double(D)
        self._T = ctypes.c_double(T)
        self._P = ctypes.c_double(P)
        self._mdot = ctypes.c_double(mdot)

    def _check_init(self):
        """ Check if model has been fully initialized. """
        try:
            msg0 = 'Missing convective heat transfer coefficient'
            msg1 = 'Missing wall temperature profile function'
            assert hasattr(self, '_htc'), msg0
            assert hasattr(self, '_cctwall'), msg1
        except AssertionError as err:
            print(f'Fallback to missing parameter:\n{err}')
            self.set_wall_conditions(0.0)

        if not hasattr(self, '_rtol') or not hasattr(self, '_atol'):
            self.set_tolerances()

        if not hasattr(self, '_initstep') or not hasattr(self, '_maxsteps'):
            self.set_steps()

    def set_tolerances(self, rtol=1.0e-06, atol=1.0e-15):
        """ Set relative and absolute tolerances.

        Parameters
        ----------
        rtol : float, optional
            Relative tolerance for solver. Default is `1.0e-06`.
        atol : float, optional
            Absolute tolerance for solver. Default is `1.0e-15`.

        Raises
        ------
        ValueError
            If any of the provided tolerances is non-positive.
        """
        try:
            assert rtol > 0, f'Found non-positive relative tolerance : {rtol}'
            assert atol > 0, f'Found non-positive absolute tolerance : {atol}'
        except AssertionError as err:
            raise ValueError(err)

        self._rtol = ctypes.c_double(rtol)
        self._atol = ctypes.c_double(atol)

    def set_steps(self, initstep=1.0e-05, maxsteps=10000):
        """ Set initial step and maximum number of steps.

        Parameters
        ----------

        Raises
        ------
        """
        self._initstep = ctypes.c_double(initstep)  # FIXME Unused!
        self._maxsteps = ctypes.c_uint(maxsteps)

    def set_wall_conditions(self, htc, tw=None):
        """ Set heat transfer coefficient and wall temperature.

        This method can build the wall temperature in some different manners.
        First, if wall temperature is not provided, the reactor wall will be
        assumed at inlet gas temperature. If a numeric wall temperature is
        provided or the previous case is used, this is employed to build a
        constant wall temperature function in the format compatible with the
        underlining library.

        Note
        ----
        It is important to keep a local class pointer to the original function
        so that garbage collector does not cause a segmentation fault by
        suppressing it.

        TODO
        ----
        - Parse sympy compatible strings!
        - Accept pure C-functions!

        Parameters
        ----------

        Raises
        ------
        """
        if math.isclose(htc, 0.0):
            print('Assuming adiabatic system')

        self._htc = ctypes.c_double(htc)

        tw = float(self._T.value) if tw is None else tw
        if isinstance(tw, (ctypes.c_double, float, int)):
            self._twptr = tw
            self._pytwall = lambda x: self._twptr
        else:
            self._pytwall = tw

        ftype = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
        self._cctwall = ftype(lambda x: self._pytwall(x))

    def solve(self, L, dx=None, saveas=None):
        """ Integrate problem.

        Parameters
        ----------

        Raises
        ------
        """
        self._nowstr = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        self._check_init()

        if dx is None:
            dx = ctypes.c_double(L / 100)

        if saveas is None:
            saveas = f'results-from-{self._nowstr}.csv'

        self._saveas = saveas
        saveas = ctypes.create_string_buffer(saveas.encode('utf-8'))
        length = ctypes.c_double(L)

        t0 = time.time()
        cresult = self._func(self._mech, self._phase, self._T, self._P,
                             self._X, self._D, length, self._mdot, self._htc,
                             self._cctwall, saveas, dx, self._rtol,
                             self._atol, self._maxsteps)
        print(f'Calculation took {time.time()-t0} s\n\n')
        return cresult

    def plot_species(self, vars, saveas=None, labels=None, **kwargs):
        if isinstance(vars, str):
            vars = [vars]
        vars = list(set(vars))

        if labels is not None:
            if isinstance(labels, str):
                labels = [labels]
            labels = list(set(labels))

        if saveas is None:
            num = 0
            while True:
                saveas = f'fig-from-{self._nowstr}-{num}.png'
                if not os.path.exists(saveas):
                    break

        data = read_csv(self._saveas, usecols=['x'] + vars)

        pyplot.close('all')
        pyplot.style.use('bmh')

        for i, v in enumerate(vars):
            l = v if labels is None else labels[i]
            pyplot.plot(data['x'], data[v], label=l)

        pyplot.xlabel('Position (m)')
        pyplot.ylabel('Mass fraction')
        pyplot.legend()
        pyplot.savefig(saveas, dpi=kwargs.get('dpi', 300))


def test_full():
    """ Print test several use cases. """
    mech = "CT-hydrocarbon-dalmazsi-2017-mech.xml"
    phase = "gas"
    L = 0.45
    D = 0.028
    T = 300.0
    P = 5000.0
    X = "N2:0.64, C2H2:0.36"
    mdot = 1.0e-05
    htc = 10.0

    def Tw(x):
        """ Wall temperature in terms of position. """
        Ta, Tc, Ts = 300.0, 1173.0, 400.0
        x1, x2, m1, m2 = 0.02492942, 0.40810172, 0.78913918, 11.91548263
        term1 = 1 - numpy.exp(-(x / x1) ** m1)
        term2 = 1 - numpy.exp(-(x / x2) ** m2)
        wallT = Ta + (Tc - Ta) * term1 - (Tc - Ts) * term2
        return 0.97 * wallT

    r = PyPlugFlowReactor(mech, phase, D, T, P, X, mdot)

    print(79 * '*')
    print('Convective with wall profile function')
    print(79 * '*')
    r.set_wall_conditions(htc, tw=Tw)
    r.solve(L)
    r.plot_species('C2H2')

    print(79 * '*')
    print('Convective with wall temperature as inlet')
    print(79 * '*')
    r.set_wall_conditions(htc, tw=None)
    r.solve(L)

    print(79 * '*')
    print('Convective with fixed wall temperature')
    print(79 * '*')
    r.set_wall_conditions(htc, tw=900.0)
    r.solve(L)

    print(79 * '*')
    print('Adiabatic and wall temperature does not matter')
    print(79 * '*')
    r.set_wall_conditions(0.0, tw=None)
    r.solve(L)


def main():
    """Call to main test."""
    test_full()


if __name__ == '__main__':
    main()
