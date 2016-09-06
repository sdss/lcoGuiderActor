"""
A relatively basic proportional-integral-derivative controller, used to compute
the axis corrections to apply to the telescope, using the errors computed from
the guider frames.
"""

import numpy

class PID(object):
    """A class to handle PID loops"""
    def __init__(self, dt, Kp, Ti, Td, Imax=-1, nfilt=1, tfilt=None, ncorr=1):

        self.dt = dt                    # Time between PID updates
        self.Kp = Kp                    # Proportional term
        self.Ti = Ti                    # Integral time
        self.Td = Td                    # Derivative time

        # NOTE: TODO: Imax is currently unused!
        self.Imax = Imax                # limit to abs(.Ix)
        self.tfilt = tfilt
        if tfilt != None:
            self.nfilt = tfilt / self.dt
        else:
            self.nfilt = nfilt              # number of inputs to smooth over.

        # The number of iterations after which a correction will be applied.
        self.ncorr = ncorr
        self.corr_count = 0  # A counter that updates each iteration.

        self.histlen = 0

        self._x = None                  # previous value of the error, x

    def __str__(self):
        return ('K_P={0:g} T_i={1:g} T_d={2:g} tfilt={3} nfilt={4:d} ncorr={5:d}'
                .format(self.Kp, self.Ti, self.Td, self.tfilt, self.nfilt, self.ncorr))

    def filterX(self, x):
        """ Apply some smoothing/predictive filter to inputs. Dumb median now; kalman maybe later.
        No adjustements being made to dt/Ti/Td.
        """

        if len(self._xhist) > 1:
            self._xhist[:-1] = self._xhist[1:]
        self._xhist[-1] = x

        self.histlen = min(self.histlen+1, len(self._xhist))
        return numpy.median(self._xhist[-self.histlen:])

    def update(self, x, dt = None):
        """GIven a new sample of the error, x, return the PID update"""

        if dt is None:
            dt = self.dt

        if self._x is None:
            self.Dx = 0
            self.Ix = 0
            self._xhist = numpy.zeros(self.nfilt, dtype='f4')
            x = self.filterX(x)
        else:
            x = self.filterX(x)
            self.Dx = (x - self._x)/dt
            self.Ix += dt*x / self.Ti # scale by Ti now

        # NOTE: TODO: ignoring this because we are scaling Ti with altitude,
        # so to do this right, we'd have to think about how to use Imax in that case.
        # if self.Imax > 0 and abs(self.Ix) > self.Imax:
        #     self.Ix = self.Imax if self.Ix > 0 else -self.Imax

        self._x = x

        self.corr_count += 1

        correction = self._x + self.Td*self.Dx
        if self.Ti:
            correction += self.Ix

        return self.Kp*correction

    def reset(self, dt=None):
        """Reset the PID loop, in particular the I term. Optionally reset the sampling rate. """
        self._x = None
        if dt:
            self.dt = dt
            if self.tfilt:
                self.nfilt = self.tfilt / dt

        # NOTE: I'm not completely sure if this is ok, but I think if we reset the loop
        # we also want to reset the counter (JSG)
        self.corr_count = 0

    def ZieglerNichols(self, Kpc, Pc, loopType="PID"):
        """Perform Ziegler-Nichols tuning of a PID loop; Kpc is the critical proportional
        gain that just starts to oscillate with period Pc"""

        Pc = float(Pc)

        if loopType == "P":
            self.Kp = 0.50*Kpc
            self.Ti = self.Td = 0
        elif loopType == "PI":
            self.Kp = 0.45*Kpc
            self.Ti = Pc/1.2
            self.Td = 0
        elif loopType == "PID":
            self.Kp = 0.60*Kpc
            self.Ti = Pc/2
            self.Td = Pc/8
        else:
            raise RuntimeError, ("I don't know how to tune a %s loop" % loopType)

    def setPID(self, dt=None, Kp=None, Ti=None, Td=None, Imax=None,
               nfilt=None, tfilt=None, ncorr=None):

        needReset = (not dt or not Imax or not nfilt or not tfilt or not ncorr)

        if dt is not None:
            try:
                1/dt
            except ZeroDivisionError:
                raise RuntimeError, "You may not specify a time sample of 0s"

            self.dt = dt

        if Kp is not None:
            self.Kp = Kp

        if Ti is not None:
            self.Ti = Ti

        if Td is not None:
            self.Td = Td

        if Imax is not None:
            self.Imax = Imax

        if tfilt:
            self.tfilt = tfilt
            self.nfilt = tfilt / self.dt
        elif nfilt:
            self.nfilt = nfilt

        if ncorr is not None:
            self.ncorr = ncorr

        if needReset:
            self.reset()
