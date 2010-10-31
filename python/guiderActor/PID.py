import numpy

class PID(object):
    """A class to handle PID loops"""
    def __init__(self, dt, Kp, Ti, Td, Imax=-1, nfilt=1):
        self.dt = dt                    # Time between PID updates
        self.Kp = Kp                    # Proportional term
        self.Ti = Ti                    # Integral time
        self.Td = Td                    # Derivative time

        self.Imax = Imax                # limit to abs(.Ix)
        self.nfilt = nfilt              # number of inputs to smooth over.
        self.histlen = 0
        
        self._x = None                  # previous value of the error, x
        
    def __str__(self):
        return "K_P = %g  T_i = %g  T_d = %g" % (self.Kp, self.Ti, self.Td)

    def filterX(self, x):
        """ Apply some smoothing/predictive filter to inputs. Dumb median now; kalman maybe later.
        No adjustements being made to dt/Ti/Td.
        """ 

        if len(self._xhist) > 1:
            self._xhist[:-1] = self._xhist[1:]
        self._xhist[-1] = x

        self.histlen = min(self.histlen+1, len(self._xhist))
        return numpy.median(self._xhist[-self.histlen:])
    
    def update(self, x):
        """GIven a new sample of the error, x, return the PID update"""
        
        if self._x is None:
            self.Dx = 0
            self.Ix = 0
            self._xhist = numpy.zeros(self.nfilt, dtype='f4')
            x = self.filterX(x)
        else:
            x = self.filterX(x)
            self.Dx = (x - self._x)/self.dt
            self.Ix += self.dt*x

        if self.Imax > 0 and abs(self.Ix) > self.Imax:
            self.Ix = self.Imax if self.Ix > 0 else -self.Imax

        self._x = x

        correction = self._x + self.Td*self.Dx
        if self.Ti:
            correction += self.Ix/self.Ti
            
        return self.Kp*correction

    def reset(self, dt=None):
        """Reset the PID loop, in particular the I term. Optionally reset the sampling rate. """
        self._x = None
        if dt:
            self.dt = dt

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

    def setPID(self, dt=None, Kp=None, Ti=None, Td=None, Imax=None, nfilt=None):

        needReset = dt != None or Imax != None or nfilt != None

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

        if nfilt:
            self.nfilt = nfilt

        if needReset:
            self.reset()
