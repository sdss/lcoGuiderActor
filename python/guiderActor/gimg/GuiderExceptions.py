from opscore.utility.qstr import qstr


__all__ = ['GuiderError', 'BadDarkError', 'FlatError', 'NoFibersFoundError',
           'BadReadError', 'raiseGuiderError', 'BadBiasError']


# for guider-related exception handling
class GuiderError(Exception):
    """A general exception related to problems processing guider data."""
    pass


class BadDarkError(GuiderError):
    """An exception due to a problem processing a dark."""
    pass


class BadBiasError(GuiderError):
    """An exception due to a problem processing a bias."""
    pass


class FlatError(GuiderError):
    """An exception due to a problem processing a flat."""
    pass


class NoFibersFoundError(FlatError):
    """An exception due to not finding any fibers in an image."""
    pass


class BadReadError(GuiderError):
    """An exception due to an incorrectly read guider image."""
    pass


def raiseGuiderError(cmd, msg):
    """Convenience function to issue `cmd.error` and raise the error."""

    cmd.error(qstr(msg))

    raise GuiderError(msg)
