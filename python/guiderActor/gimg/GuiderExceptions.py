__all__ = ['GuiderError','FlatError','NoFibersFoundError','BadReadError']

# for guider-related exception handling
class GuiderError(Exception):
    """A general exception related to problems processing guider data."""
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
