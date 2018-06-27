.. _lcoGuiderActor-changelog:

==========
Change Log
==========

This document records the main changes to the lcoGuiderActor code.


.. _changelog-1.0.2:

1.0.2 (unreleased)
------------------

Added
^^^^^
* Reimplemented ``centerUp`` to use any available guide fibre. In this case, the full offset correction is applied.
* ``guideStep`` now executes all the offsets at once making use of ``tcc guideoffset``.

Changed
^^^^^^^
* ``CCDInfo`` initialised with bias level zero given that the image passed has already been bias subtracted. Should address the problem with LED jitter.

Refactored
^^^^^^^^^^
* Applied ``autopep8``, ``isort``, and ``unify`` to all files.
* Removed references to APO.
* Removed deprecated files.

Fixed
^^^^^
* Ticket `#1433 <https://trac.sdss.org/ticket/1433>`__: lack of plateGuideOffset files should be fatal for APOGEE.

`View changes <https://github.com/sdss/lcoGuiderActor/compare/1.0.1...HEAD>`__


.. _changelog-1.0.1:

1.0.1 (2018-01-25)
------------------

Added
^^^^^
* Stores the flexures detected using the LEDs as ``xFlex`` and ``yFlex`` in the binary table of the ``proc-gimg`` images.
* LED ``Fiber`` now gets the PSF-related fields (``fwhm``, ``flux``, ``sky``) populated.
* In frame image, the LED stamp uses the saturated mask when calling PyGuide.
* Cleaned up the LED code a bit.
* ``guider flat`` now accepts a ``force`` flag. If not set and the flat cannot detect at least one LED, the flat will fail. If set, it will warn but continue.

Removed
^^^^^^^
* The ``checkTritium()`` function in ``GProbe``.
* Lots of obsolete files.

`View changes <https://github.com/sdss/lcoGuiderActor/compare/1.0.1...1.0.0>`__


.. _changelog-1.0.0:

1.0.0 (2017-12-18)
-------------------

Added
^^^^^
* Initial version of the guiderActor fork for LCO.
* Bumpversion and new version attribute.

Removed
^^^^^^^
* The ``setup.py`` file that is of no use right now.
