#!/usr/bin/env python
"""Start the guiderActor."""

from guiderActor import Msg, GuiderActor

# start a new actor
if __name__ == "__main__":
    guider = GuiderActor.GuiderActor.newActor(location='lco')
    guider.run(Msg=Msg)
