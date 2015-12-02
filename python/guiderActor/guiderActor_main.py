#!/usr/bin/env python
"""Start the guiderActor."""

from guiderActor import Msg, GuiderActor

# start a new actor
if __name__ == "__main__":
    #LCOHACK: since we don't have a proper domain name yet, 
    #LCOHACK: force the location to lco.
    guider = GuiderActor.GuiderActor.newActor(location='lco')
    guider.run(Msg=Msg)
