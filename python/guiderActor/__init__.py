import Queue

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#
# Queue names
#
MASTER = 0
GCAMERA = 1

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class Msg(object):
    # Priorities
    CRITICAL = 0
    HIGH = 2
    MEDIUM = 4
    NORMAL = 6
    
    # Command types
    START_GUIDING = 0
    EXPOSE = 1
    EXPOSURE_FINISHED = 2
    ENABLE_FIBER = 3
    LOAD_CARTRIDGE = 4
    SET_GUIDE_MODE = 5
    STATUS = 6
    ABORT_EXPOSURE = 7
    SET_TIME = 8

    def __init__(self, type, data=None, priority=NORMAL):
        self.type = type
        self.data = data
        self.priority = priority

    def __str__(self):
        return "%s : %s" % (self.type, self.data)

    def __repr__(self, rhs):
        return self.priority - rhs.priority

def flushQueue(queue):
    """flush queue"""
    
    while True:     
        try:
            msg = queue.get(timeout=0)
        except Queue.Empty:
            return

__all__ = ["MASTER", "GCAMERA", "Msg"]
