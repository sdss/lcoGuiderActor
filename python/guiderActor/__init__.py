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
    
    # Command types; use classes so that the unique IDs are automatically generated
    class START_GUIDING(): pass
    class EXPOSE(): pass
    class EXPOSURE_FINISHED(): pass
    class ENABLE_FIBER(): pass
    class LOAD_CARTRIDGE(): pass
    class SET_GUIDE_MODE(): pass
    class STATUS(): pass
    class ABORT_EXPOSURE(): pass
    class SET_TIME(): pass

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
