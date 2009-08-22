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
    class FAIL(): pass
    class LOAD_CARTRIDGE(): pass
    class SET_GUIDE_MODE(): pass
    class STATUS(): pass
    class ABORT_EXPOSURE(): pass
    class SET_SCALE(): pass
    class SET_TIME(): pass

    def __init__(self, type, cmd, **data):
        self.type = type
        self.cmd = cmd
        self.priority = Msg.NORMAL      # may be overridden by **data
        #
        # convert data[] into attributes
        #
        for k, v in data.items():
            self.__setattr__(k, v)
        self.__data = data.keys()

    def __repr__(self):
        values = []
        for k in self.__data:
            values.append("%s : %s" % (k, self.__getattribute__(k)))

        return "%s, %s: {%s}" % (self.type.__name__, self.cmd, ", ".join(values))

    def __cmp__(self, rhs):
        """Used when sorting the messages in a priority queue"""
        return self.priority - rhs.priority

def flushQueue(queue):
    """flush queue"""
    
    while True:     
        try:
            msg = queue.get(timeout=0)
        except Queue.Empty:
            return

__all__ = ["MASTER", "GCAMERA", "Msg"]