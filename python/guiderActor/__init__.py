import actorcore.Actor

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#
# Queue names
#
MASTER = 0
GCAMERA = 1
MOVIE = 2

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class Msg(actorcore.Actor.Msg):
    # Priorities
    CRITICAL = 0
    HIGH = 2
    MEDIUM = 4
    NORMAL = 6

    # Command types; use classes so that the unique IDs are automatically generated
    class START_GUIDING(): pass
    class STOP_GUIDING(): pass
    class EXPOSE(): pass
    class EXPOSURE_FINISHED(): pass
    class TAKE_FLAT(): pass
    class TAKE_BIAS(): pass
    class FLAT_FINISHED(): pass
    class TAKE_DARK(): pass
    class DARK_FINISHED(): pass
    class BIAS_FINISHED(): pass
    class ENABLE_FIBER(): pass
    class LOAD_CARTRIDGE(): pass
    class SET_GUIDE_MODE(): pass
    class STATUS(): pass
    class ABORT_EXPOSURE(): pass
    class SET_PID(): pass
    class SET_SCALE(): pass
    class SET_SPECIAL_GPROBES(): pass
    class SET_TIME(): pass
    class CHANGE_SCALE(): pass
    class REPROCESS_FILE(): pass
    class READ_PLATE_FILES(): pass
    class DECENTER(): pass
    class CENTERUP(): pass
    class SET_REFRACTION(): pass
    class STAR_IN_FIBER(): pass
    class MAKE_MOVIE(): pass
    class ONESTEP(): pass

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

__all__ = ["MASTER", "GCAMERA", "MOVIE", "Msg"]

__version__ = '1.0.1dev'
