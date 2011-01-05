import math
import numpy
import os
import subprocess
import tempfile
from opscore.utility.YPF import YPF

from guiderActor import *

def getDefGprobes(plateId, cartId):
	if not 'GUIDERACTOR_DIR' in os.environ:
		raise RuntimeError('guiderActor product must be setup')
	if not 'PLATEDB_DIR' in os.environ:
		raise RuntimeError('platedb product must be setup')

	fiberinfoFile = os.path.join(os.environ.get('GUIDERACTOR_DIR'), 'etc', 'gcamFiberInfo.par')
	scratchFile = tempfile.NamedTemporaryFile()

	args = ['catPlPlugMapM']
	if cartId:
		args.extend(['-c', str(cartId)])
	args.append(str(plateId))
	#args.extend(['>>', scratchFile.name])
	subprocess.check_call(args, stdout=scratchFile)

	return getGprobes(fiberinfoFile, scratchFile.name, cartId)

def getGprobes(fiberinfofn, plugmapfn, cartId):
	par = YPF(fiberinfofn)
	probes = par.structs['GPROBE'].asObjlist()
	probes = [p for p in probes if p.cartridgeId == cartId]
	guideprobes = [p for p in probes if p.fiberType == 'GUIDE' or p.fiberType == 'ACQUIRE']
	tritiumprobes = [p for p in probes if p.fiberType == 'TRITIUM']
        
	plugmapfile = YPF(plugmapfn)
	plugmap = plugmapfile.structs['PLUGMAPOBJ'].asObjlist()
	aholes = [p for p in plugmap if p.holeType == 'ALIGNMENT']
	gholes = [p for p in plugmap if p.holeType == 'GUIDE']

	print '%i probes' % len(probes)
	print '  %i guide/capture probes' % len(guideprobes)
	print '  %i tritium probes' % len(tritiumprobes)
	print '%i alignment holes' % len(aholes)
	print '%i guide holes' % len(gholes)

	gprobes = {}
	for (p, ghole, ahole) in zip(guideprobes + tritiumprobes, gholes + [None], aholes + [None]):
		info = GprobeInfo(p.exists, True,
                                  p.xcen, p.ycen, p.radius, p.rot,
                                  p.xferruleOffset, p.yferruleOffset,
                                  p.focusOffset, p.fiberType, 0)
                
		info.enabled = True
		info.exists = p.exists
		info.fiber_type = p.fiberType
		info.flags = 0
		info.xCenter = p.xcen
		info.yCenter = p.ycen
		info.radius = p.radius
		info.xFerruleOffset = p.xferruleOffset
		info.yFerruleOffset = p.yferruleOffset
		info.rotation = p.rot
		info.focusOffset = p.focusOffset
		if p.fiberType in ['GUIDE', 'ACQUIRE']:
			info.ra  = ghole.ra
			info.dec = ghole.dec
			info.xFocal = ghole.xFocal
			info.yFocal = ghole.yFocal
			info.phi = 90 - math.atan2(ahole.yFocal - ghole.yFocal,
                                                   ahole.xFocal - ghole.xFocal) * 180/math.pi
			info.mag = ghole.mag

		else:
			info.ra = 0
			info.dec = 0
			info.xFocal = 0
			info.yFocal = 0
			info.phi = 0
			info.rotStar2Sky = numpy.nan


		gprobes[p.gProbeId] = info
	return gprobes

if __name__ == "__main__":
	pass

