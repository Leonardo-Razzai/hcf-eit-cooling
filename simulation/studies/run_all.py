"""Regenerate every result file of the paper (~55 min in total)."""
import time

import calibration
import dynamics
import lines
import robustness
import saturation
import temperature
from _common import log

if __name__ == '__main__':
    for name, fn in (('saturation', saturation.run),
                     ('robustness', lambda: [f() for f in robustness.STUDIES.values()]),
                     ('calibration', lambda: [f() for f in calibration.STUDIES.values()]),
                     ('dynamics', lambda: [f() for f in dynamics.STUDIES.values()]),
                     ('temperature', lambda: [f() for f in temperature.STUDIES.values()]),
                     ('lines', lambda: [f() for f in lines.STUDIES.values()])):
        t0 = time.time()
        fn()
        log(name, 'done in %.0f s' % (time.time() - t0))
