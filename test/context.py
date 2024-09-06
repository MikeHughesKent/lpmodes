# -*- coding: utf-8 -*-
"""
Adjusts path so the package can be found if not installed with pip

@author: Mike Hughes, Applied Optics Group, University of Kent
"""

import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
