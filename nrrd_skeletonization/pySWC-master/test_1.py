# -*- coding: utf-8 -*-
"""
Scale and find surface area
"""
import pySWC

in_file = pySWC.Swc('test_file.swc')

in_file.halve_num_nodes()

in_file.save_file("test_scaled.swc")