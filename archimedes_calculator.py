import numpy as np
import matplotlib.pyplot as plt
import math

disk_r = 0.5
disk_A = math.pi*(disk_r*disk_r)
print("disk_area: ", disk_A)
disk_COM_y = 0.932596

tp1_y = 1.12613
tp2_y = 1.12584
# average height of triple points
tp_y = (tp1_y + tp2_y)/2
print("tp_y: ", tp_y)

dist_COM_to_tp_y = abs(disk_COM_y-tp_y)
print("dist_COM_to_tp_y: ", dist_COM_to_tp_y)

theta = math.acos(dist_COM_to_tp_y/disk_r)
segment_angle = 2*theta
print("theta: ", segment_angle)

segment_area = ((segment_angle-math.sin(segment_angle))/2)*(disk_r*disk_r)
print("segment_area: ", segment_area)
if disk_COM_y >= tp_y:
    submerged_area = segment_area
else:
    submerged_area = disk_A - segment_area
print("submerged_area: ", submerged_area)
submerged_percent = submerged_area/disk_A
print(submerged_percent)



