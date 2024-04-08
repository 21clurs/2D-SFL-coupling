import numpy as np
import matplotlib.pyplot as plt

# let's open one file first
# and just plot her!!

vList = []
fList = []

with open("square.txt","r") as curr_frame_file:
    for row in curr_frame_file:
        r = row.split(" ")
        if r[0] == 'v':
            vList.append([ float(r[1]), float(r[2])])
        elif r[0] == 'f':
            fList.append([ int(r[1]), int(r[2])])

# why? I just like dealing with numpy :')
# also, idk maybe one day will look into matplotlib .fill()
v = np.array(vList)
f = np.array(fList)

for face in f:
    [indexStart, indexEnd] = face
    xPts = np.array([v[indexStart,0], v[indexEnd,0]])
    yPts = np.array([v[indexStart,1], v[indexEnd,1]])
    plt.plot(xPts, yPts, '#50a0c8')

# set axes
ax = plt.gca()
ax.set_xlim([-2, 2])
ax.set_ylim([-2, 2])

#plt.show()
plt.savefig('./outFrames/'+'testing'+'.png')

