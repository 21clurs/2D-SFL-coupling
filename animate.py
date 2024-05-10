import numpy as np
import matplotlib.pyplot as plt
import sys
import os

outdir = "./sim/out/"
frames = os.listdir(outdir)
for i in range(len(frames)):

    vList = []
    fList = []
    vnList = []

    if os.path.isfile(outdir+str(i)+".txt"):
        with open(outdir+str(i)+".txt","r") as curr_frame_file:
            for row in curr_frame_file:
                r = row.split(" ")
                if r[0] == 'v':
                    vList.append([ float(r[1]), float(r[2])])
                elif r[0] == 'f':
                    fList.append([ int(r[1]), int(r[2])])
                elif r[0] == 'vn':
                    vnList.append([ float(r[1]), float(r[2])])

        # why? I just like dealing with numpy :')
        # also, idk maybe one day will look into matplotlib .fill()
        v = np.array(vList)
        f = np.array(fList)
        vn = np.array(vnList)

        for face in f:
            [indexStart, indexEnd] = face
            xPts = np.array([v[indexStart,0], v[indexEnd,0]])
            yPts = np.array([v[indexStart,1], v[indexEnd,1]])
            plt.plot(xPts, yPts, '#50a0c8')
        if len(sys.argv)>1 and sys.argv[1]=="-shownorms":
            for i in range(len(vn)):
                plt.arrow(v[i,0], v[i,1], vn[i,0]*0.2, vn[i,1]*0.2, head_width=.05)

        # set axes
        ax = plt.gca()
        ax.set_xlim([-2, 2])
        ax.set_ylim([-2, 2])
        ax.set_aspect('equal')

        #plt.show()
        plt.savefig('./outFrames/frame-'+str(i)+'.png')

        plt.clf()

