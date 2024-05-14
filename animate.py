import numpy as np
import matplotlib.pyplot as plt
import sys
import os

outdir = "./sim/out/"
frames = os.listdir(outdir)
for frame in range(30):

    vList = []  # vertices
    fList = []  # faces
    vnList = [] # vertex normals
    vtList = [] # vertex tangents

    if os.path.isfile(outdir+str(frame)+".txt"):
        with open(outdir+str(frame)+".txt","r") as curr_frame_file:
            for row in curr_frame_file:
                r = row.split(" ")
                if r[0] == 'v':
                    vList.append([ float(r[1]), float(r[2])])
                elif r[0] == 'f':
                    fList.append([ int(r[1]), int(r[2])])
                elif r[0] == 'vn':
                    vnList.append([ float(r[1]), float(r[2])])
                elif r[0] == 'vt':
                    vtList.append([ float(r[1]), float(r[2])])

        # why? I just like dealing with numpy :')
        # also, idk maybe one day will look into matplotlib .fill()
        v = np.array(vList)
        f = np.array(fList)
        vn = np.array(vnList)
        vt = np.array(vtList)

        for face in f:
            [indexStart, indexEnd] = face
            xPts = np.array([v[indexStart,0], v[indexEnd,0]])
            yPts = np.array([v[indexStart,1], v[indexEnd,1]])
            plt.plot(xPts, yPts, '#50a0c8')
        
        # set axes
        ax = plt.gca()
        ax.set_xlim([-2, 2])
        ax.set_ylim([-2, 2])
        ax.set_aspect('equal')

        if len(sys.argv)>1:
            if "-shownorms" in sys.argv:
                for i in range(len(vn)):
                    plt.arrow(v[i,0], v[i,1], vn[i,0]*0.2, vn[i,1]*0.2, head_width=.05)
            if "-showtangents" in sys.argv:
                for i in range(len(vn)):
                    plt.arrow(v[i,0], v[i,1], vt[i,0]*0.2, vt[i,1]*0.2, head_width=.05, color="r")

        # plt.show()
        plt.savefig('./outFrames/frame-'+str(frame)+'.png', format="png")
        
        plt.clf()

