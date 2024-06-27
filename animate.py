import numpy as np
import matplotlib.pyplot as plt
import sys
import os

outdir = "./sim/out/"
frames = os.listdir(outdir)
start = 0
end = len(frames)

if len(sys.argv)>1 and "-framerange" in sys.argv:
    i = sys.argv.index("-framerange")
    try:
        start = int(sys.argv[i+1])
        end = int(sys.argv[i+2])
    except IndexError:
        print("No valid range of frames entered. Retrieving all frames...")

for frame in range(start, end):

    vList = []  # vertices
    fList = []  # faces
    vnList = [] # vertex normals
    vtList = [] # vertex tangents
    uList = [] # velocities
    vColorList = [] # vertice colours
    mpList = [] # marker particles

    if os.path.isfile(outdir+str(frame)+".txt"):
        with open(outdir+str(frame)+".txt","r") as curr_frame_file:
            for row in curr_frame_file:
                r = row.split(" ")
                if r[0] == 'v':
                    vList.append([ float(r[1]), float(r[2])])
                    if len(r)>3:
                        if "c" in r[3]:
                            vColorList.append('g')
                        elif "a" in r[3]:
                            vColorList.append('#3ec1d5')
                        elif "s" in r[3]:
                            vColorList.append('#ff0070')
                        elif "t" in r[3]:
                            vColorList.append('#1b3481')
                        else:
                            vColorList.append('k')
                    else:
                        vColorList.append('k')
                elif r[0] == 'f':
                    fList.append([ int(r[1]), int(r[2])])
                elif r[0] == 'vn':
                    vnList.append([ float(r[1]), float(r[2])])
                elif r[0] == 'vt':
                    vtList.append([ float(r[1]), float(r[2])])
                elif r[0] == 'u':
                    uList.append([ float(r[1]), float(r[2])])
                elif r[0] == 'p':
                    mpList.append([ float(r[1]), float(r[2])])

        # why? I just like dealing with numpy :')
        # also, idk maybe one day will look into matplotlib .fill()
        v = np.array(vList)
        f = np.array(fList)
        vn = np.array(vnList)
        vt = np.array(vtList)
        u = np.array(uList)
        mp = np.array(mpList)

        for face in f:
            [indexStart, indexEnd] = face
            xPts = np.array([v[indexStart,0], v[indexEnd,0]])
            yPts = np.array([v[indexStart,1], v[indexEnd,1]])
            faceColor = "#005668"
            if(vColorList[indexStart] == 'k' or vColorList[indexEnd] == 'k' ):
                faceColor = "#50a0c8"
            plt.plot(xPts, yPts, faceColor)

        if "-showpoints" in sys.argv:
            for j in range(len(v)):
                vertex = v[j]
                if vColorList[j] != 'k':
                    plt.plot(vertex[0], vertex[1],marker=".",color=vColorList[j])
        
        if "-showmarkers" in sys.argv:
            for j in range(len(mp)):
                particle = mp[j]
                plt.plot(particle[0], particle[1], marker=".",color='r')

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
            if "-showvels" in sys.argv:
                for i in range(len(u)):
                    plt.arrow(v[i,0], v[i,1], u[i,0]*0.5, u[i,1]*0.5, head_width=.05)
            if "-shownumpoints" in sys.argv:
                plt.text(1.4, -1.9, "n: {}".format(len(v)), fontsize = 11)
            
        plt.text(-1.9, -1.9, "Frame: {}".format(frame), fontsize = 11)

        #plt.show()
        plt.savefig('./outFrames/frame-'+str(frame)+'.png', format="png")
        
        plt.clf()

