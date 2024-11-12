import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from pathlib import Path

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

simFramesDir = ""
if len(sys.argv)>1 and "-outdir" in sys.argv:
    i = sys.argv.index("-outdir") 
    simFramesDir = sys.argv[i+1] + "/"
    outdir = outdir + simFramesDir + "/"
    Path('./outFrames/' + simFramesDir).mkdir(parents=True, exist_ok=True)

for frame in range(start, end):

    vList = []  # vertices
    vLiquidList = []
    fList = []  # faces
    vnList = [] # vertex normals
    vtList = [] # vertex tangents
    uList = [] # velocities
    vColorList = [] # vertice colours
    mpList = [] # marker particles
    #mpvList = [] # marker particle velocities
    rbCOM = []
    rbV_t = []
    rbTheta = []
    dt = 0
    outFreq = 1
    bgFieldPos = []
    bgFieldVel = []
    are = 0

    if os.path.isfile(outdir+str(frame)+".txt"):
        with open(outdir+str(frame)+".txt","r") as curr_frame_file:
            for row in curr_frame_file:
                r = row.split(" ")
                if r[0] == 'v':
                    vList.append([ float(r[1]), float(r[2])])
                    if len(r)>3:
                        if "c" in r[3]:
                            # vColorList.append('g')
                            vColorList.append('#219ebc')
                            vLiquidList.append([ float(r[1]), float(r[2])])
                        elif "a" in r[3]:
                            #vColorList.append('#3ec1d5')
                            vColorList.append('#fca311')
                            vLiquidList.append([ float(r[1]), float(r[2])])
                        elif "s" in r[3]:
                            #vColorList.append('#ff0070')
                            vColorList.append('#219ebc')
                            vLiquidList.append([ float(r[1]), float(r[2])])
                        elif "t" in r[3]:
                            #vColorList.append('#1b3481')
                            vColorList.append('#0f4c5c')
                            vLiquidList.append([ float(r[1]), float(r[2])])
                        else:
                            vColorList.append('k')
                    else:
                        vColorList.append('k')
                elif r[0] == 'f':
                    fList.append([ int(r[1]), int(r[2]) ])
                elif r[0] == 'vn':
                    vnList.append([ float(r[1]), float(r[2]) ])
                elif r[0] == 'vt':
                    vtList.append([ float(r[1]), float(r[2]) ])
                elif r[0] == 'u':
                    uList.append([ float(r[1]), float(r[2]) ])
                elif r[0] == 'p':
                    mpList.append([ float(r[1]), float(r[2]) ])
                #elif r[0] == 'pv':
                    #mpvList.append([ float(r[1]), float(r[2])])
                elif r[0] == 'dt':
                    dt = float(r[1])
                elif r[0] == 'out-freq':
                    outFreq = float(r[1])
                elif r[0] == 'rb':
                    rbCOM.append([ float(r[1]), float(r[2]) ])
                    rbTheta.append(float(r[3]))
                    rbV_t.append([ float(r[1]), float(r[2]) ])
                elif r[0] == 'bg':
                    bgFieldPos.append([ float(r[1]), float(r[2]) ])
                    bgFieldVel.append([ float(r[3]), float(r[4]) ])
                elif r[0] == 'area':
                    area = float(r[1])

        # why? I just like dealing with numpy :')
        # also, idk maybe one day will look into matplotlib .fill()
        v = np.array(vList)
        f = np.array(fList)
        vn = np.array(vnList)
        vt = np.array(vtList)
        u = np.array(uList)
        mp = np.array(mpList)
        #mpv = np.array(mpvList)

        for face in f:
            [indexStart, indexEnd] = face
            xPts = np.array([v[indexStart,0], v[indexEnd,0]])
            yPts = np.array([v[indexStart,1], v[indexEnd,1]])
            faceColor = "#005668"
            faceColor = "C0"
            if(vColorList[indexStart] == 'k' or vColorList[indexEnd] == 'k' ):
                faceColor = "#50a0c8"
                faceColor = "#bababa"
                plt.plot(xPts, yPts, faceColor,zorder=1)
            else:                
                plt.plot(xPts, yPts, faceColor,zorder=2)
        for i in range(len(rbCOM)):
            rMatrix = np.array([ [np.cos(rbTheta[i]), -np.sin(rbTheta[i])] , [np.sin(rbTheta[i]), np.cos(rbTheta[i])] ])
            north = 0.1 *  np.dot(rMatrix, np.array([0.0, 1.0]))
            plt.arrow(rbCOM[i][0], rbCOM[i][1], north[0], north[1], head_width=.05, color='k')
            east = 0.1 * np.dot(rMatrix, np.array([1.0, 0.0]))
            plt.arrow(rbCOM[i][0], rbCOM[i][1], east[0], east[1], head_width=.05, color='k')

        # set axes
        ax = plt.gca()
        # for wide cup
        #ax.set_xlim([-3, 3])
        #ax.set_ylim([-1.3, 2.0])
        # for big_up cup
        ax.set_xlim([-4.0,5.0])
        ax.set_ylim([-1.5, 2.5])
        ax.set_aspect('equal')

        if len(sys.argv)>1:
            if "-showpoints" in sys.argv:
                for j in range(len(v)):
                    vertex = v[j]
                    if vColorList[j] != 'k':
                        plt.plot(vertex[0], vertex[1],marker=".",color=vColorList[j])
            if "-showmarkers" in sys.argv:
                for j in range(len(mp)):
                    particle = mp[j]
                    #vel = mpv[j]
                    plt.plot(particle[0], particle[1], marker=".", markersize=1,color='#0a9396')
                    #plt.plot(particle[0], particle[1], marker=".", markersize=15,color='paleturquoise', zorder=0)
                    #plt.arrow(particle[0], particle[1], vel[0]*0.5, vel[1]*0.5, head_width=.05)
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
                plt.text(1.4, -1.9, "n: {}".format(len(vLiquidList)), fontsize = 11)
            if "-showfield" in sys.argv:
                for i in range(len(bgFieldPos)):
                    plt.arrow(bgFieldPos[i][0], bgFieldPos[i][1], bgFieldVel[i][0]*.5, bgFieldVel[i][1]*.5, head_width=.02, color="b")
        
        plt.gca().spines['top'].set_visible(False) 
        plt.gca().spines['right'].set_visible(False) 
        plt.gca().spines['bottom'].set_visible(False) 
        plt.gca().spines['left'].set_visible(False) 
        #plt.xticks([])
        #plt.yticks([])
        
        
        #plt.text(-2.9, -1.7, "t: {curr_t:.2f}".format(curr_t = outFreq*dt*frame), fontsize = 11)
        #plt.text(-2.9, -1.9, "Frame: {}".format(frame), fontsize = 11)
        #plt.text(2.0, -1.7, "t: {curr_t:.2f}".format(curr_t = outFreq*dt*frame), fontsize = 11)
        #plt.text(2.0, -1.7, "Frame: {}".format(frame), fontsize = 11)

        # for wide cup
        #plt.text(2.4, -1.2, "t: {curr_t:.2f}".format(curr_t = outFreq*dt*frame), fontsize = 11)
        # for extra wide cup
        plt.text(4.0, -1.2, "t: {curr_t:.2f}".format(curr_t = outFreq*dt*frame), fontsize = 11)
        # for big cup
        #plt.text(1.8, -1.7, "t: {curr_t:.2f}".format(curr_t = outFreq*dt*frame), fontsize = 11)

        # area
        #plt.text(4.0, -0.9, "area: {area:.2f}".format(area=area), fontsize = 11)

        #plt.show()
        plt.savefig('./outFrames/' + simFramesDir + 'frame-'+str(frame)+'.png', format="png", bbox_inches="tight")
        
        plt.clf()

