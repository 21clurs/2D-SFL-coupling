A 2D implementation of [Surface-Only Liquids](http://www.cs.columbia.edu/cg/surfaceliquids/)
## Building
```
cd sim
make
```
# Running
Starting in `sim` directory:
```
./sim.out
cd ..
python animate.py
cd outFrames
ffmpeg -i frame-%d.png video.mp4
```

Here we use many ffmpeg defaults (i.e. framerate is 25, but of course these are changeable). I'll one day make this so that it isn't so silly...

If you want a gif instead for whatever reason, I like to first generate a colour palette for it, as otherwise the colours get a bit funky as a result of the conversions:
```
ffmpeg -i frame-%d.png -vf palettegen palette.png
ffmpeg -i frame-%d.png -i palette.png -lavfi paletteuse animation.gif
```
