# SIMPLACE-CPlantBox Docker Image

### Build the docker

```
docker build -t sp-cpb .
```

### Run the docker

With enabled DISPLAY
```
xhost +local:docker
docker run -it --rm -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix:rw sp-cpb bash
```

No DISPLAY (3D plots wont be shown iteractevly)
```
docker run -it sp-cpb
```

### Run the SIMPLACE-CPlantBox example script

```
cd ~/workspace/simplace_wrapper/python_simplace/trunk/PyPlantBox/misc/
python example_simplace_cplantbox.py
```

Output files (CSV and VTP) are saved in:
```
~/workspace/simplace_run/output/PyPlantBox/lintul/
```
