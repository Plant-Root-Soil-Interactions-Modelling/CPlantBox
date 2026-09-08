# SIMPLACE-CPlantBox

To run the coupled example `../example7_5_coupling_pbcm.py`, both SIMPLACE (v5.0) and CPlantBox must be properly installed in the system. For simplicity, we prepared a `docker` image with all the dependencies. 

Please make sure to have `docker` installed and follow the steps below.

## Run the example from docker

Navigate to `simplace_cpb/` and run the command below to lauch the `docker` container with display enabled 
```bash
xhost +local:docker
docker run -it --rm \
  -e DISPLAY=$DISPLAY \
  -v /tmp/.X11-unix:/tmp/.X11-unix:rw \
  -v $PWD/..:/home/dev/workspace/CPlantBox/tutorial/chapter7_coupled \
  -w /home/dev/workspace/CPlantBox/tutorial/chapter7_coupled \
  murilodsv/sp-cpb bash
```

Now execute the `example7_5_coupling_pbcm.py` script:

```bash
python example7_5_coupling_pbcm.py 
```

Output files (CSV and VTP) are saved in:

```bash
ls ~/workspace/simplace_run/output/PyPlantBox/lintul/
```

# Troubleshooting

### Image could not be pulled from [DockerHub](https://hub.docker.com/)

Try building locally. Navigate to the ```Dockerfile``` path and execute:

```bash
docker build -t sp-cpb .
```

### Running the docker without display device

Its normally possible to run the docker image without display device, but the plots would ***NOT*** be forwarded through X11 to your display device.

```bash
docker run -it sp-cpb
```
