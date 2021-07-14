
isoSegmenter
==============

A program for segmenting genomes into isochores

[![Nosetest workflow](https://github.com/bunop/isoSegmenter/actions/workflows/main.yml/badge.svg)](https://github.com/bunop/isoSegmenter/actions/workflows/main.yml)
[![Coverall](https://coveralls.io/repos/github/bunop/isoSegmenter/badge.svg?branch=master)](https://coveralls.io/github/bunop/isoSegmenter?branch=master)
![Docker Automated build](https://img.shields.io/docker/automated/bunop/isosegmenter)

Installing isoSegmenter using Docker
------------------------------------

Those are instructions on how tu run [isoSegmenter](https://github.com/bunop/isoSegmenter) inside a docker container. By using docker, you will start a container in which isoSegmenter is already installed (you don't need to install all requirements since they are already satisfied in image). First install [docker](http://docs.docker.com/installation/#installation) on your platform. Then you can download the pre-build docker image or modify and build a new docker image

## Running isoSegmenter under a pre builded image

You can start a isosegmenter image simply by typing:

```bash
docker run -ti bunop/isosegmenter /bin/bash
```

This will be enough to download and get an isoSegmenter running container. You can otionally mount a local directoy inside the running container by adding a local directory as a docker volume. More information can be found inside Docker tutorial on [Managing data in Container](http://docs.docker.com/userguide/dockervolumes/). `isoSegmenter.py` will be placed under `/usr/local/bin/` directory

## Build you own image

Clone isoSegmenter repository from github on your machine, enter in Docker directory and build the docker images. You can modify the image to satisfy your needs. Then you have to run the builded container with the same name you have defined previously.

```bash
$ git clone https://github.com/bunop/isoSegmenter.git
$ cd Docker
$ docker build --rm -t bunop/isosegmenter .
$ docker run -ti bunop/isosegmenter /bin/bash
```

More information can be found on [Docker documentation](http://docs.docker.com/) or on [Working with Docker Images](http://docs.docker.com/userguide/dockerimages/) tutorial
