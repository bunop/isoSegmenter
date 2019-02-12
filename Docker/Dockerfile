
#
# VERSION 0.6
# DOCKER-VERSION  1.9.1
# AUTHOR:         Paolo Cozzi <paolo.cozzi@ibba.cnr.it>
# DESCRIPTION:    A python based image with isoSegmenter installed
# TO_BUILD:       docker build --rm -t bunop/isosegmenter .
# TO_RUN:         docker run -ti bunop/isosegmenter /bin/bash
# TO_TAG:         docker tag bunop/isosegmenter:latest bunop/isosegmenter:0.5
#

FROM python:2.7

LABEL maintainer="Paolo Cozzi <paolo.cozzi@ibba.cnr.it>" \
      description="A python based image with isoSegmenter installed" \
      version="0.6"

# Set corrent working directory
WORKDIR /root

# Install deb dependancies. Then clean  packages and histories
RUN apt-get update \
 && apt-get install -y libgd-dev \
    libgif-dev \
 && apt-get clean

# now clone isoSegmenter project
RUN git clone https://github.com/bunop/isoSegmenter.git /root/isoSegmenter

# Setting working directory
WORKDIR /root/isoSegmenter

# setting isosegmenter version
ENV VERSION=16d4206

# checking out a tagged release
RUN git fetch && \
  git checkout $VERSION && \
  git checkout -b $VERSION

# Byte compiling libraries
RUN pip install .

# Adding a volume
RUN mkdir /data
VOLUME /data

# Setting /data as working directory
WORKDIR /data/

# executing a default command
CMD [ "python", "/usr/local/bin/isoSegmenter.py", "--help" ]
