
#
# VERSION 0.3
# DOCKER-VERSION  1.9.1
# AUTHOR:         Paolo Cozzi <paolo.cozzi@ptp.it>
# DESCRIPTION:    An image python base image which isoSegmenter will run
# TO_BUILD:       docker build --rm -t bunop/isosegmenter .
# TO_RUN:         docker run -ti bunop/isosegmenter /bin/bash
# TO_TAG:         docker tag bunop/isosegmenter:latest bunop/isosegmenter:0.3
#

FROM python:2.7.10

MAINTAINER Paolo Cozzi <paolo.cozzi@tecnoparco.org>

# Set corrent working directory
WORKDIR /root

# Install deb dependancies
RUN apt-get update && apt-get install -y libgd-dev \
    libgif-dev

# Clean packages and histories
RUN apt-get clean && rm -rf /var/lib/apt/lists/

# Install gdmodule manually
RUN git clone https://github.com/Solomoriah/gdmodule.git

# Set corrent working directory
WORKDIR /root/gdmodule

# Install package from github repository
RUN pip install .

# now clone isoSegmenter project
RUN git clone https://github.com/bunop/isoSegmenter.git /root/isoSegmenter

# Setting working directory
WORKDIR /root/isoSegmenter

# Byte compiling libraries
RUN pip install .

# Adding a volume
RUN mkdir /data
VOLUME /data

# Setting /data as working directory
WORKDIR /data/

# executing a default command
CMD [ "python", "/usr/local/bin/isoSegmenter.py", "--help" ]
