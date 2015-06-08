
#
# VERSION 0.1
# DOCKER-VERSION  1.6.2
# AUTHOR:         Paolo Cozzi <paolo.cozzi@tecnoparco.org>
# DESCRIPTION:    An image python base image which isochoreFinder will run
# TO_BUILD:       docker build --rm -t bunop/isochorefinder .
# TO_RUN:         docker run -ti bunop/isochorefinder /bin/bash
# TO_TAG:         docker tag bunop/isochorefinder:latest bunop/isochorefinder:0.1
#

FROM python:2.7.9

MAINTAINER Paolo Cozzi <paolo.cozzi@tecnoparco.org>

# Set corrent working directory
WORKDIR /root

# Install deb dependancies
RUN apt-get update && apt-get install -y libgd2-xpm-dev \
    libgif-dev

# Clean packages and histories
RUN apt-get clean && rm -rf /var/lib/apt/lists/

# Install isochoreFinder requirements
RUN pip install Pillow matplotlib biopython

# Install gdmodule
ADD https://pypi.python.org/packages/source/g/gdmodule/gdmodule-0.58.tar.gz /root/gdmodule-0.58.tar.gz
RUN tar -xzf /root/gdmodule-0.58.tar.gz

# Set corrent working directory
WORKDIR /root/gdmodule-0.58

# Install package manually (note the invocation of Setup - this seems to be the principal reason
# that we can't install gd libraries using pip)
RUN python Setup.py build && python Setup.py install

# now clone isochoreFinder in /opt directory
RUN git clone https://github.com/bunop/isochoreFinder.git /opt/isochoreFinder

# Setting the final working directory
WORKDIR /opt/isochoreFinder

# Byte compiling libraries
RUN python -m compileall -f GClib

# simbolick linking
RUN ln -s /opt/isochoreFinder/isochoreFinder.py /usr/local/bin && \
    ln -s /opt/isochoreFinder/isochoreFamilies.py /usr/local/bin

# Adding a volume
RUN mkdir /data
VOLUME /data

# executing a default command
CMD [ "python", "/opt/isochoreFinder/isochoreFinder.py", "--help" ]
