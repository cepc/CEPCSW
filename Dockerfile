##############################################################################
# Docker image for CEPCSW
##############################################################################

# # Instruction
# To build the docker image:
#   $ docker build -t cepc/cepcsw .
# Or with CVMFS installed also
#   $ docker build -t cepc/cepcsw-cvmfs . --build-arg CVMFSMOD=INSIDE
#
# To publish it to DockerHub:
#   $ docker push cepc/cepcsw

FROM centos:7

ARG CVMFSMOD

# Basic
RUN yum install -y sudo
RUN sudo yum install -y redhat-lsb wget

# Install singularity
RUN sudo yum update -y && \
    sudo yum install -y epel-release && \
    sudo yum update -y && \
    sudo yum install -y singularity

# If the CVMFS is installed in the host, just mount the corresonding /CVMFS directories
# $ docker run --privileged --rm -i -t \
#           -v /cvmfs/sft.cern.ch:/cvmfs/sft.cern.ch \
#           -v /cvmfs/cepcsw.ihep.ac.cn:/cvmfs/cepcsw.ihep.ac.cn \
#           -v /cvmfs/container.ihep.ac.cn:/cvmfs/container.ihep.ac.cn \
#           cepc/cepcsw /bin/bash
# Inside the Docker container, we could start the Singularity provided by IHEP
# $ export SINGULARITY_BINDPATH=/cvmfs
# $ singularity shell /cvmfs/container.ihep.ac.cn/singularity/image/SL69/sl69worknode20200729.sif

##############################################################################
# Install CVMFS
# Configure IHEP
# Enable /cvmfs/cepcsw.ihep.ac.cn and /cvmfs/container.ihep.ac.cn

RUN if [ "$CVMFSMOD" = "INSIDE" ]; then \
  sudo yum install -y https://ecsft.cern.ch/dist/cvmfs/cvmfs-release/cvmfs-release-latest.noarch.rpm \
  && sudo yum install -y cvmfs \
  && sudo mkdir /etc/cvmfs/keys/ihep.ac.cn \
  && sudo curl -o /etc/cvmfs/keys/ihep.ac.cn/ihep.ac.cn.pub http://cvmfs-stratum-one.ihep.ac.cn/cvmfs/software/client_configure/ihep.ac.cn/ihep.ac.cn.pub \
  && sudo curl -o /etc/cvmfs/domain.d/ihep.ac.cn.conf http://cvmfs-stratum-one.ihep.ac.cn/cvmfs/software/client_configure/ihep.ac.cn.conf \
  && echo "CVMFS_REPOSITORIES='sft.cern.ch,cepcsw.ihep.ac.cn,container.ihep.ac.cn'" | sudo tee    /etc/cvmfs/default.local \
  && echo "CVMFS_HTTP_PROXY=DIRECT"                                                 | sudo tee -a /etc/cvmfs/default.local \
  && cat /etc/cvmfs/default.local \
  && sudo mkdir -p /cvmfs/sft.cern.ch \
  && sudo mkdir -p /cvmfs/cepcsw.ihep.ac.cn \
  && sudo mkdir -p /cvmfs/container.ihep.ac.cn; \
  fi

# START Container:
# # docker run  --privileged  --rm -i -t cepc/cepcsw-cvmfs /bin/bash
# Due to the fuse issue, following commands need to be run inside container when --privileged is specified
# $ mount -t cvmfs sft.cern.ch /cvmfs/sft.cern.ch
# $ mount -t cvmfs container.ihep.ac.cn /cvmfs/container.ihep.ac.cn
# $ mount -t cvmfs cepcsw.ihep.ac.cn /cvmfs/cepcsw.ihep.ac.cn

##############################################################################
# Install necessary packages
##############################################################################

RUN yum install -y git
RUN yum install -y libglvnd-devel
RUN yum install -y mesa-libGLU-devel
RUN yum install -y libXmu-devel
RUN yum install -y motif-devel

# For runtime
RUN yum install -y compat-db47
