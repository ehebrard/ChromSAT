FROM ubuntu:17.10

RUN requirements="make git hgsubversion ninja-build python3 python3-pip python3-dev gcc-7 clang-5.0 openssh-client" && \
    apt-get update && \
    apt-get install --no-install-recommends -y $requirements

RUN pip3 install setuptools wheel  # C'est un peu stupide mais c'est pas par défault...
RUN pip3 install meson

WORKDIR /tmp/build/

RUN hg clone http://hg@bitbucket.org/gkatsi/minicsp

ADD ./tclap /tmp/build/minicsp/tclap
#
ADD ./src /tmp/build/minicsp/src
ADD ./meson.build /tmp/build/minicsp

RUN ls -all minicsp


RUN CXX=clang++-5.0 meson  --buildtype=debugoptimized minicsp minicsp/debugopt
RUN ninja -C minicsp/debugopt/

RUN ln -s /tmp/build/minicsp/debugopt/gc /usr/bin