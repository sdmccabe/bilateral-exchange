FROM ubuntu:14.04

RUN apt-get clean
RUN rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN echo "deb http://dl.bintray.com/sbt/debian /" | tee -a /etc/apt/sources.list.d/sbt.list
RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv 642AC823
RUN apt-get update
RUN apt-get install -yqq time git curl build-essential checkinstall autotools-dev wget cmake software-properties-common
RUN apt-get clean
RUN rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN mkdir model



COPY RNG.cpp model/
COPY RNG.h model/
COPY main.cpp model/
COPY main.h model/
COPY parameters.cfg model/
COPY easylogging++.h model/
COPY Makefile model/

WORKDIR /model/

RUN wget http://www.hyperrealm.com/libconfig/libconfig-1.5.tar.gz
RUN tar xvf libconfig-1.5.tar.gz
RUN wget https://github.com/gflags/gflags/archive/v2.1.2.tar.gz
RUN tar xvf v2.1.2.tar.gz
WORKDIR /model/libconfig-1.5
RUN ./configure && make && make install
WORKDIR /model/gflags-2.1.2/build
RUN cmake .. && make && make install
RUN ldconfig

WORKDIR /model/
CMD ["bash"]

#CMD /
#CMD go build zi-traders.go
#CMD echo "Built Go model"
#CMD bash -c "/usr/bin/time -f '1,%e,%U,%S' ./zi-traders"

#docker-machine create -d virtualbox --virtualbox-boot2docker-url file://$HOME/Dropbox/boot2docker-v1.9.1-fix1.iso --virtualbox-memory 1536 --virtualbox-disk-size 10000 fixedjava
#docker-machine create -d virtualbox --virtualbox-boot2docker-url https://github.com/tianon/boot2docker-legacy/releases/download/v1.10.0-rc1/boot2docker.iso --virtualbox-memory 1536 --virtualbox-disk-size 10000 --virtualbox-cpu-count 2 fixedjava
