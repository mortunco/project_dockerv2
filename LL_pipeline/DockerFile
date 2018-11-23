FROM ubuntu:16.04
MAINTAINER Tunc Morova, tmorova15@ku.edu.tr

RUN \
  sed -i 's/# \(.*multiverse$\)/\1/g' /etc/apt/sources.list && \
  apt-get update && \
  apt-get -y upgrade && \
  apt-get install -y build-essential && \
  apt-get install -y software-properties-common && \
  apt-get install -y curl man unzip vim wget && \
  apt-get -y install python2.7 python-pip && \
  apt-get install -y aptitude && \
  apt-get install -y libcurl4-openssl-dev && \
  apt-get install -y libxml2-dev && \
  apt-get install -y libssl-dev && \
  apt-get install -y libmariadb-client-lgpl-dev && \
  rm -rf /var/lib/apt/lists/*

### Install java 7
RUN add-apt-repository ppa:openjdk-r/ppa && \
    apt-get update && \
    apt-get install -y openjdk-7-jdk && \
    apt-get install -y ant && \
    apt-get clean;


### Install snpEFF and snpSIFT
RUN wget https://sourceforge.net/projects/snpeff/files/snpEff_v4_2_core.zip
RUN unzip snpEff_v4_2_core.zip

RUN apt-get update && apt-get install -y bedtools apt-transport-https ca-certificates


RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu xenial/" >>/etc/apt/sources.list
#RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
RUN apt-get update &&  apt-get install -y --allow-unauthenticated  r-base r-base-dev

ADD ./myscript.R /usr/myscript.R
RUN Rscript /usr/myscript.R
RUN echo "Java Check..." && java -version
RUN echo "Python Check" && python --version && pip --version
RUN echo "R check" && R --version
WORKDIR /project/Batch1/
#CMD sh firestarter.sh