FROM ubuntu:16.04
MAINTAINER Kevin Murray <spam@kdmurray.id.au>

ADD . /experiments
RUN bash /experiments/docker/build.sh
WORKDIR /experiments
