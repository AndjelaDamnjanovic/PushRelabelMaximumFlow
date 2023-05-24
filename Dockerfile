FROM ubuntu:latest

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y 
RUN apt-get install qtcreator -y
RUN apt-get install build-essential -y 
RUN apt install qtbase5-dev qt5-qmake -y
RUN apt-get install xcb

COPY . /usr/src/graph 
WORKDIR /usr/src/graph

RUN mkdir app
WORKDIR /usr/src/graph/app

RUN qmake ../GraphVisual
RUN make

ENTRYPOINT ["./graph"]