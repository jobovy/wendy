FROM continuumio/miniconda3

MAINTAINER Jo Bovy

# Add ffmpeg dependency for notebook movies
USER main
RUN wget https://johnvansickle.com/ffmpeg/builds/ffmpeg-git-amd64-static.tar.xz -O ffmpeg-git-amd64-static.tar.xz
USER root
RUN FFMPEG_DIR=$(tar -tJf ffmpeg-git-amd64-static.tar.xz | egrep '^[^/]+/?$') && tar xvJf ffmpeg-git-amd64-static.tar.xz && ln $FFMPEG_DIR/ffmpeg /usr/local/bin/ffmpeg

USER main

ADD requirements.txt requirements.txt
RUN conda install pip
RUN pip install -r requirements.txt
RUN pip install -U --no-deps git+git://github.com/jobovy/wendy.git#egg=wendy
