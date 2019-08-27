FROM continuumio/miniconda3

MAINTAINER Jo Bovy

### create user with a home directory
ARG NB_USER
ARG NB_UID
ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}
WORKDIR ${HOME}

# Add ffmpeg dependency for notebook movies
USER ${NB_UID}
RUN wget https://johnvansickle.com/ffmpeg/builds/ffmpeg-git-amd64-static.tar.xz -O ffmpeg-git-amd64-static.tar.xz
USER root
RUN FFMPEG_DIR=$(tar -tJf ffmpeg-git-amd64-static.tar.xz | egrep '^[^/]+/?$') && tar xvJf ffmpeg-git-amd64-static.tar.xz && ln $FFMPEG_DIR/ffmpeg /usr/local/bin/ffmpeg

USER ${NB_UID}

ADD requirements.txt requirements.txt
RUN conda install pip
RUN pip install -r requirements.txt
RUN pip install -U --no-deps git+git://github.com/jobovy/wendy.git#egg=wendy
