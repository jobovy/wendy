FROM andrewosh/binder-base

MAINTAINER Jo Bovy

# Add ffmpeg dependency for notebook movies
USER main
RUN wget https://johnvansickle.com/ffmpeg/builds/ffmpeg-git-amd64-static.tar.xz -O ffmpeg-git-amd64-static.tar.xz
RUN file ffmpeg-git-amd64-static.tar.xz
RUN FFMPEG_DIR=$(tar -ztf ffmpeg-git-amd64-static.tar.xz | egrep '^[^/]+/?$') && echo $FFMPEG_DIR
USER root
RUN FFMPEG_DIR=$(tar -ztf ffmpeg-git-amd64-static.tar.xz | egrep '^[^/]+/?$') && tar xvf ffmpeg-git-amd64-static.tar.xz && ln $FFMPEG_DIR/ffmpeg /usr/local/bin/ffmpeg

USER main

# Install requirements for Python 2 and 3
ADD requirements.txt requirements.txt
RUN pip install -r requirements.txt
RUN /home/main/anaconda/envs/python3/bin/pip install -r requirements.txt
RUN pip install -U --no-deps git+git://github.com/jobovy/wendy.git#egg=wendy
RUN /home/main/anaconda/envs/python3/bin/pip install -U --no-deps git+git://github.com/jobovy/wendy.git#egg=wendy
