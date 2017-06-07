FROM andrewosh/binder-base

MAINTAINER Jo Bovy

# Add ffmpeg dependency for notebook movies
USER main
RUN wget http://johnvansickle.com/ffmpeg/releases/ffmpeg-release-64bit-static.tar.xz
RUN tar xvfJ ffmpeg-release-64bit-static.tar.xz
USER root
RUN ln ffmpeg-3.1.2-64bit-static/ffmpeg /usr/local/bin/ffmpeg

USER main

# Install requirements for Python 2 and 3
ADD requirements.txt requirements.txt
RUN pip install -r requirements.txt
RUN /home/main/anaconda/envs/python3/bin/pip install -r requirements.txt
RUN pip install -U --no-deps git+git://github.com/jobovy/wendy.git#egg=wendy
RUN /home/main/anaconda/envs/python3/bin/pip install -U --no-deps git+git://github.com/jobovy/wendy.git#egg=wendy
