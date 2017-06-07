FROM andrewosh/binder-base

MAINTAINER Jo Bovy

USER root

# Add ffmpeg dependency for notebook movies
RUN apt-get update
RUN apt-get install -y ffmpeg

USER main

# Install requirements for Python 2 and 3
ADD requirements.txt requirements.txt
RUN pip install -r requirements.txt
RUN /home/main/anaconda/envs/python3/bin/pip install -r requirements.txt
RUN pip install -U --no-deps git+git://github.com/jobovy/wendy.git#egg=wendy
RUN /home/main/anaconda/envs/python3/bin/pip install -U --no-deps git+git://github.com/jobovy/wendy.git#egg=wendy
