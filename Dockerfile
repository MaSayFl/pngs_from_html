FROM python:3.9-slim-bullseye

LABEL base.image="python:3.9-slim-bullseye"
LABEL description="image for test"
LABEL build.date=$builddate

RUN apt-get update && apt-get -y upgrade && apt-get clean && rm -rf /var/lib/apt/lists/*
WORKDIR /usr/src/app
COPY . .