FROM alpine

RUN apk update
RUN apk add gcc g++ make cmake git
RUN apk add python3-dev py3-virtualenv

ADD . /code/
COPY . /code/fastpli
WORKDIR /code/fastpli
RUN git clean -d -f -x

CMD make build