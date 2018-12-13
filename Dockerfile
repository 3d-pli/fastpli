# docker build --add-host dl-cdn.alpinelinux.org:151.101.112.249 -t fastpli .
# docker run fastpli

FROM alpine

RUN apk update
RUN apk add make
RUN apk add gcc
RUN apk add g++
RUN apk add cmake
RUN apk add git
RUN apk add python3
RUN apk add python3-dev
RUN apk add py3-virtualenv

ADD . /build/f.matuschke
COPY . /build/f.matuschke/fastpli
WORKDIR /build/f.matuschke/fastpli

CMD make build