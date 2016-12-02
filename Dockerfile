FROM apiaryio/emcc:1.36


RUN apt-get install wget -y
WORKDIR /lib/
RUN wget -q -S -O - 'http://downloads.sourceforge.net/project/boost/boost/1.61.0/boost_1_61_0.tar.gz' | tar xz


RUN apt-get install mercurial -y
RUN hg clone https://bitbucket.org/eigen/eigen
WORKDIR /src