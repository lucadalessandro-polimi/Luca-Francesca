FROM aldoclemente/fdapde-docker

# install miniconda
RUN curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
 && bash /tmp/miniconda.sh -b -p /opt/conda \
 && rm /tmp/miniconda.sh

ENV PATH=/opt/conda/bin:$PATH

# install nlohmann/json single-header
RUN mkdir -p /usr/include/external/nlohmann && \
    curl -L https://raw.githubusercontent.com/nlohmann/json/develop/single_include/nlohmann/json.hpp \
    -o /usr/include/external/nlohmann/json.hpp

# accetta ToS
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main \
 && conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# crea env py312
RUN conda create -y -n py312 python=3.12 numpy matplotlib pandas seaborn

WORKDIR /root/progetto
