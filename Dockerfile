# take z3 image as base
FROM ghcr.io/z3prover/z3:ubuntu-20.04-bare-z3-sha-d66609e

# Install git
RUN apt-get update && apt-get install -y git

# Install python3.8 and pip
RUN apt-get install -y build-essential python3.8 python3-pip python3-dev

RUN pip3 -q install pip

RUN mkdir MCC

WORKDIR /MCC/

RUN git clone https://github.com/Biomathsys/MassChargeCuration.git
RUN pip3 install ./MassChargeCuration

RUN pip3 install cobra
RUN pip3 install jupyter

# set up visual studio code extensions
RUN curl -fsSL https://code-server.dev/install.sh | sh
RUN code-server --install-extension ms-python.python
RUN code-server --install-extension ms-toolsai.jupyter

# overwrite z3 default entrypoint
ENTRYPOINT ["/bin/bash"]
