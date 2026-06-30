FROM wmoore28/geant4:10.5.0_with_source

# --- agent-proxy + CentOS-7-EOL accommodations -----------------------------
# CentOS 7 is EOL: mirrorlists are gone, so point yum at vault.centos.org and
# route through the session HTTPS proxy (TLS re-terminated -> trust its CA).
ARG PROXY=http://127.0.0.1:35179
COPY ca-bundle.crt /etc/pki/ca-trust/source/anchors/agent-proxy.crt
RUN update-ca-trust extract || true
RUN sed -i -e 's/^mirrorlist=/#mirrorlist=/g' \
           -e 's|^#\?baseurl=http://mirror.centos.org|baseurl=https://vault.centos.org|g' \
           /etc/yum.repos.d/CentOS-Base.repo \
 && printf 'proxy=%s\nsslcacert=/etc/pki/ca-trust/source/anchors/agent-proxy.crt\n' "$PROXY" >> /etc/yum.conf

COPY tarballs/G4beamline-3.06-source.tgz /tmp/G4beamline.tgz
COPY test/g4bl-run.sh /share/g4bl-run.sh
COPY test/g4bl-test.g4bl /share/g4bl-test.g4bl

# Dependencies (gsl/gsl-devel already in base; add FFTW).
RUN yum install -q -y fftw-devel && yum clean all
ENV GSL_DIR  "/usr"
ENV FFTW_DIR "/usr"

ENV Geant4_DIR /opt/geant4
ENV G4BL_DIR /opt/G4beamline
WORKDIR $G4BL_DIR
RUN mkdir /opt/G4beamline-source \
    && tar xzf /tmp/G4beamline.tgz -C /opt/G4beamline-source --strip-components 1 \
    && find /opt/G4beamline-source -name CMakelists.txt -exec rename lists.txt Lists.txt '{}' \; \
    && cmake ../G4beamline-source > /dev/null \
    && make -j `nproc` install > /dev/null \
    && sed -i 's/echo/#echo/g' bin/g4bl-setup.sh \
    && echo $GEANT4_DIR/share/*/data > $G4BL_DIR/.data
ENV PATH "$G4BL_DIR/bin:$PATH"

WORKDIR /root
RUN rm -rf /opt/G4beamline-source $G4BL_DIR/*.tgz /tmp/G4beamline.tgz \
    && echo "source $G4BL_DIR/bin/g4bl-setup.sh" >> ~/.bashrc \
    && systemd-machine-id-setup \
    && sed -i 's/\r//g' /share/g4bl-run.sh

CMD [ "/bin/bash" ]
