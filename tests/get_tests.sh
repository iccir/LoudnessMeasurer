#!/bin/sh

curl https://tech.ebu.ch/docs/testmaterial/ebu-loudness-test-setv03.zip > tests.zip
unzip tests.zip
rm -rf tests.zip
