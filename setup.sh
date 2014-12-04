#!/bin/bash

python setup.py build
python setup.py clean
mv -v build/lib*/*.so utilities/
rm -rf build