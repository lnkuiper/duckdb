#!/bin/bash
rm -rf build  dist  duckdb.cpp  duckdb.egg-info  duckdb.hpp
rm -rf ~/.local/lib/python3.7/site-packages/duckdb*
python3 setup.py install --user
