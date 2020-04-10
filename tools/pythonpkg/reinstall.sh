#!/bin/bash
rm -rf build  dist  duckdb.cpp  duckdb.egg-info  duckdb.hpp
python3 setup.py install --user
