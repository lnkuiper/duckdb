#!/bin/bash
rm -rf build  dist  duckdb.cpp  duckdb.egg-info  duckdb.hpp
python setup.py install --user
