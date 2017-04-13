# -*- makefile -*-
# this makefile contains rules for manaing TSIG code
# Copyright 2017 MIT
# Licensed under terms of the GPLv3

all: help

help:
	@echo "options include:"
	@echo "    test  run the unit tests and examples"

test:
	./setup.py test

clean:
	find . -name "*.pyc" -exec rm {} \;
	rm -rf *.egg-info
