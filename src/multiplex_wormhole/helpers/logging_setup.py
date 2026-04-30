#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: LOGGING SETUP
Sets of log files to run within functions
Created on Wed Apr 29 16:28:40 2026
@author: maggiehallerud
"""

import sys
import logging


def setup_logging(LOGFILE, VERBOSE, NAME="a"):
    # initialize new logger
    logger = logging.getLogger(NAME)
    # remove & close existing handlers
    for handler in logger.handlers[:]:
        handler.close()
        logger.removeHandler(handler)
    logger.propagate=False
    # configure
    level = logging.DEBUG if VERBOSE else logging.INFO
    logger.setLevel(level)
    formatter = logging.Formatter('%(message)s')
    # set up file handler
    file_handler = logging.FileHandler(LOGFILE, encoding='utf-8')
    file_handler.setLevel(level)
    file_handler.setFormatter(formatter)
    # set up console handler
    console_handler = logging.StreamHandler(sys.__stdout__)
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)
    # add both handlers
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    # clear handlers
    for h in logger.handlers:
        h.flush()
    return logger
