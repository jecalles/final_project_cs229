'''
context.py: a file that automagically adds the parent file to the absolute path. Makes importing file easier.

Ex:
---
# /scripts/a_file.py
from context import src, res
'''
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
