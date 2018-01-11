#!/usr/bin/python

'''Replaces spaces in filenames with underscores. Thanks to http://littlebrain.org/2008/10/19/replace-space-with-undescore-in-filename/.'''
 
import os
import sys
 
files = os.listdir(sys.argv[1])
for f in files:
	if f[0] == '.':
		print('Skipping '+f)
		continue
	if ' ' in f:
		skip = raw_input('Replacing '+f+'. Hit enter to continue or type s to skip:  ')
		if not skip:
			print(f+' -> '+f.replace(' ', '_'))
			os.rename(f, f.replace(' ', '_'))
		else:
			print('Skipping '+f)
			
			
