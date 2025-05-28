import os
import subprocess

'''
This script will execute a python script inside subdirectories. 
This assumes the file to be executed is within those sub-directories.

#####Replace directory_to_check with "/PATH/TO/PARENT_DIR"
'''

directory_to_check = ""

# Get all the subdirectories of directory_to_check recursively and store them in a list:
directories = [os.path.abspath(x[0]) for x in os.walk(directory_to_check)]
directories.remove(os.path.abspath(directory_to_check)) # If you don't want your main directory included

for i in directories:
      os.chdir(i)
      subprocess.run(["python", "SCRIPT_TO_RUN.py"])
