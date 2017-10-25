#!/usr/bin/env python2
import os
import shutil as sh

dirlist = os.listdir('./')

deletelist = []
keeplist = []

for f in dirlist:

    if (f == '0'):
        continue

    # Delete if filename can be turned into a float and is not 0
    try:

        testfloat = float(f)
        deletelist.append(f)

    except ValueError:

        if ('processor' in f):
            deletelist.append(f)

keeplist = [i for i in dirlist if i not in deletelist]


print('Keeping the following files:')
print(keeplist)

print('Deleting the following files: ')
print(deletelist)

answer = raw_input('Proceed? [y]/n: ')
if (answer == 'y' or answer == 'Y'):
    for f in deletelist:
        sh.rmtree(f)
else:
    quit('Cancelled deletion.')

def listdir_fullpath(d):
        return [os.path.join(d, f) for f in os.listdir(d)]

keep_files = ["0/p", "0/U", "constant/polyMesh/blockMeshDict"]
filelist = listdir_fullpath('0') +  listdir_fullpath('constant/polyMesh')

for f in  [x for x in filelist if x not in keep_files]:
    os.remove(f)
