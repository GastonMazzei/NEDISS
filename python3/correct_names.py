import os

BASEDIR = 'graphic/rendered-program-output'
DEST = 'graphic/video'
files = os.listdir(BASEDIR)
L = len(files)
D = len(str(L))

for f in files:
	newname = str(int(f.split('.')[1])).zfill(D)+'.png'
	os.system(f'mv "{BASEDIR}/{f}" "{BASEDIR}/{newname}"')



os.system(f'ffmpeg  -i "{BASEDIR}/%0{D}d.png" -vcodec mpeg4 "{DEST}/result.avi"')
