import os, sys

BASEDIR = 'graphic/rendered-program-output'
DEST = 'graphic/video'
files = os.listdir(BASEDIR)
L = len(files)
D = len(str(L))

for f in files:
	try:
		newname = str(int(f.split('.')[1])).zfill(D)+'.png'
		os.system(f'mv "{BASEDIR}/{f}" "{BASEDIR}/{newname}"')
	except:
		for _ in range(100): 
			print('Fatal error processing names... directory {BASEDIR} must contain trash')
			print(f, D)
		sys.exit(1)


os.system(f'ffmpeg -framerate 20 -i "{BASEDIR}/%0{D}d.png" -vcodec mpeg4 "{DEST}/result.avi"')
