

with open('graphviz_test_new.dot','r') as f:
	d = f.readlines()

with open('graphviz_test_new.edited.dot','w') as f:
	for l in d:
		if '->' in l:
			f.write(l[:-2] + ' [dir=none];\n')
		else:
			f.write(l)
		
