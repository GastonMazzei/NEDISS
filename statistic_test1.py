
with open('tmp/trash.txt','r') as f: d = f.readlines()
 
cs = {x:0 for x in range(11)}

for _ in d:
     for k in cs.keys():
             #if f'Added ix {k} to CHECKED' in _: cs[k]+=1;
             if f'starting to locate this ix: {k}' in _: cs[k]+=1;
 
for k in cs.keys():
	print(k,cs[k])
