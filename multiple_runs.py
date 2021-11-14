import os, sys


command = sys.argv[1]


errors = 0
success = 0

while (errors == 0):
    response = os.system(command + " >  tmp/trash.txt")
    if response == 0: success += 1
    else: errors += 1
    try: os.remove('tmp/trash.txt')
    except: pass
    print(f'So far no errors after {success} runs ;-)')

print(f'\n--------------\nFinal report: observed {errors} errors in {success+errors} runs\n------------------\n')
