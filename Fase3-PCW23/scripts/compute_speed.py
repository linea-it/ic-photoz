
import datetime 
import time 


name = ['50 files', '100 files', '150 files', '314 files']
obj = [88895872, 178386176, 266835897, 556636904] 
fzboost = ['0:11:07', '0:24:47', '0:30:26', '0:57:25'] 
bpz = ['0:07:52', '0:14:48', '0:26:24', '0:42:52']

print()
print('FlexZBoost')
for i in range(4):
    x = time.strptime(fzboost[i],'%H:%M:%S')
    sec = datetime.timedelta(hours=x.tm_hour,minutes=x.tm_min,seconds=x.tm_sec).total_seconds()
    speed = (sec*1000)/obj[i]
    print(f'{name[i]}: {round(speed,3)}ms/obj')

print()
print('BPZ')
for i in range(4):
    x = time.strptime(bpz[i],'%H:%M:%S')
    sec = datetime.timedelta(hours=x.tm_hour,minutes=x.tm_min,seconds=x.tm_sec).total_seconds()
    speed = (sec*1000)/obj[i]
    print(f'{name[i]}: {round(speed,3)}ms/obj')

