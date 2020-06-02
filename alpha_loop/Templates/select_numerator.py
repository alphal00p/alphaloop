#!/usr/bin/env python3
import sys
import glob
mode=sys.argv[1]

all_yaml_files = glob.glob('*.yaml')
all_yaml_files_filtered = []
for p in all_yaml_files:
    try:
        int(p.split('_')[-1].split('.')[0])
        all_yaml_files_filtered.append(p)
    except:
        continue

for p in all_yaml_files_filtered:
 print('Processing %s ...'%p)
 tt=open(p,'r').read()
 new_lines = []
 if mode.upper().startswith('F'):
    for i_line, line in enumerate(tt.split('\n')):
        if i_line==0 and line[-3:]==' {}':
            new_lines.append(line[:-3])
            continue
        if i_line==3 and line[-2:]!='{}':
            new_lines.append('%s {}'%line)
            continue
        if i_line in [1,2] and line[0]=='#':
            new_lines.append(line[1:])
            continue
        if i_line in [4,5,6,7] and line[0]!='#':
            new_lines.append('#%s'%line)
            continue
        new_lines.append(line)
 elif mode.upper().startswith('M'):
    for i_line, line in enumerate(tt.split('\n')):
        if i_line==0 and line[-2:]!='{}':
            new_lines.append('%s {}'%line)
            continue
        if i_line==3 and line[-3:]==' {}':
            new_lines.append(line[:-3])
            continue
        if i_line in [1,2] and line[0]!='#':
            new_lines.append('#%s'%line)
            continue
        if i_line in [4,5,6,7] and line[0]=='#':
            new_lines.append(line[1:])
            continue
        new_lines.append(line)
 else:
    for i_line, line in enumerate(tt.split('\n')):
        if i_line==0 and line[-3:]==' {}':
            new_lines.append(line[:-3])
            continue
        if i_line==3 and line[-3:]==' {}':
            new_lines.append(line[:-3])
            continue
        if i_line in [1,2] and line[0]=='#':
            new_lines.append(line[1:])
            continue
        if i_line in [4,5,6,7] and line[0]=='#':
            new_lines.append(line[1:])
            continue
        new_lines.append(line)

 open(p,'w').write('\n'.join(new_lines))
