#!/usr/bin/env python3
import sys
mode=sys.argv[1]

tt=open(sys.argv[2],'r').read()
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
else:
    for i_line, line in enumerate(tt.split('\n')):
        if i_line==0 and line[-2:]!='{}':
            new_lines.append('%s {}'%line)
            continue
        if i_line==3 and line[-3:]!=' {}':
            new_lines.append(line[:-3])
            continue
        if i_line in [1,2] and line[0]!='#':
            new_lines.append('#%s'%line)
            continue
        if i_line in [4,5,6,7] and line[0]=='#':
            new_lines.append(line[1:])
            continue
        new_lines.append(line)
open(sys.argv[2],'w').write('\n'.join(new_lines))
