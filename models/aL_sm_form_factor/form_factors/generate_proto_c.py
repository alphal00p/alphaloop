import re

file_names = ["astu", "atsu", "aust", "bstu", "btsu", "bust", "cstu"]

form_template = """ 
Symbols s,t,u,xs,xu,xt,MT;
CFunctions A02, B02, C02, D02;

Symbols a02, b02s, b02t, b02u, c02s, c02t, c02u, d02su, d02st, d02tu, d02us, d02ts, d02ut, b020;
ExtraSymbols, array, Z;
Format O2;
"""

for name in file_names:
    fortran_expr = open(name + ".dat", "r").read()
    form_expr = fortran_expr.replace("**","^") + ";"
    form_expr = form_expr.replace("+", "\n\t+")
    form_expr = form_expr.replace("-", "\n\t-")
    form_template += "\nLocal " + name + " = " + form_expr + "\n"

form_template += """ 
id A02(MT?) = a02;
id B02(s, ?MT) = b02s;
id B02(t, ?MT) = b02t;
id B02(u, ?MT) = b02u;
id C02(0,0,s,?MT) = c02s;
id C02(0,0,t,?MT) = c02t;
id C02(0,0,u,?MT) = c02u;
id D02(0,0,0,0,s,u,?MT) = d02su;
id D02(0,0,0,0,s,t,?MT) = d02st;
id D02(0,0,0,0,t,u,?MT) = d02tu;
id D02(0,0,0,0,u,s,?MT) = d02us;
id D02(0,0,0,0,t,s,?MT) = d02ts;
id D02(0,0,0,0,u,t,?MT) = d02ut;
id B02(0,?MT) = b020;
id d02us = d02su;
id d02ts = d02st;
id d02ut = d02tu;

Format C;
.sort
"""

for name in file_names:
    form_template += """
#optimize {0}
#write <{0}.proto_c> \"NUM_VAR=`optimmaxvar_' \"
#write <{0}.proto_c> \"%O\"
#write <{0}.proto_c> \"{0} = %e \", {0}
#clearoptimize
.sort
""".format(name)

import os 
os.remove("form_factors.frm")
form_file = open("form_factors.frm", "a")
form_file.write(form_template)
form_file.close()
os.system("form form_factors.frm") 


