1: In order to build the form factors, first compile OneLOop (I recommend 
commenting out the part that prints the banner first, so it doesn't break 
the rust dashboard).

2: copy libavh_olo.a to this folder. 

3: run the mathematica notebook OneLoopFormFactor4Mathijs.nb 

4: run generate_proto_c.py 

5: run process_proto_c.py