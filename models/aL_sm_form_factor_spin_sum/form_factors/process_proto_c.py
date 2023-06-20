import re
import os

"""Some global data"""

file_names = ["astu", "atsu", "aust", "bstu", "btsu", "bust", "cstu"]
function_names = ["APHOAMPFFSTU", "APHOAMPFFTSU", "APHOAMPFFUST", "BPHOAMPFFSTU", "BPHOAMPFFTSU", "BPHOAMPFFUST", "CPHOAMPFFSTU"]
common_denominators = ["(96.0*xs*xs*xt*xt*xu*xu*E1*E2*E3*E42*s*s)",
"(96.0*xs*xs*xt*xt*xu*xu*E1*E2*E3*E42*t*t)",
"(96.0*xs*xs*xt*xt*xu*xu*E1*E2*E3*E42*u*u)",
 "(48.0*xs*xs*xs*xt*xt*xt*xu*xu*xu*(E1*E2*E32*E42*s*s*t*t*u))",
 "(48.0*xs*xs*xs*xt*xt*xt*xu*xu*xu*(E12*E2*E3*E42*s*s*t*t*u))",
 "(48.0*xs*xs*xs*xt*xt*xt*xu*xu*xu*(E1*E22*E3*E42*s*s*t*u*u))",
 "(3.0*xs*xs*xs*xs*xt*xt*xt*xt*xu*xu*xu*xu*(E1*E2*E3*E4*s*t*t*u))"]

# Define the replacement function
def replace_fn_array(match):
    n = int(match.group(1))
    return f"Z[{n-1}]"

def replace_fn_pow(match):
    x = match.group(1)
    n = int(match.group(2))
    return '*'.join([x]*n)


template_code = open("form_factors_template.cpp", "r").read()

for i in range(len(file_names)):
    main_code = ""
    file_name = file_names[i]
    function_name = function_names[i]
    common_denominator = common_denominators[i]
    proto_c_file = open(file_name + ".proto_c", "r")
    proto_c = proto_c_file.read()

    find_num_var_pattern = r"NUM_VAR=(\d+)"
    num_var = re.match(find_num_var_pattern, proto_c).group(1)
    proto_c = proto_c.replace("NUM_VAR={0}".format(num_var), "")
    # Define the regular expression pattern
    pattern_array = r"Z\[(\d+)\]"
    pattern_pow = r"pow\((\w+),(\d+)\)"

    # Replace every instance of Z[n] with Z[n-1]
    proto_c= re.sub(pattern_array, replace_fn_array, proto_c)
    # Replace pow(x, n) with x^n using only multiplication
    proto_c= re.sub(pattern_pow, replace_fn_pow, proto_c)

    #GPT regex magic
    pattern = r"(?<!\[)\d+(\.\d+)?(?!\.)(?!\])|\.\d+(?!\])"
    replacement = lambda match: f"{float(match.group()):.1f}"
    proto_c = re.sub(pattern, replacement, proto_c)
    
    #needs a patch
    proto_c = proto_c.replace("d2.0", "d02")
    proto_c = proto_c.replace("c2.0", "c02")
    proto_c = proto_c.replace("b2.0", "b02")
    proto_c = proto_c.replace("a2.0", "a02")
    proto_c = proto_c.replace("b20.0", "b020")
    
    main_code += "\n\t" + "complex<double> "  + file_name + ";"
    main_code += "\n\t" + "complex<double>" + " Z[{0}];".format(num_var)
    main_code += "\n\t" + proto_c

    main_code += "return " + file_name  + ";"
        #if chosen_type == "f64":
        #    main_code += "*out = 0.5*(*out + conj(*out));"
        #else:
        #    main_code += "\n*out = 0.5*(*out + (*out).conj());"
    main_code = main_code.replace("MT*", "")

    main_code = main_code.replace("xs", "p_a.xs")
    main_code = main_code.replace("xt", "p_a.xt")
    main_code = main_code.replace("xu", "p_a.xu")
    
    main_code = main_code.replace("a02", "p_a.a02")
    
    main_code = main_code.replace("b02s", "p_a.b02s")
    main_code = main_code.replace("b02t", "p_a.b02t")
    main_code = main_code.replace("b02u", "p_a.b02u")
    main_code = main_code.replace("b020", "p_a.b020")

    main_code = main_code.replace("c02s", "p_a.c02s")
    main_code = main_code.replace("c02t", "p_a.c02t")
    main_code = main_code.replace("c02u", "p_a.c02u")

    main_code = main_code.replace("d02su", "p_a.d02su")
    main_code = main_code.replace("d02st", "p_a.d02st")
    main_code = main_code.replace("d02tu", "p_a.d02tu")
    
    template_code = template_code.replace("int {};".format(file_name), main_code)

cpp_file = open("form_factors_f64.cpp", "w")
cpp_file.write(template_code)
cpp_file.close()