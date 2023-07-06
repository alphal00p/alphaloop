import re
import os

file_names = ["astu", "atsu", "aust", "bstu", "btsu", "bust", "cstu"]
function_names = ["APHOAMPFFSTU", "APHOAMPFFTSU", "APHOAMPFFUST", "BPHOAMPFFSTU", "BPHOAMPFFTSU", "BPHOAMPFFUST", "CPHOAMPFFSTU"]

common_denominators = ["(96.0*xs*xs*xt*xt*xu*xu*E1*E2*E3*E42*s*s)","(96.0*xs*xs*xt*xt*xu*xu*E1*E2*E3*E42*t*t)","(96.0*xs*xs*xt*xt*xu*xu*E1*E2*E3*E42*u*u)", "(48.0*xs*xs*xs*xt*xt*xt*xu*xu*xu*(E1*E2*E32*E42*s*s*t*t*u))","(48.0*xs*xs*xs*xt*xt*xt*xu*xu*xu*(E12*E2*E3*E42*s*s*t*t*u))","(48.0*xs*xs*xs*xt*xt*xt*xu*xu*xu*(E1*E22*E3*E42*s*s*t*u*u))","(3.0*xs*xs*xs*xs*xt*xt*xt*xt*xu*xu*xu*xu*(E1*E2*E3*E4*s*t*t*u))"]
types = ["complex<double>", "complex128", "mppp::complex"]
type_labels = ["f64", "f128", "mpfr"]
one_loop_functions = ["a02", "b02s", "b02t", "b02u","c02s", "c02t", "c02u","b020", "d02su", "d02st", "d02tu", "d02us", "d02ut", "d02ts"]

function_dict = {
    "a02": "OLO_A0c(&rslt, &M2_f64)", 
    "b02s": "OLO_B0cc(&rslt, &s_f64, &M2_f64, &M2_f64)",
    "b02t": "OLO_B0cc(&rslt, &t_f64, &M2_f64, &M2_f64)",
    "b02u": "OLO_B0cc(&rslt, &u_f64, &M2_f64, &M2_f64)",
    "c02s": "OLO_C0cc(&rslt, &zero_f64, &zero_f64, &s_f64, &M2_f64, &M2_f64, &M2_f64)",
    "c02t": "OLO_C0cc(&rslt, &zero_f64, &zero_f64, &t_f64, &M2_f64, &M2_f64, &M2_f64)",
    "c02u": "OLO_C0cc(&rslt, &zero_f64, &zero_f64, &u_f64, &M2_f64, &M2_f64, &M2_f64)",
    "b020": "OLO_B0cc(&rslt, &zero_f64, &M2_f64, &M2_f64)",
    "d02su": "OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &s_f64, &u_f64, &M2_f64, &M2_f64, &M2_f64, &M2_f64)",
    "d02st": "OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &s_f64, &t_f64, &M2_f64, &M2_f64, &M2_f64, &M2_f64)",
    "d02tu": "OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &t_f64, &u_f64, &M2_f64, &M2_f64, &M2_f64, &M2_f64)",
    "d02us": "OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &u_f64, &s_f64, &M2_f64, &M2_f64, &M2_f64, &M2_f64)",
    "d02ut": "OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &u_f64, &t_f64, &M2_f64, &M2_f64, &M2_f64, &M2_f64)",
    "d02ts": "OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &t_f64, &s_f64, &M2_f64, &M2_f64, &M2_f64, &M2_f64)"
}

function_dict_massless = {
    "a02": "OLO_A0c(&rslt, &zero_f64)", 
    "b02s": "OLO_B0cc(&rslt, &s_f64, &zero_f64, &zero_f64)",
    "b02t": "OLO_B0cc(&rslt, &t_f64, &zero_f64, &zero_f64)",
    "b02u": "OLO_B0cc(&rslt, &u_f64, &zero_f64, &zero_f64)",
    "c02s": "OLO_C0cc(&rslt, &zero_f64, &zero_f64, &s_f64, &zero_f64, &zero_f64, &zero_f64)",
    "c02t": "OLO_C0cc(&rslt, &zero_f64, &zero_f64, &t_f64, &zero_f64, &zero_f64, &zero_f64)",
    "c02u": "OLO_C0cc(&rslt, &zero_f64, &zero_f64, &u_f64, &zero_f64, &zero_f64, &zero_f64)",
    "b020": "OLO_B0cc(&rslt, &zero_f64, &zero_f64, &zero_f64)",
    "d02su": "OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &s_f64, &u_f64, &zero_f64, &zero_f64, &zero_f64, &zero_f64)",
    "d02st": "OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &s_f64, &t_f64, &zero_f64, &zero_f64, &zero_f64, &zero_f64)",
    "d02tu": "OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &t_f64, &u_f64, &zero_f64, &zero_f64, &zero_f64, &zero_f64)",
    "d02us": "OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &u_f64, &s_f64, &zero_f64, &zero_f64, &zero_f64, &zero_f64)",
    "d02ut": "OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &u_f64, &t_f64, &zero_f64, &zero_f64, &zero_f64, &zero_f64)",
    "d02ts": "OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &t_f64, &s_f64, &zero_f64, &zero_f64, &zero_f64, &zero_f64)"
}

scaling_dict = {
    "a02": "1.0 / (M_f64 * M_f64)", 
    "b02s": "1.0",
    "b02t": "1.0",
    "b02u": "1.0",
    "c02s": "M_f64*M_f64",
    "c02t": "M_f64*M_f64",
    "c02u": "M_f64*M_f64",
    "b020": "1.0",
    "d02su": "M_f64*M_f64*M_f64*M_f64",
    "d02st": "M_f64*M_f64*M_f64*M_f64",
    "d02tu": "M_f64*M_f64*M_f64*M_f64",
    "d02us": "M_f64*M_f64*M_f64*M_f64",
    "d02ut": "M_f64*M_f64*M_f64*M_f64",
    "d02ts": "M_f64*M_f64*M_f64*M_f64"
}

# Define the replacement function
def replace_fn_array(match):
    n = int(match.group(1))
    return f"Z[{n-1}]"

def replace_fn_pow(match):
    x = match.group(1)
    n = int(match.group(2))
    return '*'.join([x]*n)

source_file = open("form_factors_stu.cpp", "r").read()
main_code = source_file

for i in range(len(file_names)):
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
    
    for float_type in zip(types, type_labels):
        main_code += """void {0}_{1}({2} E1,
                       {2} E2,
                       {2} E3,
                       {2} p1_p2,
                       {2} p1_p3,
                       {2} p2_p3,
                       {2}* out)""".format(function_name, float_type[1], float_type[0])
        main_code += "\n{\n"
        main_code += "\t {0} s,t,u,E4,E12,E22,E32,E42;".format(float_type[0])

        main_code += "\n\t complex<double> s_f64, t_f64, u_f64;"
        main_code += "\n\t {0} xs,xt,xu;".format(float_type[0])

        main_code += "\n"

        main_code += "\n\t s = 2.0*p1_p2;"
        main_code += "\n\t t = 2.0*p2_p3;"
        main_code += "\n\t u = 2.0*p1_p3;"
        main_code += "\n\t E4 = -E1-E2-E3;"
        main_code += "\n\t E12 = E1*E1;"
        main_code += "\n\t E22 = E2*E2;"
        main_code += "\n\t E32 = E3*E3;"
        main_code += "\n\t E42 = E4*E4;"

        main_code += "\n"
        
        main_code += "\n\t s_f64 = (complex<double>)s;"
        main_code += "\n\t t_f64 = (complex<double>)t;"
        main_code += "\n\t u_f64 = (complex<double>)u;"

        main_code += "\n\t OLO_SCALE(&scale);"

        main_code += "\n"

        main_code += "\n\t xs = s / (M_{0} * M_{0});".format(float_type[1])
        main_code += "\n\t xt = t / (M_{0} * M_{0});".format(float_type[1])
        main_code += "\n\t xu = u / (M_{0} * M_{0});".format(float_type[1])
        main_code += "\n"
        main_code += "\n\t complex<double> M2_f64 = M_f64*M_f64;"
        main_code += "\n"
        main_code += "\n\t complex<double> rslt[3];" 
        for function in one_loop_functions:
            if function in proto_c:
                main_code += "\n\t" + function_dict[function].replace("type", float_type[1]) + ";"
                main_code += "\n\t" + float_type[0] + " " + function + " = ({0}){1}".format(float_type[0], "rslt[0]*{0};".format(scaling_dict[function]))
                main_code += "\n"

        main_code += "\n\t" + float_type[0] + " " + file_name + ";"
        main_code += "\n\t" + float_type[0] + " Z[{0}];".format(num_var)
        main_code += "\n\t" + proto_c

        main_code += "*out =  " + "2.0*3.0*top_Q8_{0}/(M_{0}*M_{0}*M_{0}*M_{0}*M_PI*M_PI) * ".format(float_type[1]) + file_name + "/" + common_denominators[i] + ";"
        #if float_type[1] == "f64":
        #    main_code += "*out = 0.5*(*out + conj(*out));"
        #else:
        #    main_code += "\n*out = 0.5*(*out + (*out).conj());"
        main_code += "\n}\n"    
        main_code = main_code.replace("MT", "M_{0}".format(float_type[1]))

os.remove("form_factors.cpp")
cpp_file = open("form_factors.cpp", "a")
cpp_file.write(main_code)
cpp_file.close()
os.system("g++ -c form_factors.cpp -L./ -lavh_olo -lgfortran -fPIC -O2")
os.system("ar -rf form_factors.a libavh_olo.a form_factors.o")



