import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,currentdir + "/rust_backend")
sys.path.insert(0,currentdir + "/../..")

import integrand
import deformation
import random

import loop_momenta_generator
import nloop_integrands

import madgraph.integrator.vectors as vectors


channel = 1
mu_sq = -1e-5
qs = [[0., 0., 0., 0.], [-495., 0., 0., 495.], [ 5., 26.97011261, 42.00346171, 992.50208264], [505., 0., 0., 495.]]

qs_l = [vectors.LorentzVector(q) for q in qs]
external_momenta = [qs_l[i] - qs_l[i - 1] for i in range(4)]
masses = [0, 0, 0, 0]

# fixed like this in C++
e_cm_sq = (external_momenta[0] + external_momenta[1]).square()

cppi = loop_momenta_generator.DeformationCPPinterface(qs_l)
cppi.set_option("Channel_ID", channel)
cppi.set_option("Mapping", 2) # 2 = Weinzierl
cppi.set_option("mu_UV_sq", mu_sq)

cpp_integrand = nloop_integrands.IntegrandCPPinterface("box1L_direct_integration_subtracted", qs)
cpp_integrand.set_option("CHANNEL_ID", channel)
cpp_integrand.set_option("UVSQ", mu_sq)

# defaults to Weinzierl mapping
parameterization = deformation.Parameterization(e_cm_sq=e_cm_sq, region=0, channel=channel, qs_py=qs)

deformation = deformation.Deformation(e_cm_sq=e_cm_sq, mu_sq=mu_sq, region=0, qs_py=qs, masses=masses)

box = integrand.Integrand("box1L_direct_integration_subtracted", channel=channel, region=0, mu_sq=mu_sq)

box.set_qs(qs)

random.seed(123)
for _ in range(10):
    point = [random.random(), random.random(), random.random(), random.random()]
    print("IN:", point)

    r, jac = parameterization.map(point)
    print("Mapped:", r, "jac:", jac)

    cpp_r, cpp_jac = cppi.hypcub_mapping(point)
    print("CPP", cpp_r, cpp_jac)

    print("Diff", r - cpp_r, jac - cpp_jac)  # should be 0

    cppi.deform_loop_momentum(r)
    cpp_out_r = cppi.get_deformed_loop_momentum()
    cpp_jac = cppi.get_jacobian()

    mapped_r, jac_real, jac_imag = deformation.deform(r)
    out_r = [ complex(real, imag) for real, imag in mapped_r ]

    print("Deformed:", out_r, "jac:", complex(jac_real, jac_imag))
    print("DeformedCPP:", cpp_out_r, "jac:", cpp_jac)

    print("Diff", [a - b for a, b in zip(out_r, cpp_out_r)], complex(jac_real, jac_imag) - cpp_jac)  # should be 0

    out_real, out_imag = box.evaluate(out_r)
    print("Result:", complex(out_real, out_imag))

    out_cpp = cpp_integrand.evaluate(out_r)
    print("ResultCPP:", out_cpp)
    print("Diff", out_cpp - complex(out_real, out_imag))
