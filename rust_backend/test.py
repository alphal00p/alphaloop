import integrand
import deformation
import random

qs = [[0., 0., 0., 0.], [-495., 0., 0., 495.], [ 5., 26.97011261, 42.00346171, 992.50208264], [505., 0., 0., 495.]]
masses = [0, 0, 0, 0]

# defaults to Weinzierl mapping
parameterization = deformation.Parameterization(e_cm_sq=100, region=0, channel=1, qs_py=qs)

deformation = deformation.Deformation(e_cm_sq=100, mu_sq=-1000, region=0, qs_py=qs, masses=masses)

box = integrand.Integrand("box1L_direct_integration_subtracted", channel=1, region=0, mu_sq=1e-5)

box.set_qs(qs)

for _ in range(1000):
    point = [random.random(), random.random(), random.random(), random.random()]
    print("IN:", point)

    r, jac = parameterization.map(point)
    print("Mapped:", r, "jac:", jac)

    mapped_r, jac_real, jac_imag = deformation.deform(r)
    out_r = [ complex(real, imag) for real, imag in mapped_r ]

    print("Deformed:", out_r, "jac:", complex(jac_real, jac_imag))

    out_real, out_imag = box.evaluate(out_r)
    print("Result:", complex(out_real, out_imag))
