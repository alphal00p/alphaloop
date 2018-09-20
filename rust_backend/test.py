import integrand

box = integrand.Integrand("box1L_direct_integration_subtracted", channel=1, region=0, mu_sq=1e-5)

box.set_qs([[0., 0., 0., 0.], [-495., 0., 0., 495.], [ 5., 26.97011261, 42.00346171, 992.50208264], [505., 0., 0., 495.]])

print(box.evaluate([3.0+2j, 60.0+8j,2.0+2j,3.0+0j]))
