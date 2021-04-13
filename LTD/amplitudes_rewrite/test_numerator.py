import amplitudes as amp

runcard = "/home/armin/my_programs/MG5_py_3_alphaloop_master/PLUGIN/alphaloop/LTD/amplitudes_rewrite/runcard_template_numerator.yaml"

exporter = amp.AmpExporter.from_runcard(runcard)
exporter.export()