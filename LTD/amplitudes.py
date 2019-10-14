#!/usr/bin/env python3

# This file contains the description of some aplitudes
# import topologies
import vectors
import os
import sys
import numpy as np
import copy
import yaml

params = {
    'g_f':   1.166390e-05,
    'alpha_s':   1.180000e-01,
    'alpha_ew':   1./1.325070e+02,
    'C_F':   4./3.,
    'q_u': 2. / 3.,
    'q_d': -1. / 3.,
}


class qqbar_diagram(object):
    def __init__(self, name, propagators, pows, chain, factor, ct, uv, positions=[]):
        """
            name:     given name of the diagram, usefull for debugging
            dens:     list with all the denominator for this diagram
            pows:     power of each propagator
            chain:    gamma chain contraction
            position: position of the vectors that depend on the loop momentum
            sing:     sign rising from taking the derivative wrt k0
            uv:       True/False if needs uv counterterms
            ct:       True/False if it is a counterterm and therefore comes with a minus sing
        """
        self.name = name
        self.dens = propagators
        self.pows = pows
        self.chain = chain
        self.positions = positions
        self.signs = 0
        self.factor = factor
        self.uv = uv
        self.ct = ct

    def to_flat_format(self):
        """ Turn this instance into a flat dictionary made out of simple lists or dictionaries only."""
        res = {}

        res["name"] = str(self.name)
        res["denominators"] = self.dens
        res["pows"] = self.pows
        res["chain"] = self.chain
        res["positions"] = self.positions
        res["loop_signature"] = self.signs
        res["factor"] = [float(self.factor.real), float(self.factor.imag)]
        res["ct"] = self.ct
        return res


class Amplitude(object):
    def __init__(self, topology, amp_type, polarizations=None, uv_pos=-1, mu_r_sq=1e2):
        """
            uv_pos:       position of the uv propagator
            mu_r_sq:  mu renormalization for the CT
        """

        with open("topologies.yaml", 'r') as stream:
            try:
                topologies = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)
        for top in topologies:
            if top['name'] == topology:
                mytop = top
                break
        else:
            raise AssertionError("Could not find topology %s" % topology)
        self.topology_name = topology
        self.type = amp_type
        self.ps = [vectors.LorentzVector(p)
                   for p in mytop['external_kinematics']]
        self.uv_pos = uv_pos
        if uv_pos >= 0:
            del self.ps[uv_pos]
        self.compute_invariants(polarizations)
        self.mu_r_sq = mu_r_sq
        self.diags = []

    def compute_invariants(self, polarizations):
        self.sij = {}
        for i, pi in enumerate(self.ps):
            for j, pj in enumerate(self.ps):
                self.sij['s%d%d' % (i+1, j+1)] = (pi+pj).square()

        # Polaizatoin vector contains a vector with all the polarization kinds
        # None, [vbar,+/-], [v,+/-], [ubar,+/-], [u,+/-], [a,+/-],
        if polarizations == None:
            self.polarizations = []
        else:
            self.polarizations = []
            for pol in polarizations:
                if any(pol[0] == p for p in ["u", "ubar", "v", "vbar", "a"]):
                    if any(pol[1] == s for s in ["+", "-"]):
                        self.polarizations += [pol[0]+pol[1]]
                    else:
                        raise BaseException(
                            "No polarization known with this sign %s" % str(pol[1]))
                else:
                    raise BaseException(
                        "No polarization known with this name %s" % pol[0])

    def add_born(self, chain, factor):
        if self.type == 'qqbar_photons':
            diag = qqbar_diagram("born", [], [], chain, factor, False, False)
            diag.positions = []
            diag.signs = 0
            self.diags += [diag]
        else:
            raise BaseException(
                "Not knonw diagram structure for %s" % self.type)

    def create_amplitude(self, name, loops, diags, vectors):
        self.name = name
        self.vectors = vectors
        # Find position in chain of loop dependent pieces
        loop_dependent = [
            i+1 for (i, v) in enumerate(vectors) if len(v[0]) > 0]
        for diag in diags:
            diag.positions = [pos for (pos, v_id) in enumerate(diag.chain) if any(v_id == i for i in loop_dependent)]
            diag.loop_signature = vectors[diag.positions[0]-1][0][0]

        # Add new diagram to the amplitude
        self.diags.extend(diags)

        if self.type == 'qqbar_photons':
            # Create UV CT
            for diag in diags:
                if diag.uv and self.uv_pos >= 0:
                    # Remove Leading UV pole
                    uv_diag = copy.deepcopy(diag)
                    uv_diag.name = "UV" + uv_diag.name
                    uv_diag.dens = [self.uv_pos]
                    uv_diag.ct = not uv_diag.ct
                    # Rise power of the unique propagator
                    uv_diag.pows = [len(diag.pows)]
                    for i in uv_diag.positions:
                        uv_diag.chain[i] = self.uv_pos
                    self.diags += [uv_diag]

                    # Check if a bubble
                    if len(diag.dens) == 2:
                        pair = [uv_diag.chain[uv_diag.positions[0]],
                                uv_diag.chain[uv_diag.positions[0]+2]]
                        rpair = [i for i in reversed(pair)]
                        for (pp, s) in zip([pair, rpair], ["a", "b"]):
                            uv_log_diag = copy.deepcopy(uv_diag)
                            uv_log_diag.name += "_LOG"+s
                            uv_log_diag.ct = not uv_log_diag.ct
                            # Insert momenta in the gamma chain
                            for v in pp:
                                uv_log_diag.chain.insert(
                                    uv_log_diag.positions[0], v)
                            # Update positions
                            if s == "b":
                                uv_log_diag.positions = [
                                    uv_log_diag.positions[0], uv_log_diag.positions[0]+2]
                            else:
                                # uv_log_diag.positions = [
                                #    uv_log_diag.positions[0]+1, uv_log_diag.positions[0]+2]
                                continue
                            # Rise power of the unique propagator
                            uv_log_diag.pows = [len(diag.pows)+1]

                            self.diags += [uv_log_diag]

            # Create UV Approximations

            # Number of diagrams for the approximation
            n_UV_LO_diags = 0
            for diag in diags:
                if self.uv_pos >= 0:
                    # Approximate Leading UV pole
                    uv_diag = copy.deepcopy(diag)
                    uv_diag.name += "UV_LO"
                    uv_diag.dens = [self.uv_pos-1]
                    uv_diag.pows = [len(diag.pows)]
                    self.diags += [uv_diag]
                    n_UV_LO_diags += 1

                    # Check if a bubble
                    if len(diag.dens) == 2:
                        pair = [self.uv_pos,
                                uv_diag.chain[uv_diag.positions[0]+2]]
                        rpair = [i for i in reversed(pair)]
                        for (pp, s) in zip([pair, rpair], ["a", "b"]):
                            uv_log_diag = copy.deepcopy(uv_diag)
                            uv_log_diag.name += "_LOG"+s
                            uv_log_diag.ct = not uv_log_diag.ct
                            # Insert momenta in the gamma chain
                            for v in pp:
                                uv_log_diag.chain.insert(
                                    uv_log_diag.positions[0], v)
                            # Update positions
                            if s == "b":
                                uv_log_diag.positions = [
                                    uv_log_diag.positions[0], uv_log_diag.positions[0]+2]
                            else:
                                uv_log_diag.positions = [
                                    uv_log_diag.positions[0]+1, uv_log_diag.positions[0]+2]
                            # Rise power of the unique propagator
                            uv_log_diag.pows = [len(diag.pows)+1]

                            self.diags += [uv_log_diag]
                            n_UV_LO_diags += 1

            # Set of diagrams to compute the amplitude and to do the same using the UV approximation
            set_amp = [d.name for d in self.diags if not "UV_LO" in d.name and not "born" in d.name]
            set_uv = [d.name for d in self.diags if "UV" in d.name]
            self.sets = [set_amp, set_uv]

            # Print all diagrams
            # for diag in self.diags:
            #    print diag.to_flat_format()
        else:
            print("NO SUCH AMPLITUDE! %s" % name)

    def to_flat_format(self):
        """ Turn this instance into a flat dictionary made out of simple lists or dictionaries only."""

        res = {}

        res['name'] = self.name
        res['amp_type'] = self.type
        res['topology'] = self.topology_name
        res['diagrams'] = [diag.to_flat_format() for diag in self.diags]
        res['vectors'] = [
            [l, [float(v) for v in vec]] for (l, vec) in self.vectors]
        res['ps'] = [[float(v) for v in vec] for vec in self.ps]
        res['pols_type'] = self.polarizations
        res['sets'] = self.sets
        res['mu_r_sq'] = self.mu_r_sq

        return res


class AmplitudesCollection(dict):
    def add(self, amplitude, entry_name=None):
        assert(isinstance(amplitude, Amplitude))
        if entry_name == None:
            if amplitude.name == None:
                raise BaseException(
                    "Specify the name of the hard-coded amplitude to add to the collection.")
            self[amplitude.name] = amplitude
        else:
            self[entry_name] = amplitude

    def export_to(self, output_path, format='yaml'):
        """ Exports this amplitude to a given format."""

        export_format = format.lower()
        allowed_export_format = ['yaml']
        if export_format not in ['yaml']:
            raise BaseException("Amplitude can only be exported in the following formats: %s" % (
                ', '.join(allowed_export_format)))

        if export_format == 'yaml':
            try:
                import yaml
                from yaml import Loader, Dumper
            except ImportError:
                raise BaseException(
                    "Install yaml python module in order to export amplitudes to yaml.")

        flat_record = [self[k].to_flat_format() for k in sorted(self.keys())]
        if output_path is not None:
            open(output_path, 'w').write(
                yaml.dump(flat_record, Dumper=Dumper))
        else:
            return yaml.dump(flat_record, Dumper=Dumper)


if __name__ == "__main__":
    amplitudes_collection = AmplitudesCollection()

    # =================== add amplitude ddAAA ======================= #
    tree_factor = params['alpha_ew']**1.5*params['q_d']**3 * (4.0 * np.pi)**1.5
    # Initialize
    amp = Amplitude("manual_uuWWZ_amplitude_P2",      # name of topology
                    'qqbar_photons',                  # process type
                    zip(["u", "vbar", "a", "a", "a"], # polarizations
                        ["+", "-", "+", "+", "+"]),
                    1,                                # uv_pos
                    91.188)                           # mu_r_sq
    # Born Level Diagram
    born_factor = tree_factor / amp.sij['s23'] / amp.sij['s15']
    amp.add_born([8, 6, 9, 7, 10], born_factor)
    # Diagrams and vectors
    factor = -1j * params['C_F'] * tree_factor\
        * params['alpha_s'] * (4.0 * np.pi)
    amp.create_amplitude(
        'ddAAA',
        1,
        [qqbar_diagram("D1", [0, 3, 4, 5], [1, 1, 1, 1], [8, 6, -1, 3, 9, 4, 10, 5, -1], factor/amp.sij["s23"], False, False),
         qqbar_diagram("D2", [0, 4, 5], [1, 1, 1], [8, 6, 9, 7, -1, 4, 10, 5, -1], factor/amp.sij["s23"]/amp.sij["s15"], False, True,),
         qqbar_diagram("D3", [0, 4], [1, 1], [8, 6, 9, 7, -1, 4, -1, 7, 10], factor/amp.sij["s23"]/amp.sij["s15"]**2, False, True,),
         qqbar_diagram("D4", [0, 3, 4], [1, 1, 1], [8, 6, -1, 3, 9, 4, -1, 7, 10], factor/amp.sij["s23"]/amp.sij["s15"], False, True,),
         qqbar_diagram("D5", [0, 2, 3, 4], [1, 1, 1, 1], [-1, 2, 8, 3, 9, 4, -1, 7, 10], factor/amp.sij["s15"], False, False,),
         qqbar_diagram("D6", [0, 2, 3], [1, 1, 1], [-1, 2, 8, 3, -1, 6, 9, 7, 10], factor/amp.sij["s23"]/amp.sij["s15"], False, True,),
         qqbar_diagram("D7", [0, 3], [1, 1], [8, 6, -1, 3, -1, 6, 9, 7, 10], factor/amp.sij["s23"]**2/amp.sij["s15"], False, True,),
         qqbar_diagram("D8", [0, 2, 3, 4, 5], [1, 1, 1, 1, 1], [-1, 2, 8, 3, 9, 4, 10, 5, -1], factor, False, False,),
         qqbar_diagram("IR", [0, 2, 5], [1, 1, 1], [-1, 2, 8, 6, 9, 7, 10, 5, -1], factor/amp.sij["s23"]/amp.sij["s15"], True, True,),
         ],
        # Vectors [loopmomenta, externals]
        [
            [[-1], -amp.ps[0]],
            [[-1], -amp.ps[0]-amp.ps[1]],
            [[-1], amp.ps[3]+amp.ps[4]],
            [[-1], amp.ps[4]],
            [[-1], vectors.LorentzVector([0, 0, 0, 0])],
            [[], -amp.ps[1]-amp.ps[2]],
            [[], amp.ps[0]+amp.ps[4]],
        ],
    )
    # Store
    amplitudes_collection.add(amp)

    # =================== add amplitude ddAA ======================= #
    tree_factor = params['alpha_ew'] * params['q_d']**2 * (4.0 * np.pi)
    # Initialize
    amp = Amplitude("manual_eeAA_amplitude_P2",  # name of topology
                    'qqbar_photons',             # process type
                    zip(["u", "vbar", "a", "a"], # polarizations
                        ["+", "-", "+", "+"]),
                    1,                           # uv_pos
                    91.188)                      # mu_r_sq
    # Born Level Diagram
    born_factor = tree_factor / amp.sij['s23']
    amp.add_born([6, 5, 7], born_factor)
    # Diagrams and vectors
    factor = -1j*params['C_F']*tree_factor * params['alpha_s'] * (4.0 * np.pi)
    amp.create_amplitude(
        'ddAA',
        1,
        [qqbar_diagram("D1", [0, 3, 4], [1, 1, 1], [6, 5, -1, 3, 7, 4, -1], factor/amp.sij["s23"], False, True,),
         qqbar_diagram("D2", [0, 2, 3], [1, 1, 1], [-1, 2, 6, 3, -1, 5, 7], factor/amp.sij["s23"], False, True),
         qqbar_diagram("D3", [0, 2, 3, 4], [1, 1, 1, 1], [-1, 2, 6, 3, 7, 4, -1], factor, False, False,),
         qqbar_diagram("D4", [0, 3], [1, 1], [6, 5, -1, 3, -1, 5, 7], factor/amp.sij["s23"]/amp.sij["s23"], False, True),
         qqbar_diagram("IR", [0, 2, 4], [1, 1, 1], [-1, 2, 6, 5, 7, 4, -1], factor / amp.sij["s23"], True, True),
         ],
        # Vectors [loopmomenta, externals]
        [
            [[-1], -amp.ps[0]],
            [[-1], -amp.ps[0]-amp.ps[1]],
            [[-1], amp.ps[3]],
            [[-1], vectors.LorentzVector([0, 0, 0, 0])],
            [[], -amp.ps[1]-amp.ps[2]],
        ],
    )
    # Store
    amplitudes_collection.add(amp)

    # Export amplitudes
    root_dir = ''
    amplitudes_collection.export_to(os.path.join(root_dir, 'amplitudes.yaml'))
