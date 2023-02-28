#!/usr/bin/env python3

# Example : 
# python3 gen_gm_id_plots.py sky130_fd_pr__nfet_01v8 ../skywater-pdk-libs-sky130_fd_pr/cells/nfet_01v8/sky130_fd_pr__nfet_01v8.bins.csv sky130_fd_pr__nfet_01v8_%s.png

import PySpice.Logging.Logging as Logging
logger = Logging.setup_logging()
from PySpice.Spice.Netlist import Circuit
from PySpice.Unit import *
import matplotlib.pyplot as plt
import numpy as np
import os.path
import argparse

def create_test_circuit(fet_type, iparam, fet_L, fet_W, coner_path):
    c=Circuit('gm_id')
    c.include(coner_path)

    # create the circuit
    c.V('gg', 1, c.gnd, 0@u_V)
    c.V('dd', 2, c.gnd, 1.8@u_V)
    c.X('M1', fet_type, 2, 1, c.gnd, c.gnd, L=fet_L, W=fet_W, ad="'W*0.29'", pd="'2*(W+0.29)'", as_="'W*0.29'", ps="'2*(W+0.29)'", nrd="'0.29/W'", nrs="'0.29/W'", sa=0, sb=0, sd=0, nf=1, mult=1)
    return c

def run_sim(c, iparam, fet_W):
    sim = c.simulator()
    sim.save_internal_parameters(iparam%'gm', iparam%'id', iparam%'gds', iparam%'cgg')

    # run the dc simulation
    an = sim.dc(Vgg=slice(0, 1.8, 0.01))

    # calculate needed values..need as_ndarray() since most of these have None as the unit and that causes an error
    gm = an.internal_parameters[iparam%'gm'].as_ndarray()
    id = an.internal_parameters[iparam%'id'].as_ndarray()
    gm_id = gm / id
    cgg = an.internal_parameters[iparam%'cgg'].as_ndarray()
    ft = gm / cgg
    id_W = id / fet_W
    gds = an.internal_parameters[iparam%'gds'].as_ndarray()
    gm_gds = gm / gds

    return gm_id, ft, id_W, gm_gds, an.nodes['v-sweep'], gm, id, cgg, gds

def init_plots():
    figs = [plt.figure(), plt.figure(), plt.figure(), plt.figure()]
    plts = [f.subplots() for f in figs]
    figs[0].suptitle('Id/W vs gm/Id')
    plts[0].set_xlabel("gm/Id")
    plts[0].set_ylabel("Id/W")
    figs[1].suptitle('fT vs gm/Id')
    plts[1].set_xlabel("gm/Id")
    plts[1].set_ylabel("f_T")
    figs[2].suptitle('gm/gds vs gm/Id')
    plts[2].set_xlabel("gm/Id")
    plts[2].set_ylabel("gm/gds")
    figs[3].suptitle('gm/Id vs Vgg')
    plts[3].set_xlabel("Vgg")
    plts[3].set_ylabel("gm/Id")
    return figs, plts

def gen_plots(gm_id, id_W, ft, gm_gds, vsweep, fet_W, fet_L, plts):
    # plot some interesting things
    plts[0].plot(gm_id, id_W, label=f'W {fet_W} x L {fet_L}')
    plts[1].plot(gm_id, ft, label=f'W {fet_W} x L {fet_L}')
    plts[2].plot(gm_id, gm_gds, label=f'W {fet_W} x L {fet_L}')
    plts[3].plot(vsweep, gm_id, label=f'W {fet_W} x L {fet_L}')

def read_bins(fname):
    import csv
    r=csv.reader(open(fname, 'r'))
    return r

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate gm/id plots for a given fet type.')
    parser.add_argument('--fet_type', help='FET type to simulate')
    parser.add_argument('--bins_csv', help='bins csv file')
    parser.add_argument('--out_file', help='template with 1 %%s which will contain the plot name. 4 are generated per LxW combo. (Ex:sky130_fd_pr__nfet_01v8_%%s.png)')
    parser.add_argument('--width', help='If specified, only W/L pairs for that width are processed.', required=False)
    parser.add_argument('--corner_path', help='the path to the corner file to use')

    args = parser.parse_args()

    fet_type = args.fet_type
    bins_fname = args.bins_csv
    figname = args.out_file
    only_W = None
    if args.width:
        only_W = float(args.width)
    print(f'Simulating {fet_type} with bins {bins_fname}')

    corner_path = args.corner_path

    iparam = f'@m.xm1.m{fet_type}[%s]'
    fet_L = 0.15
    fet_W = 1
    c = create_test_circuit(fet_type, iparam, fet_L, fet_W, corner_path)
    bins = read_bins(bins_fname)
    next(bins)

    figtitles = ['Id_w', 'fT', 'gm_gds', 'gm_id']
    figs, plts = init_plots()
    import h5py
    h5name = os.path.splitext(figname % 'data')[0] + '.h5'
    out = h5py.File(h5name, "w")
    bins_d = out.create_dataset('bins', (0, 2), maxshape=(None,2))
    gm_d = out.create_dataset('gm', (0, 0), maxshape=(None,None))
    id_d = out.create_dataset('id', (0, 0), maxshape=(None,None))
    cgg_d = out.create_dataset('cgg', (0, 0), maxshape=(None,None))
    gds_d = out.create_dataset('gds', (0, 0), maxshape=(None,None))
    vsweep_d = out.create_dataset('vsweep', (0,0), maxshape=(None,None))
    idx = 0
    for dev, bin, fet_W, fet_L in bins:
        fet_W, fet_L = float(fet_W), float(fet_L)
        if only_W is not None and fet_W != only_W:
            continue
        print(f'{bin}: {dev}  W {fet_W} x L {fet_L}')
        c.element('XM1').parameters['W'] = fet_W
        c.element('XM1').parameters['L'] = fet_L
        gm_id, ft, id_W, gm_gds, vsweep, gm, id, cgg, gds = run_sim(c, iparam, fet_W)
        if idx == 0:
            gm_d.resize(len(gm_id), 1)
            id_d.resize(len(id_W), 1)
            cgg_d.resize(len(ft), 1)
            gds_d.resize(len(gm_gds), 1)
            vsweep_d.resize(len(vsweep), 1)
        bins_d.resize(idx+1, 0)
        gm_d.resize(idx+1, 0)
        id_d.resize(idx+1, 0)
        cgg_d.resize(idx+1, 0)
        gds_d.resize(idx+1, 0)
        vsweep_d.resize(idx+1, 0)
        bins_d[idx,:] = [fet_W, fet_L]
        gm_d[idx, :] = gm
        id_d[idx, :] = id
        cgg_d[idx, :] = cgg
        gds_d[idx, :] = gds
        vsweep_d[idx, :] = vsweep  # should be the same for every row
        idx += 1
        gen_plots(gm_id, id_W, ft, gm_gds, vsweep, fet_W, fet_L, plts)
    for f,nm in zip(figs, figtitles):
        f.legend()
        f.tight_layout()
        f.savefig(figname % nm)
