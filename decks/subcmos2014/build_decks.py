#!/usr/bin/env python3
"""Generate isolated NMOS/PMOS process cuts from RIT Sub-CMOS slides.
All lengths are microns. Mask coordinates are assumed, not extracted artwork.
"""
import argparse
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
CHECKPOINTS = ['initial', 'pad', 'well_mask', 'wells', 'locos_mask', 'locos',
               'threshold', 'gate_oxide', 'poly_doped', 'gate', 'ldd',
               'spacers', 'source_drain', 'annealed', 'contacts', 'metal1']

class Deck:
    def __init__(self, device, pwell, post_locos_model, pmos_vt_dose=2.6e12,
                 pmos_vt_energy=35, blanket_vt_dose=7e11):
        self.device, self.pwell = device, pwell
        self.pmos_vt_dose, self.pmos_vt_energy = pmos_vt_dose, pmos_vt_energy
        self.blanket_vt_dose = blanket_vt_dose
        self.post_locos_model = post_locos_model
        self.lines = [f'# RIT Sub-CMOS 2014: {device.upper()}, isolated well cross-section.',
                      '# Source: rit_process/SubCmos2014_0p5um.pdf; see README.md.',
                      '# 1 um drawn gate (0.5 um layout lambda); assumed 2D masks.',
                      'option quiet', 'method fermi vertical grid.oxide=0.04',
                      # P diffusivity in SiO2: TSUPREM-4 2002.4 default (Table A-3,
                      # 4.338e10 um^2/min, 4.44 eV; Ghezzo & Brown, JECS 120,
                      # 146, 1973). Stock SUPREM-IV (7.6e-3 cm^2/s, 3.5 eV) is
                      # 5.5x faster at 1000 C, 27x at 800 C, and leaked n+ poly
                      # P through the 15 nm gate oxide into the NMOS channel.
                      # Segregation/transport already equal TSUPREM-4's.
                      'phosphorus oxide Dix.0=7.23 Dix.E=4.44']
        self.recipes = []
    def add(self, *lines): self.lines.extend(lines)
    def checkpoint(self, name):
        self.add(f'structure outfile={name}.str')
    def cut(self, material, left, right, top=-5, bottom=1):
        self.add(f'etch {material} start x={left:g} y={top:g}',
                 f'etch continue x={left:g} y={bottom:g}',
                 f'etch continue x={right:g} y={bottom:g}',
                 f'etch done x={right:g} y={top:g}')
    def resist(self, windows, planar=False):
        self.add('deposit photoresist file=resist_top.dat divisions=8' if planar else
                 'deposit photoresist thick=1.5 divisions=6')
        for lo, hi in windows: self.cut('photoresist',lo,hi)
    def implant(self, ion, dose, energy):
        self.add(f'implant {ion} dose={dose:g} energy={energy:g} pearson')
    def recipe(self, name, page, target, ramp_up, soak, cool, gas,
               anneal=5, flood=0, ramp_gas=None, pull=15):
        self.add(f'# PDF p.{page}: {name}; ramps discretized in <=50 C intervals.')
        segments = [(12,800,'nitrogen'),(15,800,'nitrogen')]
        def ramp(a,b,duration,ambient):
            count=max(1,math.ceil(abs(b-a)/50))
            segments.extend((duration/count,a+(b-a)*(i+.5)/count,ambient)
                            for i in range(count))
        ramp(800,target,ramp_up,ramp_gas or gas)
        if flood: segments.append((flood,target,'dryo2'))
        segments.append((soak,target,gas))
        if anneal: segments.append((anneal,target,'nitrogen'))
        ramp(target,800,cool,'nitrogen')
        segments.append((pull,800,'nitrogen'))
        previous=None
        for duration,temp,ambient in segments:
            continuation=previous is not None and not (
                ambient in ('dryo2','weto2') and previous=='nitrogen')
            self.add(f'diffuse time={duration:.10g} temp={temp:.10g} {ambient}'+
                     (' continue' if continuation else ''))
            previous=ambient
        self.recipes.append({'name':name,'pdf_page':page,'segments':segments})
    def build(self):
        for x,spacing,tag in [(0,.2,'left'),(.5,.08,''),(2,.08,''),(3,.08,''),
                              (5,.06,''),(5.5,.04,''),(6,.06,''),(6.5,.04,''),
                              (7,.06,''),(9,.08,''),(10,.2,'right')]:
            self.add(f'line x loc={x:g} spac={spacing:g}'+(f' tag={tag}' if tag else ''))
        for y,spacing,tag in [(0,.01,'top'),(.1,.025,''),(.4,.08,''),
                              (1,.2,''),(3,.4,''),(8,.8,'bottom')]:
            self.add(f'line y loc={y:g} spac={spacing:g}'+(f' tag={tag}' if tag else ''))
        self.add('region silicon xlo=left xhi=right ylo=top yhi=bottom',
                 'boundary exposed xlo=left xhi=right ylo=top yhi=top',
                 'boundary backside xlo=left xhi=right ylo=bottom yhi=bottom',
                 'initialize boron conc=1e15 orient=100')
        self.checkpoint('initial')
        self.recipe('250 pad oxide',24,1000,20,40,40,'dryo2',pull=12)
        self.checkpoint('pad')
        self.add('deposit nitride thick=0.15 divisions=5')
        if self.device=='pmos':
            self.add('etch nitride all')
            self.implant('phosphorus',9.5e12,125)
        self.recipe('330 n-well masking oxide',43,1000,20,50,40,'weto2',
                    flood=5,ramp_gas='nitrogen',pull=12)
        if self.device=='nmos': self.add('etch nitride all')
        self.add('# PDF p.48: 25 nm oxide loss in pre-strip BOE.',
                 'etch oxide dry thick=0.025')
        if self.pwell=='bf2': self.implant('bf2',8e12,110)
        else: self.implant('boron',2e13,35)
        self.checkpoint('well_mask')
        self.recipe('10 well drive',56,1100,30,600,60,'nitrogen',pull=12)
        self.checkpoint('wells')
        self.add('etch oxide all')
        self.recipe('250 second pad oxide',59,1000,20,40,40,'dryo2',pull=12)
        self.add('deposit nitride thick=0.15 divisions=5')
        # Keep nitride over the main active and a separate body-contact island.
        for lo,hi in [(-.1,.5),(2,3),(9,10.1)]: self.cut('nitride',lo,hi)
        if self.device=='nmos':
            self.resist([(-.1,.5),(2,3),(9,10.1)])
            self.implant('boron',8e13,100)
            self.add('etch photoresist all')
        self.checkpoint('locos_mask')
        self.add('# Stress coefficients from SUPREM example 9, not PDF calibration.',
                 'oxide Vc=300 Vr=30 Vd=25 stress.dep=true',
                 'method viscous oxide.rel=1e-2 grid.oxide=0.04')
        self.recipe('406 field oxidation',68,1100,30,65,55,'weto2',
                    flood=5,ramp_gas='nitrogen')
        self.add('etch nitride all','etch oxide dry thick=0.05')
        self.checkpoint('locos')
        if self.post_locos_model=='vertical':
            self.add('# Thin post-LOCOS oxides: vertical growth approximation.',
                     'method fermi vertical grid.oxide=0.04')
        self.recipe('311 sacrificial Kooi oxide',73,900,10,40,15,'weto2',
                    flood=5,ramp_gas='dryo2')
        self.implant('boron',self.blanket_vt_dose,60)  # PDF p.75: 7e11
        # PDF p.76: 2.6e12 at 60 keV. The p.74 sizing assumes half the dose
        # stops in the ~1000 A Kooi oxide; at 60 keV the simulation keeps ~95 %
        # and the PMOS is a normally-on buried channel (Vt +0.35 V). 35 keV
        # keeps the PDF dose (~67 % retained): Vt -0.885 V, SS 97 mV/dec.
        # See suprem4/doc/rit_cmos_full_flow_repair_20261119.md.
        if self.device=='pmos': self.implant('boron',self.pmos_vt_dose,self.pmos_vt_energy)
        self.checkpoint('threshold')
        # PDF p.78-79: the Kooi oxide (~1000 A target, 104 nm simulated) is
        # etched to clear and HF-dipped twice, so the gate oxide grows on bare
        # silicon. A timed strip with 10% overetch clears it; the nominal
        # 0.1 um etch left 4-5 nm of boron-bearing oxide under the gate oxide.
        self.add('etch oxide dry thick=0.11')
        # Soak calibrated (45 -> 40 min) so the as-grown gate oxide meets the
        # PDF's 150 A target on bare silicon (NMOS 14.8 / PMOS 15.1 nm).
        self.recipe('215 gate oxide',81,900,10,40,20,'dryo2',anneal=0)
        self.checkpoint('gate_oxide')
        # Resolve the 50 keV P profile (sigma=28.4 nm); spacing is a mesh
        # constraint, not an incident-dose adjustment. Keep a pre-anneal
        # checkpoint so implantation and subsequent transport can be audited.
        self.add('deposit poly thick=0.5 divisions=10 space=0.01')
        self.implant('phosphorus',2e16,50)
        self.add('structure outfile=poly_implanted.str')
        self.recipe('120 poly activation',88,1000,20,15,40,'nitrogen',anneal=0)
        self.checkpoint('poly_doped')
        self.add('etch poly left p1.x=5.5','etch poly right p1.x=6.5')
        self.checkpoint('gate')
        self.resist([(3,9)])
        self.implant('phosphorus' if self.device=='nmos' else 'boron',
                     2.5e13 if self.device=='nmos' else 4e13,
                     60 if self.device=='nmos' else 50)
        self.add('etch photoresist all')
        self.checkpoint('ldd')
        self.add('deposit oxide thick=0.5 divisions=10 space=0.04',
                 'etch oxide dry thick=0.5')
        self.checkpoint('spacers')
        self.resist([(3,9)],planar=True)
        self.implant('phosphorus' if self.device=='nmos' else 'boron',2e15,
                     60 if self.device=='nmos' else 50)
        self.add('etch photoresist all')
        self.resist([(.5,2)],planar=True)
        self.implant('boron' if self.device=='nmos' else 'phosphorus',2e15,
                     50 if self.device=='nmos' else 60)
        self.add('etch photoresist all')
        self.checkpoint('source_drain')
        self.recipe('280 source/drain activation',116,1000,20,20,40,'dryo2')
        self.checkpoint('annealed')
        self.add('# PDF p.117/118: use 0.5 um ILD, despite p.17 summary 0.4 um.',
                 'deposit oxide thick=0.5 divisions=10 space=0.06')
        self.resist([(.75,1.75),(3.75,4.75),(5.55,6.45),(7.25,8.25)],planar=True)
        self.add('etch oxide dry thick=1.5','etch photoresist all')
        self.checkpoint('contacts')
        self.add('deposit aluminum thick=0.75 divisions=10 space=0.06')
        for lo,hi in [(-.1,.5),(2,3.4),(5,5.25),(6.75,7),(8.8,10.1)]:
            self.cut('aluminum',lo,hi,top=-6,bottom=1)
        self.checkpoint('metal1')
        self.add(f'structure outfile={self.device}.str','quit')
        return '\n'.join(self.lines)+'\n'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pwell',choices=['bf2','boron'],default='bf2',
                        help='PDF p.54 improved BF2 or p.50 legacy boron process')
    parser.add_argument('--post-locos-model',choices=['vertical','viscous'],default='vertical',
                        help='thin-oxide growth model after stress-dependent LOCOS')
    parser.add_argument('--pmos-vt-dose',type=float,default=2.6e12,
                        help='PMOS-only Vt adjust boron dose (PDF p.76: 2.6e12)')
    parser.add_argument('--blanket-vt-dose',type=float,default=7e11,
                        help='blanket NMOS+PMOS Vt adjust boron dose at 60 keV (PDF p.75: 7e11)')
    parser.add_argument('--pmos-vt-energy',type=float,default=35,
                        help='PMOS-only Vt adjust boron energy in keV (PDF p.76: 60; '
                             'default 35 calibrated to the simulated Kooi retention)')
    parser.add_argument('--output',type=Path,default=HERE,
                        help='write a separate deck variant without replacing the main decks')
    args=parser.parse_args()
    output=args.output.resolve();output.mkdir(parents=True,exist_ok=True)
    metadata={'source':'rit_process/SubCmos2014_0p5um.pdf','pwell_variant':args.pwell,
              'gate_length_um':1.0,'post_locos_model':args.post_locos_model,
              'pmos_vt_dose_cm2':args.pmos_vt_dose,'pmos_vt_energy_keV':args.pmos_vt_energy,
              'blanket_vt_dose_cm2':args.blanket_vt_dose,
              'checkpoints':CHECKPOINTS,'devices':{}}
    for device in ['nmos','pmos']:
        deck=Deck(device,args.pwell,args.post_locos_model,args.pmos_vt_dose,
                  args.pmos_vt_energy,args.blanket_vt_dose)
        (output/(device+'.in')).write_text(deck.build())
        metadata['devices'][device]={'recipes':deck.recipes}
    (output/'resist_top.dat').write_text('0 -2\n10 -2\n')
    (output/'process.json').write_text(json.dumps(metadata,indent=2)+'\n')

if __name__=='__main__':main()
