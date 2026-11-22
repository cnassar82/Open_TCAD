#!/usr/bin/env python3
"""Generate isolated NMOS/PMOS process cuts from the RIT Advanced CMOS slides.

Source: rit_process/AdvCmos2014_0p25um.pdf (Lynn Fuller, August 17, 2014).
lambda = 0.25 um, Lpoly = 0.5 um. All lengths are microns. Mask coordinates
are assumed, not extracted artwork. See README.md for every PDF reference and
every assumption the PDF does not settle.
"""
import argparse
import os
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
CHECKPOINTS = ['initial', 'pad', 'active_mask', 'first_oxide', 'recessed_oxide',
               'wells', 'threshold', 'gate_oxide', 'poly', 'gate', 'reox', 'ldd',
               'spacers', 'source_drain', 'annealed', 'contacts', 'metal1']

# Layout (um). Field windows open the nitride; everything else is active.
# The device active is 7 um wide: with a 3 um active the field oxide blocked
# enough of the well implants, and the 6 h drive spread them far enough
# sideways, to cut the channel n-well by ~25 % (doc rit_advcmos_first_flow).
FIELD = [(-.1, .3), (1.5, 1.9), (8.5, 9.6)]
BODY_ACTIVE = (.3, 1.5)   # 1.2 um: a 0.8 um window is nearly filled by the 0.35 um spacer nitride
DEVICE_ACTIVE = (1.9, 8.5)
# PDF p.4: Lmin (drawn) 0.5 um; Lpoly 0.40 um after etch, 0.35 um after re-oxidation.
# ADVCMOS_LPOLY overrides the as-etched poly length (centered on x = 5.0).
_LPOLY = float(os.environ.get('ADVCMOS_LPOLY', '0.5'))
GATE = (round(5.0 - _LPOLY / 2, 4), round(5.0 + _LPOLY / 2, 4))
CONTACTS = {'body': (.6, 1.2), 'source': (3.6, 4.2), 'gate': (4.85, 5.15),
            'drain': (5.8, 6.4)}
METAL_GAPS = [(-.1, .3), (1.55, 2.0), (4.45, 4.6), (5.4, 5.55), (8.5, 9.6)]
PROBES = {'body': .9, 'source': 3.9, 'channel': 5.0, 'drain': 6.1}
RIGHT = 9.5


class Deck:
    def __init__(self, device, rox_soak, gate_soak, n2o_as, pad_soak, poly_reox='none',
                 poly_reox_ox=0.05):
        self.device = device
        self.poly_reox = poly_reox
        self.poly_reox_ox = poly_reox_ox
        self.rox_soak, self.gate_soak, self.n2o_as = rox_soak, gate_soak, n2o_as
        self.pad_soak = pad_soak
        self.lines = [f'# RIT Advanced CMOS 2014: {device.upper()}, isolated well cross-section.',
                      '# Source: rit_process/AdvCmos2014_0p25um.pdf; see README.md.',
                      '# 0.5 um drawn poly gate (0.25 um lambda); assumed 2D masks.',
                      'option quiet', 'method fermi vertical grid.oxide=0.02',
                      # TSUPREM-4 2002.4 default P diffusivity in SiO2; see
                      # decks/subcmos2014/build_decks.py for the derivation.
                      'phosphorus oxide Dix.0=7.23 Dix.E=4.44']
        self.recipes = []

    def add(self, *lines): self.lines.extend(lines)
    def checkpoint(self, name): self.add(f'structure outfile={name}.str')

    def cut(self, material, left, right, top=-5, bottom=1):
        self.add(f'etch {material} start x={left:g} y={top:g}',
                 f'etch continue x={left:g} y={bottom:g}',
                 f'etch continue x={right:g} y={bottom:g}',
                 f'etch done x={right:g} y={top:g}')

    def resist(self, windows, planar=False):
        self.add('deposit photoresist file=resist_top.dat divisions=8' if planar else
                 'deposit photoresist thick=1.5 divisions=6')
        for lo, hi in windows: self.cut('photoresist', lo, hi)

    def implant(self, ion, dose, energy):
        self.add(f'implant {ion} dose={dose:g} energy={energy:g} pearson')

    def recipe(self, name, page, target, intervals, push=12, stabilize=15, pull=12):
        """Bruce furnace recipe: push/stabilize at 800 C, then (minutes,
        'ramp'|'hold', ambient) intervals, then pull at 800 C. Ramps are
        discretized in <= 50 C steps."""
        self.add(f'# PDF p.{page}: {name}; ramps discretized in <=50 C intervals.')
        segments = [(push, 800, 'nitrogen'), (stabilize, 800, 'nitrogen')]
        temp = 800
        for duration, kind, ambient in intervals:
            if duration <= 0: continue
            if kind == 'hold':
                segments.append((duration, temp, ambient)); continue
            end = target if kind == 'up' else 800
            count = max(1, math.ceil(abs(end - temp) / 50))
            segments.extend((duration / count, temp + (end - temp) * (i + .5) / count, ambient)
                            for i in range(count))
            temp = end
        segments.append((pull, 800, 'nitrogen'))
        previous = None
        for duration, temp, ambient in segments:
            continuation = previous is not None and not (
                ambient in ('dryo2', 'weto2') and previous == 'nitrogen')
            self.add(f'diffuse time={duration:.10g} temp={temp:.10g} {ambient}' +
                     (' continue' if continuation else ''))
            previous = ambient
        self.recipes.append({'name': name, 'pdf_page': page, 'segments': segments})

    def pad_oxide(self, name, page):
        # PDF p.23 recipe 250: ramp 20 min O2, soak 40 min O2, anneal 5 min N2,
        # ramp-down 40 min N2 (1000 C), push 12 / stabilize 15 / pull 12.
        # The 40 min soak grows 414 A here; the soak is calibrated to 500 A.
        self.recipe(name, page, 1000, [(20, 'up', 'dryo2'), (self.pad_soak, 'hold', 'dryo2'),
                                       (5, 'hold', 'nitrogen'), (40, 'down', 'nitrogen')])

    def recessed_oxide(self, name):
        # PDF p.36-39: thermal oxide, 3650 A target (Athena, p.39), "OX04 Tube 1".
        # The recipe is not given; assume the 0.5 um process field-oxide form
        # (SubCmos2014 p.68, recipe 406, 1100 C wet) with the soak calibrated
        # to 3650 A on bare silicon.
        self.recipe(name, 37, 1100, [(30, 'up', 'nitrogen'), (5, 'hold', 'dryo2'),
                                     (self.rox_soak, 'hold', 'weto2'),
                                     (5, 'hold', 'nitrogen'), (55, 'down', 'nitrogen')])

    def build(self):
        sign = 1 if self.device == 'nmos' else -1
        for x, spacing, tag in [(0, .1, 'left'), (.3, .05, ''), (1.5, .05, ''), (1.9, .04, ''),
                                (3, .08, ''), (4.4, .02, ''), (4.75, .01, ''), (5, .02, ''),
                                (5.25, .01, ''), (5.6, .02, ''), (7, .08, ''), (8.5, .04, ''),
                                (RIGHT, .1, 'right')]:
            self.add(f'line x loc={x:g} spac={spacing:g}' + (f' tag={tag}' if tag else ''))
        for y, spacing, tag in [(0, .005, 'top'), (.05, .01, ''), (.2, .02, ''), (.5, .05, ''),
                                (1, .15, ''), (3, .4, ''), (6, .8, 'bottom')]:
            self.add(f'line y loc={y:g} spac={spacing:g}' + (f' tag={tag}' if tag else ''))
        self.add('region silicon xlo=left xhi=right ylo=top yhi=bottom',
                 'boundary exposed xlo=left xhi=right ylo=top yhi=top',
                 'boundary backside xlo=left xhi=right ylo=bottom yhi=bottom',
                 '# PDF p.6/20: p-type 10 ohm-cm (100) substrate.',
                 'initialize boron conc=1.3e15 orient=100')
        self.checkpoint('initial')
        self.pad_oxide('500 A pad oxide, recipe 250', 23)
        self.checkpoint('pad')
        self.add('# PDF p.24: LPCVD Si3N4 1500 A.', 'deposit nitride thick=0.15 divisions=5')
        for lo, hi in FIELD: self.cut('nitride', lo, hi)
        self.checkpoint('active_mask')
        self.recessed_oxide('first recessed oxide, 3650 A')
        self.checkpoint('first_oxide')
        # PDF p.37: strip the first oxide (wet). A vertical timed strip keeps
        # the pad oxide under the nitride; 10% overetch.
        self.add('etch oxide dry thick=0.40')
        self.recessed_oxide('second recessed oxide, 3650 A')
        self.checkpoint('recessed_oxide')
        # PDF p.20 steps 10-14: well implants through nitride/pad oxide (active)
        # and the recessed oxide (field). Each cut is a single-well device.
        if self.device == 'pmos': self.implant('phosphorus', 3e13, 170)
        else: self.implant('boron', 8e13, 80)
        self.add('# PDF p.40: hot phosphoric acid nitride strip.', 'etch nitride all')
        # PDF p.42-43: recipe 11, 1100 C N2, ramp 30, soak 360, anneal 5, ramp-down 60.
        self.recipe('well drive, recipe 11', 43, 1100,
                    [(30, 'up', 'nitrogen'), (360, 'hold', 'nitrogen'),
                     (5, 'hold', 'nitrogen'), (60, 'down', 'nitrogen')])
        self.checkpoint('wells')
        # PDF p.50-53: Vt adjust through the 500 A pad oxide. p.53 and p.48
        # give PMOS 3.5e12 P31 60 keV; the p.20 step list (1.75e12 B11) is
        # inconsistent with the p+ poly PMOS and with RIT's own simulation.
        if self.device == 'nmos': self.implant('boron', 3e12, 30)
        else: self.implant('phosphorus', 3.5e12, 60)
        self.checkpoint('threshold')
        # PDF p.55: 10:1 BOE 45 s clears the pad oxide; timed with overetch.
        self.add('etch oxide dry thick=0.06')
        # PDF p.56-58: recipe 213, 900 C: ramp 30 min O2, 30 min N2O, 30 min O2,
        # ramp-down 30 min N2; target 100 A. SUPREM-IV has no N2O ambient.
        self.recipe('100 A N2O gate oxide, recipe 213', 58, 900,
                    [(30, 'up', 'dryo2'), (30, 'hold', self.n2o_as),
                     (self.gate_soak, 'hold', 'dryo2'), (30, 'down', 'nitrogen')],
                    stabilize=30, pull=15)
        self.checkpoint('gate_oxide')
        # Geometric poly re-oxidation (the vertical oxidation model does not
        # oxidize poly): each exposed poly face loses 0.44 x the poly oxide,
        # so the poly is deposited and etched to its post-re-ox size, and the
        # sidewall oxide is formed by conformal deposition plus vertical
        # etch-back (as the nitride spacers are).  PDF p.4: Lpoly 0.40 um after
        # etch -> 0.35 um after re-oxidation.  The poly oxide thickness (500 A
        # for recipe 250 at the <110> poly rate) comes from the compress-model
        # poly oxidation tests (suprem4/experiments/poly_reox_20261122).
        trim = 0.44 * self.poly_reox_ox if self.poly_reox == 'geometric' else 0.0
        gate = (round(GATE[0] + trim, 4), round(GATE[1] - trim, 4))
        self.add('# PDF p.63: LPCVD poly 4000 A, undoped; doped by the D/S implants.',
                 f'deposit poly thick={0.4 - trim:g} divisions=10 space=0.01')
        self.checkpoint('poly')
        self.add(f'etch poly left p1.x={gate[0]:g}', f'etch poly right p1.x={gate[1]:g}')
        self.checkpoint('gate')
        if self.poly_reox == 'geometric':
            self.add(f'# Poly re-ox (geometric): {1e4*self.poly_reox_ox:.0f} A poly oxide, '
                     f'{1e4*trim:.0f} A poly consumed per face (drawn {GATE[1]-GATE[0]:g} um '
                     f'-> {gate[1]-gate[0]:g} um).',
                     f'deposit oxide thick={self.poly_reox_ox:g} divisions=5 space=0.01',
                     f'etch oxide dry thick={self.poly_reox_ox:g}')
        self.pad_oxide('500 A poly re-oxidation, recipe 250', 74)
        self.checkpoint('reox')
        self.resist([DEVICE_ACTIVE])
        # PDF p.76/80: p-LDD B11 4e13 50 keV, n-LDD P31 4e13 60 keV.
        if self.device == 'nmos': self.implant('phosphorus', 4e13, 60)
        else: self.implant('boron', 4e13, 50)
        self.add('etch photoresist all')
        self.checkpoint('ldd')
        # PDF p.83-85: LPCVD nitride 3500 A, anisotropic etch, ~0.3 um spacers.
        self.add('deposit nitride thick=0.35 divisions=10 space=0.03',
                 'etch nitride dry thick=0.37')
        self.checkpoint('spacers')
        # PDF p.88-92: N+ D/S P31 4e15 60 keV, P+ D/S B11 4e15 50 keV. The
        # gate is open in its own D/S window: n+ poly NMOS, p+ poly PMOS.
        self.resist([DEVICE_ACTIVE], planar=True)
        self.implant('phosphorus' if self.device == 'nmos' else 'boron', 4e15,
                     60 if self.device == 'nmos' else 50)
        self.add('etch photoresist all')
        self.resist([BODY_ACTIVE], planar=True)
        self.implant('boron' if self.device == 'nmos' else 'phosphorus', 4e15,
                     50 if self.device == 'nmos' else 60)
        self.add('etch photoresist all')
        self.checkpoint('source_drain')
        # PDF p.95: recipe 284, 1000 C N2: ramp 20, soak 20, anneal 5, ramp-down 40.
        self.recipe('D/S anneal, recipe 284', 95, 1000,
                    [(20, 'up', 'nitrogen'), (20, 'hold', 'nitrogen'),
                     (5, 'hold', 'nitrogen'), (40, 'down', 'nitrogen')], pull=15)
        self.checkpoint('annealed')
        # PDF p.96-103 TiSi2 salicide is not modeled (no silicide material in
        # this simulator build); it does not change the channel.
        self.add('# PDF p.105: PECVD TEOS ILD 5000 A.',
                 'deposit oxide thick=0.5 divisions=10 space=0.05')
        self.resist(list(CONTACTS.values()), planar=True)
        self.add('etch oxide dry thick=1.5', 'etch photoresist all')
        self.checkpoint('contacts')
        # PDF p.113-114: W plugs + Al; modeled as aluminum filling the cuts.
        self.add('deposit aluminum thick=0.6 divisions=10 space=0.05')
        for lo, hi in METAL_GAPS: self.cut('aluminum', lo, hi, top=-6, bottom=1)
        self.checkpoint('metal1')
        self.add(f'structure outfile={self.device}.str', 'quit')
        return '\n'.join(self.lines) + '\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pad-soak', type=float, default=51,
                        help='dry soak (min) of recipe 250 (PDF: 40), calibrated to 500 A')
    parser.add_argument('--rox-soak', type=float, default=22,
                        help='wet soak (min) of each recessed oxidation, calibrated to 3650 A')
    parser.add_argument('--gate-soak', type=float, default=12,
                        help='final dry O2 soak (min) of recipe 213 (PDF: 30); 30/10/0 min '
                             'grow 144/94/63 A with N2O as nitrogen, 12 min meets 100 A')
    parser.add_argument('--n2o-as', choices=['nitrogen', 'dryo2'], default='nitrogen',
                        help='ambient standing in for the 30 min N2O interval')
    parser.add_argument('--poly-reox', choices=['none', 'geometric'], default='none',
                        help='poly re-oxidation: none (vertical model, poly not oxidized) or '
                             'geometric (poly consumption and sidewall oxide from the poly oxide thickness)')
    parser.add_argument('--poly-reox-ox', type=float, default=0.05,
                        help='oxide grown on poly by the re-oxidation (um), for --poly-reox geometric')
    parser.add_argument('--output', type=Path, default=HERE)
    args = parser.parse_args()
    output = args.output.resolve(); output.mkdir(parents=True, exist_ok=True)
    metadata = {'source': 'rit_process/AdvCmos2014_0p25um.pdf', 'gate_length_um': GATE[1] - GATE[0],
                'pad_soak_min': args.pad_soak, 'rox_soak_min': args.rox_soak, 'gate_soak_min': args.gate_soak,
                'n2o_as': args.n2o_as, 'poly_reox': args.poly_reox,
                'poly_reox_ox_um': args.poly_reox_ox, 'checkpoints': CHECKPOINTS,
                'layout': {'field': FIELD, 'gate': GATE, 'contacts': CONTACTS,
                           'metal_gaps': METAL_GAPS, 'probes': PROBES}, 'devices': {}}
    for device in ['nmos', 'pmos']:
        deck = Deck(device, args.rox_soak, args.gate_soak, args.n2o_as, args.pad_soak, args.poly_reox,
                    args.poly_reox_ox)
        (output / (device + '.in')).write_text(deck.build())
        metadata['devices'][device] = {'recipes': deck.recipes}
    (output / 'resist_top.dat').write_text(f'0 -2\n{RIGHT:g} -2\n')
    (output / 'process.json').write_text(json.dumps(metadata, indent=2) + '\n')


if __name__ == '__main__': main()
