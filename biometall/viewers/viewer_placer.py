# **************************************************************************
# *
# * Authors:     Eduardo Rivas Tortuero
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

"""
viewer_placer.py  ?  ViewerMetalPlacer
Viewer for ProtMetalPlacer (Path 1).

Multi-metal support: handles the new 'placements' list format in
placer_result.json.  Each metal type is coloured differently.

"""

import json
import os
import shutil
from subprocess import Popen
from pwem.protocols import EMProtocol

import pyworkflow.viewer as pwviewer
import pyworkflow.protocol.params as params
from pwem.viewers.viewer_chimera import ChimeraView
from pwchem.viewers.viewers_data import PyMolViewer
from pwem.objects import AtomStruct
from biometall.protocols.protocol_metal_placer import ProtMetalPlacer

METAL_COLORS_CX = {
    'ZN': '#4472C4', 'CA': '#70AD47', 'MG': '#00B050', 'FE': '#FF6600',
    'CU': '#C55A11', 'MN': '#7030A0', 'NI': '#00B0F0',
}

class ViewerMetalPlacer(pwviewer.ProtocolViewer, PyMolViewer):
    _label        = 'Viewer Metal Placer (homology + MetalKB)'
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets      = [ProtMetalPlacer]

    def __init__(self, **kwargs):
        pwviewer.ProtocolViewer.__init__(self, **kwargs)

    def _defineParams(self, form):
        form.addSection(label='Metal Placement Results')

        group = form.addGroup('Placement summary')
        group.addParam('showSummary', params.LabelParam,
                       label='Show result summary',
                       help='Display homologue PDB, sequence identity, '
                            'identified metal(s), and ion coordinates.')

        group2 = form.addGroup('3D Visualization')
        group2.addParam('displaySoftware', params.EnumParam,
                        choices=['PyMOL', 'ChimeraX'], default=0,
                        label='Visualize pseudo-holo structure with:',
                        help='Open pseudo_holo.pdb with the selected viewer.\n'
                             'Each metal type is shown in its own colour.\n'
                             'Coordinating residues (within 3.5 Å of any ion) '
                             'are highlighted automatically.')
        group2.addParam('displayStructure', params.LabelParam,
                        label='Open pseudo-holo structure',
                        help='Launch the viewer.')

    def _getVisualizeDict(self):
        return {
            'showSummary':      self._showSummary,
            'displayStructure': self._showStructure,
        }


    def _showSummary(self, paramName=None):
        resultJson = self.protocol._getExtraPath('placer_result.json')
        if not os.path.exists(resultJson):
            self.showInfo('Protocol not yet executed or no result was produced.\n\n'
                          'Possible reasons:\n'
                          '  - No homologue with >= 30 % identity found in RCSB\n'
                          '  - MetalKB found no binding site above the energy threshold')
            return []

        with open(resultJson) as fh:
            r = json.load(fh)

        homolog  = r.get('homolog_pdb') or '-'
        identity = r.get('identity')    or 0.0

        placements = r.get('placements') or []

        # Build per-metal block
        metal_blocks = ''
        for p in placements:
            metal          = p.get('metal', '?')
            n_placed       = p.get('count_placed', 0)
            n_expected     = p.get('count_expected', 1)
            positions      = p.get('positions', [])
            coord_residues = p.get('coord_residues') or {}

            warn = ''
            if n_placed < n_expected:
                warn = (' [WARNING: only %d/%d found]' % (n_placed, n_expected))

            pos_lines = ''
            for i, pos in enumerate(positions, start=1):
                pos_lines += '\n    Site %d:  X=%.3f  Y=%.3f  Z=%.3f A' % (
                    i, pos[0], pos[1], pos[2])
            if not pos_lines:
                pos_lines = '\n    -'

            coord_lines = ''
            for site_idx, residues in sorted(coord_residues.items(),
                                             key=lambda x: int(x[0])):
                res_str = ', '.join('%s%s' % (rn, rnum) for rn, rnum in residues)
                coord_lines += '\n    Site %s: %s' % (site_idx, res_str)
            if not coord_lines:
                coord_lines = '\n    -'

            metal_blocks += (
                '\n'
                '  Metal ion:            %s (%d/%d placed)%s\n'
                '  Ion positions:        %s\n'
                '  Coordinating res.:    %s\n'
            ) % (metal, n_placed, n_expected, warn, pos_lines, coord_lines)

        if not metal_blocks:
            metal_blocks = '\n  -  (no metals placed)'

        msg = (
            'Homologue PDB:      %s\n'
            'Sequence identity:  %.1f %%\n'
            'Source:             Experimental (X-ray / cryo-EM)\n'
            '%s'
        ) % (homolog, identity * 100, metal_blocks)

        self.showInfo(msg)
        return []

    def getAtomStruct(self):
      obj = self.protocol
      # If the input is a protocol (Analyze results was used), extract the AtomStruct obj
      if issubclass(type(obj), EMProtocol):
        for output in self.protocol.iterOutputAttributes(outputClass=AtomStruct):
          obj = output[1]
      return obj

    def _showStructure(self, paramName=None):
        pdbPath = self.getAtomStruct().getFileName()
        resultJson = self.protocol._getExtraPath('placer_result.json')

        if not os.path.exists(pdbPath):
            self.showInfo('No pseudo-holo structure found.\nRun the protocol first.')
            return []

        placements = []
        if os.path.exists(resultJson):
            with open(resultJson) as fh:
                r = json.load(fh)
            placements = [p for p in (r.get('placements') or [])
                          if p.get('count_placed', 0) > 0]

        if not placements:
            self.showInfo('No metals were placed ? nothing to visualise.')
            return []

        if self.displaySoftware.get() == 0:
            pymolV = PyMolViewer(project=self.getProject())
            return pymolV._visualize(pdbPath)
        else:
            return self._openChimeraX(pdbPath, placements)

    def _openChimeraX(self, pdbPath, placements):
        cxcPath = self.protocol._getExtraPath('viewer_chimerax.cxc')
        absPdb = os.path.abspath(pdbPath)

        lines = [
            'open "%s"' % absPdb,
            'hide atoms',
            'show cartoons',
            'color gray cartoons',
        ]

        for p in placements:
            metal = p['metal']
            color = METAL_COLORS_CX.get(metal.upper(), '#FFFF00')
            ionSpec = ':%s' % metal

            lines += [
                'show %s' % ionSpec,
                'style %s sphere' % ionSpec,
                'color %s %s' % (ionSpec, color),
            ]

        lines.append('view')

        with open(cxcPath, 'w') as fh:
            fh.write('\n'.join(lines) + '\n')

        return [ChimeraView(cxcPath)]