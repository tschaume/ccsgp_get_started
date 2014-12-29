#!/bin/bash

python -m pyana.examples.gp_rapp
python -m pyana.examples.gp_stack --med LatestPatrickJieYi
python -m pyana.examples.gp_panel LatestPatrickJieYi
python -m pyana.examples.gp_rdiff --noxerr LatestPatrickJieYi
python -m pyana.examples.gp_rdiff --diffRel --noxerr LatestPatrickJieYi
python -m pyana.examples.gp_rdiff --diffRel LatestPatrickJieYi
python -m pyana.examples.gp_rdiff --divdNdy LatestPatrickJieYi
python -m pyana.examples.gp_ptspec
open \
  pyanaDir/examples/gp_rapp/output/rapp_overview_panel.pdf \
  pyanaDir/examples/gp_stack/output/stackLatestPatrickJieYiInclMed.pdf \
  pyanaDir/examples/gp_panel/output/panelLatestPatrickJieYi.pdf \
  pyanaDir/examples/gp_rdiff/output/diffAbsLatestPatrickJieYiNoXErr.pdf \
  pyanaDir/examples/gp_rdiff/output/diffRelLatestPatrickJieYiNoXErr.pdf \
  pyanaDir/examples/gp_rdiff/output/enhanceLatestPatrickJieYi.pdf \
  pyanaDir/examples/gp_rdiff/output/excessLatestPatrickJieYiDivdNdy.pdf \
  pyanaDir/examples/gp_ptspec/output/ptspec.pdf
