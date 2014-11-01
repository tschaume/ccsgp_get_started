#!/bin/bash

python -m pyana.examples.gp_stack --log --med LatestPatrickJieYi
python -m pyana.examples.gp_panel --log LatestPatrickJieYi
python -m pyana.examples.gp_rdiff --log --noxerr LatestPatrickJieYi
python -m pyana.examples.gp_rdiff --log --diffRel --noxerr LatestPatrickJieYi
python -m pyana.examples.gp_rdiff --log --diffRel LatestPatrickJieYi
python -m pyana.examples.gp_rdiff --log --divdNdy LatestPatrickJieYi
python -m pyana.examples.gp_ptspec --log
open \
  pyanaDir/examples/gp_stack/output/stackLatestPatrickJieYiInclMed.pdf \
  pyanaDir/examples/gp_panel/output/panelLatestPatrickJieYi.pdf \
  pyanaDir/examples/gp_rdiff/output/diffAbsLatestPatrickJieYiNoXErr.pdf \
  pyanaDir/examples/gp_rdiff/output/diffRelLatestPatrickJieYiNoXErr.pdf \
  pyanaDir/examples/gp_rdiff/output/enhanceLatestPatrickJieYi.pdf \
  pyanaDir/examples/gp_rdiff/output/excessLatestPatrickJieYiDivdNdy.pdf \
  pyanaDir/examples/gp_ptspec/output/ptspec.pdf
