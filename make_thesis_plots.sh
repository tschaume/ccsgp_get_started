#!/bin/bash

python -m pyana.examples.gp_stack --log --med LatestPatrickJieYi
python -m pyana.examples.gp_panel --log LatestPatrickJieYi
python -m pyana.examples.gp_rdiff --log --noxerr LatestPatrickJieYi
python -m pyana.examples.gp_rdiff --log --diffRel --noxerr LatestPatrickJieYi
python -m pyana.examples.gp_rdiff --log --diffRel LatestPatrickJieYi
python -m pyana.examples.gp_rdiff --log --divdNdy LatestPatrickJieYi
python -m pyana.examples.gp_ptspec --log
#open \
#  gp_stack/output/stackLatestPatrickJieYiInclMed.pdf \
#  gp_panel/output/panelLatestPatrickJieYi.pdf \
#  gp_rdiff/output/diffAbsLatestPatrickJieYiNoXErr.pdf \
#  gp_rdiff/output/diffRelLatestPatrickJieYiNoXErr.pdf \
#  gp_rdiff/output/enhanceLatestPatrickJieYi.pdf \
#  gp_rdiff/output/excessLatestPatrickJieYiDivdNdy.pdf \
#  gp_ptspec/output/ptspec.pdf \
#  gp_ptspec/output/meanPtLMR.pdf
