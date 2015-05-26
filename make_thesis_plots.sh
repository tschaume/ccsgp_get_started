#!/bin/bash

python -m ccsgp_get_started.examples.gp_rapp
python -m ccsgp_get_started.examples.gp_stack --med LatestPatrickJieYi
python -m ccsgp_get_started.examples.gp_panel LatestPatrickJieYi
python -m ccsgp_get_started.examples.gp_rdiff --noxerr LatestPatrickJieYi
python -m ccsgp_get_started.examples.gp_rdiff --diffRel --noxerr LatestPatrickJieYi
python -m ccsgp_get_started.examples.gp_rdiff --diffRel LatestPatrickJieYi
python -m ccsgp_get_started.examples.gp_rdiff --divdNdy LatestPatrickJieYi
python -m ccsgp_get_started.examples.gp_ptspec
open \
  output/examples/gp_rapp/rapp_overview_panel.pdf \
  output/examples/gp_stack/stackLatestPatrickJieYiInclMed.pdf \
  output/examples/gp_panel/panelLatestPatrickJieYi.pdf \
  output/examples/gp_rdiff/diffAbsLatestPatrickJieYiNoXErr.pdf \
  output/examples/gp_rdiff/diffRelLatestPatrickJieYiNoXErr.pdf \
  output/examples/gp_rdiff/enhanceLatestPatrickJieYi.pdf \
  output/examples/gp_rdiff/excessLatestPatrickJieYiDivdNdy.pdf \
  output/examples/gp_ptspec/ptspec.pdf
