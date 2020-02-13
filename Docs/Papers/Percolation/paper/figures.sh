

PDFRESOLUTION=200
WIDTH=3000
HORIZONTALPADDING=10
VERTICALPADDING=10

#convert -density $PDFRESOLUTION figuresraw/principle.pdf -resize "$WIDTH"x Fig1.png
cp figuresraw/principle.png Fig1.png

montage figuresraw/abssize_nodes.png figuresraw/relsize_nodes.png -resize "$((WIDTH / 2))"x -tile 2x1 -geometry +"$HORIZONTALPADDING"+"$VERTICALPADDING" Fig2.png

convert figuresraw/fractaldimension.png -resize "$WIDTH"x Fig3.png

montage figuresraw/totalPop4183694_00056402_ecount850_radius8000.png figuresraw/totalPop2219780_36719597_vcount378_radius8000.png figuresraw/totalPop1474347_36891685_vcount595_radius50000.png figuresraw/totalPop1474347_36891685_euclPerf0_00507404608502498_radius54000.png -resize "$((WIDTH / 2))"x -tile 2x2 -geometry +"$HORIZONTALPADDING"+"$VERTICALPADDING" Fig4.png

convert figuresraw/full_effective_pareto.png -resize "$WIDTH"x Fig5.png

convert figuresraw/aggreg_morpho_pc1-emissions_targeted.png -resize "$WIDTH"x Fig6.png

convert figuresraw/aggreg_morpho_relemissions-relefficiency_colpc1_logscale_targeted.png -resize "$WIDTH"x Fig7.png
