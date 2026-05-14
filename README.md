# CosmoTadpole2
CosmoTadpole v0.2

Code associated with the paper Climatic controls on soil production, transport and chemical erosion: Insights from modelling topography, soils and cosmogenic nuclides at Little Lake, Oregon (https://onlinelibrary.wiley.com/doi/pdf/10.1002/esp.70197).

Here, we implemented some climate-sensitive earth surface processes into a landscape evolution model with an integrated 10Be tracer. In the paper, we compared modeled inferred denudation rates from 10Be to field-derived data from a paleo-lake core.

Other relevant (and great) papers on the cores at Little Lake:
https://pubs.geoscienceworld.org/gsa/gsabulletin/article-abstract/129/5-6/715/207851/Late-Quaternary-climatic-controls-on-erosion-rates
https://pubs.geoscienceworld.org/gsa/geology/article-abstract/47/7/613/570316/The-interplay-between-physical-and-chemical

This is MATLAB-based model, but may be able to run on GNU Octave if you can get the MEX C files compiled. We include Windows and MacOSX compiled binaries for the C code. The MEX files do compile and run on Linux. The C code is from Taylor Perron's Tadpole landscape evolution model. This model is a highly modified version of it.

To run the model from the paper, 

a) load the g_p_little_lake_steady_state.mat file,
b) [~, output] = Tadpole(g,p).

If you want to run the model on your own topography, send me an email at mreed5 [at] mail.wvu.edu / miles [at] hillslope.org, and we can get it set up.

