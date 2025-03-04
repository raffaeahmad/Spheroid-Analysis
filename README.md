# Spheroid-Analysis
Raffae N. Ahmad, Virginia Tech


  If using this script, please cite: "Characterization of Glioma Spheroid Viability and Metastatic Potential Following Monophasic and Biphasic Pulsed Electric Fields, Bioelectrochemistry, 2025"

This code was created to analyze binary masks of spheroids. Two masks should be made for the outgrowth region and spheroid body. The file naming structure is samplename.tif for original image, samplename_mask_outgrowth_dayX.tif for outgrowth mask and samplename_mask_body_dayX.tif for body masks. Output is in pixel units and must be converted to um and um^2 if needed. Pixel units can bekept if you decide to use relative area between original and post-treatment images. 
