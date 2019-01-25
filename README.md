# sugar-sweet
Exploring the role of storage compounds with Karthik and Stephie

last edit: Jen, 2019 Jan 25
commit: create table of contents for repo sharing


Repository contents:

1. glycogen_analysis.m
/ takes image data (.tiff) and applies image processing to recognize cells and extract property data. cells are linked through time via a particle tracking step, and then the linked data is assembled and saved in a data structure, D. To limit errors, there is a quality control step (well, actually 4) that produces the final data structure, D5, from which further analyses can be performed.

/ this script uses the following functions from this repository:
/		i. ParticleTrim_glycogen.m
/	   ii. Particle_Track_glycogen.m



2. dynamicOutlines_glycogen_phase.m
   dynamicOutlines_glycogen_cfp.m
   dynamicOutlines_glycogen_yfp.m
/  takes output data, D5, and raw images to visualize quality of particle tracking and intensity threshold-based strain identification

/  this script uses the following function from this repository:
/		i. buildDM_glycogen.m



3. quantifyGrowth.m
/  takes output data, D5, and manipulates it to quantify and visualize growth

/  this script uses the following functions from this repository:
/ 		i. buildDM_glycogen.m
/	   ii. calculateGrowthRate_glycogen.m



4. buildDM_glycogen.m
/  takes user-set variables to assemble an easy-to-use data matrix from organized but hairier D5 data structure



5. glycogen_analysis_scrap_pile.m
/  raw code used during testing and development of analysis, including:
/		i. tested methods to detect features from phase images
/	   ii. determination of intensity threshold to separate YFP and CFP




