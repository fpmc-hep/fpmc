This folder contains different photo-induced processes with intact protons in the final state which are embedded in the FPMC.

Though the name of the folder is excl_aaaa, it includes reactions AA->AZ, AA->ZZ, AA->WW, AA->Dijet, AA->HH, AA->AA (EFT), AA->AA (Spin0, Spin2 decay). These matrix elements were computed by Sylvain Fichet and Gero von Gersdorff and implemented by Cristian Baldenegro and Matthias Saimpert.

You can follow the structure of previous implementations to implement yours. You basically have to modify the excl_wrapper for the new process, the more general wrapper in the main directory, the Makefile on the top directory, the input parameters for the FPMC (Modify files in Examples directory, look for the files that have parameters like AAF0, AAM, AAF0, etc. These are related to these photo-induced processes.) and ***take care that the helicity_amplitudes.cpp and helicity_amplitudes.h files are identical among all the folders in excl_aaaa**. For some reason FPMC doesn't compile otherwise.

Contact Cristian Baldenegro crisx.baldenegro@gmail.com if you have questions about how to implement this or if there's a bug.
