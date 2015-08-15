# Cover letter

Dear Editors,

I am submitting the attached manuscript for your consideration for publication in JCTC:

A simple method for automated equilibration detection in molecular simulations.
by John D. Chodera.

While the practice of discarding some initial data from molecular simulations (for both molecular dyanmics or Monte Carlo) to "equilibration" is essentially universally practiced in our field, there is yet no standard best practice for how this process is carried out.  While some previous published work has suggested schemes for selecting the equilibrated region of molecular simulations automatically, these methods are typically complex to implement, make incorrect assumptions about the statistics of the data, or are otherwise not universally accessible and applicable to researchers.

In this manuscript, I present an extremely simple and yet surprisingly useful automated approach to selecting the equilibrated region of molecular simulations.  In short, by selecting the region that maximizes the number of uncorrelated samples in the collected dataset, one achieves near-optimal performance, where optimality is demonstrated here by quantifying statistical performance on hundreds of independent replications of the same experiment.

In hopes of supporting reproducibility in science, I have also made all of the source code used to generate the data and produce the figures used in the manuscript freely available in the public source code repository GitHub [https://github.com/choderalab/automatic-equilibration-detection].  There is an accompanying "reproduce.sh" bash script that, when run, retrieves all of the exact versions of software dependencies required, runs all simulations, analyzes the data, and produces the figures.  It was run to produce the figures appearing in this manuscript, and others have already indicated that they have been able to run this script to reproduce the figures on their own machines.

The automatic equilibration detection scheme I describe is also made freely available in conveniently installable form through the Python 'pymbar' package [http://pymbar.org], which can be installed easily through the conda scientific Python package manager in a single line.

A version of this manuscript has been posted to the bioRxiv preprint server [http://biorxiv.org/content/early/2015/07/04/021659] and has already attracted substantial interest, having been downloaded over 278 times in the few weeks it has been posted online.

All authors (myself alone) consent to the submission of this manuscript.

All the best,

John D. Chodera
Assistant Member, Computational Biology Program
Sloan Kettering Institute
Memorial Sloan Kettering Cancer Center
