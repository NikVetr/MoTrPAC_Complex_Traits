# The impact of exercise on gene regulation in association with complex trait genetics

This is the GitHub repository associated with the publication: [insert publication here]

It primarily hosts scripts written in the R programming language that were used to generate all of the paper figures. The names of these scripts are given in figure order, i.e. fig*_. In several cases, they assume other software is installed or otherwise in your system PATH, including various R packages, [plink2](https://www.cog-genomics.org/plink/2.0/), [GCTA](https://yanglab.westlake.edu.cn/software/gcta/#Overview), [MESC](https://github.com/douglasyao/mesc), [LDSC](https://github.com/bulik/ldsc), and [Stan](https://mc-stan.org/cmdstanr/). In other cases, it uses datasets pulled from other sources -- I've rehosted these in the /data/ folder where possible, noting the original data source if you'd prefer to acquire them directly.

Additionally, this GitHub plays host to several supplementary figures, datasets, MCMC output, and other files generated by these scripts.

In the interests of reproducibility, I've endeavored that, assuming all dependencies are properly installed, these scripts function as end-to-end pipelines. That said, just because they ran on my machine (macOS 10.15.7) does not mean they'll run seamlessly on yours, at least if you're reading this before I've had the chance to dockerize everything. Please contact me at nikgvetr@stanford.edu if you encounter difficulties and I'll do me best to help.
