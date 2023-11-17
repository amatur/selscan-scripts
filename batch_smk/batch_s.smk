# reps = range(4)
# reps2 = [x+1 for x in reps]
reps = [1,2,3,4]
Ls=[1000]

#Ls=[1500,2000]
NUM_THREADS = 8

rule all:
    input: 
        expand("s{r}.{l}.vcf",l=Ls, r=reps),
        expand("s{r}.{l}.map",l=Ls, r=reps),
        expand("s{r}.{l}.vcf.rehh",r=reps, l=Ls),
        expand("outihs/s{r}.{l}.ihs.out",r=reps, l=Ls),
        
        #expand("outihs/s{r}.{l}.ihs.alt.out",r=reps, l=Ls),
        # expand("outihs/out.{effective_pop}_{l}_{r}.hapbin",r=reps, l=Ls, effective_pop=Ne),
        expand("outihs/s{r}.{l}.selbin",r=reps, l=Ls),

                
# for r in reps:
#     for l in Ls:
#         rule:
#             input: f"sim_params.yaml", 
#             output: 
#                 hap=expand("simulations_sweep/out.{effective_pop}_{leng}_{rep}.hap",  leng=l, effective_pop=Ne, rep=r),
#                 #map=expand("simulations/out.{effective_pop}_{leng}_{rep}.map", leng=l, effective_pop=Ne, rep=r),
#                 vcf=expand("simulations_sweep/out.{effective_pop}_{leng}_{rep}.vcf",leng=l, effective_pop=Ne, rep=r)
#             params: effective_pop = Ne, rep=r, leng=l
#             benchmark: 
#                 f"benchmarks/sim.{Ne}_{l}_{r}",
#             shell: f"python msprime_sel.py {{params.rep}} {{params.leng}} {{params.effective_pop}}  > {{output.hap}}"


rule vcf2map:
    input: "s{r}.{l}.vcf"
    output: "s{r}.{l}.map"
    shell:
        """
        grep -v '#' {input} | cut -f 2 | awk '{{ print 1, NR-1, $1/10000000.0, $1 }}' > {output}
        """

rule rehh:
    input: "s{r}.{l}.vcf"
    output: 
        out="s{r}.{l}.vcf.rehh",
        mem="mem/s{r}.{l}.rehh"
    # resources:
    #     mem_mb = 1000
    benchmark:
        "benchmarks/rehh.{r}.{l}",
    params:
        NUM_THREADS=NUM_THREADS
    resources:
        mem_mb=80000,
        tmpdir="~/s/tmp",
        disk_mb=10000
    threads:
        NUM_THREADS
    shell:
        """
        module load r/4.3.2; 
        set +e
        /usr/bin/time  -f "%M %e" --output-file={output.mem} Rscript ~/usr/local/bin/rehh.R {input} {params.NUM_THREADS}
        exitcode=$?
        if [ $exitcode -eq 1200 ]
        then
            exit 1
        else
            exit 0
        fi
        """

rule selscan_no_alt:
    input: 
        vcf="s{r}.{l}.vcf",
        map="s{r}.{l}.map",

    output: 
        out="outihs/s{r}.{l}.ihs.out",
        mem="mem/s{r}.{l}.selscan"
    
    # resources:
    #     mem_mb = 1000
    benchmark:
        "benchmarks/selscan.{r}.{l}",
    params:
        NUM_THREADS=NUM_THREADS
    resources:
        mem_mb=30000,
        tmpdir="~/s/tmp",
        disk_mb=10000

    threads:
        NUM_THREADS
    shell:
        """
        /usr/bin/time  -f "%M %e" --output-file={output.mem} selscan --vcf {input.vcf} --ihs --map {input.map} --trunc-ok --threads {params.NUM_THREADS} --out outihs/s{wildcards.r}.{wildcards.l} 
        """

rule selscan:
    input: "s{r}.{l}.vcf"
    output: "outihs/s{r}.{l}.ihs.alt.out"
    
    # resources:
    #     mem_mb = 1000
    benchmark:
        "benchmarks/selscan.{r}.{l}",
    params:
        NUM_THREADS=NUM_THREADS
    resources:
        mem_mb=30000,
        tmpdir="~/s/tmp",
        disk_mb=10000

    threads:
        NUM_THREADS
    shell:
        """
        /usr/bin/time  -f "%M %e" --output-file=mem/s{wildcards.r}.{wildcards.l}.selscan selscan --vcf {input} --ihs --pmap --trunc-ok --threads {params.NUM_THREADS} --out outihs/s{wildcards.r}.{wildcards.l} --alt
        """


# rule hapbin:
#     input: 
#         hap="simulations_sweep/out.{effective_pop}_{leng}_{rep}.hap",
#         map="simulations_sweep/out.{effective_pop}_{leng}_{rep}.map"
#     output: 
#         "outihs/out.{effective_pop}_{leng}_{rep}.hapbin"
#     benchmark:
#         "benchmarks/hapbin.{effective_pop}_{leng}_{rep}",
#     params:
#         NUM_THREADS=NUM_THREADS
#     shell:
#         """
#         export OMP_NUM_THREADS={params.NUM_THREADS}; /usr/bin/time  -f "%M %e" --output-file=mem_{wildcards.effective_pop}_{wildcards.leng}_{wildcards.rep}_{params.NUM_THREADS}.hapbin ihsbin --hap  {input.hap} --map {input.map} --out {output} --minmaf 0.05 --cutoff 0.05 --scale 20000000
#         """


rule selbin:
    input: 
        "s{r}.{l}.vcf"
    output: 
        out="outihs/s{r}.{l}.selbin",
        mem="mem/s{r}.{l}.selbin"
    benchmark:
        "benchmarks/selbin.{r}.{l}",
    resources:
        mem_mb=14000,
        tmpdir="~/s/tmp",
        disk_mb=10000
    params:
        NUM_THREADS=NUM_THREADS
    threads:
        NUM_THREADS
    shell:
        """
        /usr/bin/time  -f "%M %e" --output-file={output.mem} selbin ehh --all --vcf {input} --threads {params.NUM_THREADS}  --out {output.out} --alt
        """



# rule plot:

# $ cat mem_selbin_1000_10000000_1
# 7604 23.54
# (base) [aur1111@submit02 run_devel]$ cat mem_hapbin_1000_10000000_1
# 140496 5.65
# (base) [aur1111@submit02 run_devel]$ cat mem_selscan_1000_10000000_1


# wc -l

# mem_selbin mem_hapbin mem_selscan

# time_selbin time_hapbin time_selscan


# .selbin  | wc -l 
# .selscan | wc -l
# tail -n +2 {aa}.hapbin | wc -l

# nthread 1




#prefilter
#bcftools view -i 'AF>0.05'

