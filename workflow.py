'''
------------------------------------------------------------------------------------------------------------------------
This workflow is for to produce the phylogeny for the Coryphoideae
------------------------------------------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------------------------------------
Author: Oscar Wrisberg
Date: 10/11/2021
------------------------------------------------------------------------------------------------------------------------
'''

from gwf import Workflow
# import math
# import glob

gwf = Workflow()

########################################################################################################################
################################################---- Secapr Quality Check Raw----#######################################
########################################################################################################################
# def secapr_quality_check_raw(ref, output,path, done):
#     """Quality checking raw data"""
#     inputs = [path+ref]
#     outputs = [path + output, path+done]
#     options = {'cores': 1, 'memory': "10g", 'walltime': "03:00:00", 'account':"Coryphoideae"}


#     spec = """
#     source activate secapr_env

#     secapr quality_check --input {input} --output {output}
    
#     touch {done}
#     """.format(input = path + ref, output = path+output, done = path+done)

#     return (inputs, outputs, options, spec)


# ########################################################################################################################
# ################################################---- Trimmomatic ----###################################################
# ########################################################################################################################
# def trimmomatic(ref, output, path, done):
#     """Trimming raw data"""
#     inputs = [path+ref]
#     outputs = [path + output, path+done]
#     options = {'cores': 1, 'memory': "10g", 'walltime': "00:30:00", 'account':"Coryphoideae"}

#     spec = """
#     source activate trimmomatic_env

#     trimmomatic PE -threads 16 -phred33 {input}_R1.fastq {input}_R2.fastq -baseout {output}.fastq\
#     ILLUMINACLIP:/home/owrisberg/miniconda3/pkgs/trimmomatic-0.39-1/share/trimmomatic-0.39-1/adapters/TruSeq3-PE-2.fa:2:30:10:1:true\
#     LEADING:3\
#     TRAILING:3\
#     MAXINFO:40:0.8\
#     MINLEN:36\
#     2>> stderr_trim_loop_output.txt
    
#     touch {done}
#     """.format(input = path + ref, output = path+output, done = path+done)

#     return (inputs, outputs, options, spec)

########################################################################################################################
################################################---- Hybpiper ----######################################################
########################################################################################################################
def hybpiper(species, p1, p2, un, path_out, path_in, done):
    """Hybpiper."""
    inputs = [path_in + species +p1, path_in + species + p2, path_in + species + un] # The files which the job will look for before it runs
    outputs = [path_out + "/" + species, done] # The files which will have to be created in order for the job to be "completed"
    options = {'cores': 1, 'memory': "20g", 'walltime': "100:00:00", 'account':"Coryphoideae"} #Slurm commands

    spec = """
    source activate base

    cd {out}
        
    /home/owrisberg/Coryphoideae/github_code/HybPiper/reads_first.py --cpu 1 --readfiles {p1} {p2} --unpaired {un} -b /home/owrisberg/Coryphoideae/target_sequence/PhyloPalms_loci_renamed_794-176_HEYcorrected.fasta --prefix {species} --bwa

    touch {done}
    """.format(species=species, p1 = path_in + species + p1,p2 = path_in + species + p2, un = path_in + species + un , out = path_out, done = done)


    return (inputs, outputs, options, spec)

########################################################################################################################
#############################################---- Paralogs ----#########################################################
########################################################################################################################
# def get_consSeqMask(ref, bam, bed, cov, out, path, done):
#     """Get consensus fasta."""
#     inputs = [path+bed]
#     outputs = [out, path+done]
#     options = {'cores': 4, 'memory': "20g", 'walltime': "18:00:00", 'account':"megaFauna"}

#     spec = """
#     source activate workflowMap
        
#     bcftools mpileup -R {bed} -C50 -f {ref} {bam} | bcftools call -c - | vcfutils.pl vcf2fq -d {cov1} -D {cov2} -l 5 | gzip > {out}
    
#     touch {done}
#     """.format(ref = path + ref, bam = path+bam, bed = path+bed, cov1=str(float(cov)/3), cov2=str(2*float(cov)),
#                out = out, done = path+done)

#     return (inputs, outputs, options, spec)

# ########################################################################################################################
# #############################################---- Intronerate ----######################################################
# ########################################################################################################################
# def get_consSeqMask(ref, bam, bed, cov, out, path, done):
#     """Get consensus fasta."""
#     inputs = [path+bed]
#     outputs = [out, path+done]
#     options = {'cores': 4, 'memory': "20g", 'walltime': "18:00:00", 'account':"megaFauna"}

#     spec = """
#     source activate workflowMap
        
#     bcftools mpileup -R {bed} -C50 -f {ref} {bam} | bcftools call -c - | vcfutils.pl vcf2fq -d {cov1} -D {cov2} -l 5 | gzip > {out}
    
#     touch {done}
#     """.format(ref = path + ref, bam = path+bam, bed = path+bed, cov1=str(float(cov)/3), cov2=str(2*float(cov)),
#                out = out, done = path+done)

#     return (inputs, outputs, options, spec)

# ########################################################################################################################
# #############################################---- Coverage ----#########################################################
# ########################################################################################################################
# def get_consSeqMask(ref, bam, bed, cov, out, path, done):
#     """Get consensus fasta."""
#     inputs = [path+bed]
#     outputs = [out, path+done]
#     options = {'cores': 4, 'memory': "20g", 'walltime': "18:00:00", 'account':"megaFauna"}

#     spec = """
#     source activate workflowMap
        
#     bcftools mpileup -R {bed} -C50 -f {ref} {bam} | bcftools call -c - | vcfutils.pl vcf2fq -d {cov1} -D {cov2} -l 5 | gzip > {out}
    
#     touch {done}
#     """.format(ref = path + ref, bam = path+bam, bed = path+bed, cov1=str(float(cov)/3), cov2=str(2*float(cov)),
#                out = out, done = path+done)

#     return (inputs, outputs, options, spec)

########################################################################################################################
######################################################---- RUN ----#####################################################
########################################################################################################################
sp = ["1001","1002","1003","1004","1005","1006","1007","1008","1009","1010","1011","1012","1013","1014","1015","1016","1017","1018","1019","1020","1021","1022","1023","1024","1025","1026","1027","1028","1029","1030","1031","1032","1033","1034","1035","1036","1037","1038","1039","1040","1041","1042","1043","1044","1045","1046","1047","1048","1049","1050","1051","1052","1053","1054","1055","1056","1057","1058","1059","1060","1061","1062","1063","1064","1065","1066","1067","1068","1069","1070","1071","1072","1073","1074","1075","1076","1077","1078","1079","1080","1081","1082","1193","1194","1195","1197","1222","1224","1238","1239","1244","1245","1246","1247","1248","1249","1250","1251","1252","1256","1259","1266","1267","1268","1269","1273","1274","1275","1276","1344","1345","1346","1347","1348","1349","1350","1351","1352","1353","1354","1355","1362","1363","1364","1411","1412","1413","1414","1415","1416","1417","1418","1419","1421","1430","1431","1432","1433","1434","1435","1437","1438","1439","1440","1442","1445","1449","1459","1460","1461","1462","1463","1465","1469","1470","1475","1476","1477","1478","1488","1489","1492","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2018","2019","2028","2035","2037","2038","2039","2040","2041","2042","2043","2044","2045","2046","2047","2048","2049","2050","2055","2058","2059","2060","2061","2062","2063","2064","2065","2066","2067","2068","2069","3000","3002","3004","3006","3008","3010","3012","3014","3016","3018","3020","3022","3024","3026","3028","3030","3032","3034","3036","3038","3040","3042","3044","3046","3048","3050","3052","3054","3056","3058","3060","3062","3064","3066","3068","3070","3072","3074","3076","3078","3080","3082","3084","3086","3088","3090","3092","3094","3096","3098","3100","3102","3104","3106","3108","3110","3112","3114","3116","3118","3120","3122","3124","3126","3128","3130","3132","3134","3136","3138","3140","3142","3144","3146","3148","3150","3152","3154","3156","3158","3160","3162","3164","3166","3168","3170","3172","3174","3176","3178","3180","3182","3184","3186","3188","3190","3192","3194","3196","3198","3200","3202","3204","3206","3208","3210","3212","3214","3216","3218","3220","3222","3224","3226","3228","3230","3232","3234","3236","3238","3240","3242","3244","3246","3248","3250","3252","3254","3256","3258","3260","3262","3264","3266","3268","3270","3272","3274","3276","3278","3280","3282","3284","3286","3288","3290","3292","3294","3296","3298","3300","3302","3304","3306","3308","3310","3312","3314","3316","3318","3320","3322","3324","3326","3328","3330","3332","3334","3336","3338","3340","3342","3344","3346","3348","3350","3352","3354","3356","3358","3360","3362","3364","3366","3368","3370","3372","3374","3376","3378","3380","3382","3384","3386","3388","3390","3392","3394","3396","3398","3400","3402","3404","3406","3408","3410","3412","3414","3416","3418","3420","3422","3424","3426","3428","3430","3432","3434","3436","3438","3440","3442","3444",]


for i in range(len(sp)):
    #### Running Hybpiper
    gwf.target_from_template('Hybpiper_'+sp[i], hybpiper(species = sp[i],
                                                        p1 = "_1P.fastq",
                                                        p2 = "_2P.fastq",
                                                        un = "_UN.fastq",
                                                        path_out= "/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/",
                                                        path_in = "/home/owrisberg/Coryphoideae/work_flow/02_trimmed/",
                                                        done = "/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/done/"+sp[i]))
                                                                      
    


    # #### make mask bed
    # gwf.target_from_template('makeMaskBed'+sp[i], make_mask_bed(ref = "ref/mask_35_90.fa",
    #                                                               bed = "ref/mask_35_90.bed",
    #                                                               path = "/home/juraj/megaFauna/data/"+sp[i]+"/",
    #                                                               done = "done/makeMaskBed"+sp[i]))
    # #### get diploid consensus sequence
    # gwf.target_from_template('intersect_beds'+sp[i], intersect_beds(bed1 = "ref/"+bed1[i],
    #                                                                 bed2 = "ref/mask_35_90.bed",
    #                                                                 out = "ref/intersect_35_90.bed",
    #                                                                 path = "/home/juraj/megaFauna/data/"+sp[i]+"/",
    #                                                                 done = "done/intersect_beds"+sp[i]))

    # cc = pd.read_table("/home/juraj/megaFauna/data/"+sp[i]+"/cov/covAll100kb"+sp[i]+".cov", header=None)
    # gwf.target_from_template('consSeq100kb_masked'+sp[i], get_consSeq(ref = "ref/"+ref[i],
    #                                                                 bam = "bam/"+bam[i],
    #                                                                 bed = "ref/intersect_35_90.bed",
    #                                                                 cov = str(cc.iloc[0,0]),
    #                                                                 out = "/home/juraj/megaFauna/results/"+sp[i]+"/psmc/diploid100kb_masked.fq.gz",
    #                                                                 path = "/home/juraj/megaFauna/data/"+sp[i]+"/",
    #                                                                 done= "done/consSeq100kb_masked"+sp[i]))


    # #### make psmc input
    # gwf.target_from_template('makePSMCInput100kb_masked'+sp[i], make_input(inp = "psmc/diploid100kb_masked.fq.gz",
    #                                                                      out = "psmc/diploid100kb_masked.psmcfa",
    #                                                                      path = "/home/juraj/megaFauna/results/"+sp[i]+"/",
    #                                                                      done = "done/makePSMCInput100kb_masked"+sp[i]))

    # #### run psmc
    # gwf.target_from_template('runPSMC100kb_masked'+sp[i], run_psmc(inp = "psmc/diploid100kb_masked.psmcfa",
    #                                                              out = "psmc/diploid100kb_masked.psmc",
    #                                                              path = "/home/juraj/megaFauna/results/"+sp[i]+"/",
    #                                                              done = "done/runPSMC100kb_masked"+sp[i]))
