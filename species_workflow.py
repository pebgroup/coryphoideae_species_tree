'''
------------------------------------------------------------------------------------------------------------------------
This workflow is used to transform the raw sequence data into sequences ready for alignment.
Workflow is the following

1: Secapr quality check of the sequences
2: Trimming of the sequences using trimmomatic
3: Hybpipering the species in order to create the exons of the baits for each species
4: Checking for Paralogs
5: Running Intronerate again to get the Introns for each species
6: Trimming for Coverage of sequencing and joining the exon and intron of each species to create supercontigs

------------------------------------------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------------------------------------
Author: Oscar Wrisberg
Date: 10/11/2021
------------------------------------------------------------------------------------------------------------------------
'''

from os import O_SYNC, name
from gwf import Workflow
import os.path
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
    outputs = [path_out + species, done] # The files which will have to be created in order for the job to be "completed"
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

def paralogs(species,path_in, done):
    """Find Paralog genes and write them in the file called paralog.txt"""
    inputs = [path_in + species]
    outputs = [done]
    options = {'cores': 2, 'memory': "10g", 'walltime': "8:00:00", 'account':"Coryphoideae"}

    spec = """
    source activate base
    
    
    cd {path_in}
        
    python /home/owrisberg/Coryphoideae/github_code/HybPiper/paralog_investigator.py {sp} 2>> paralog.txt

    
    touch {done}


    """.format(sp = species, done = done, path_in = path_in)
    return (inputs, outputs, options, spec)

def no_paralogs(species, path_in, done, no_paralogs):
    """Wrapper script to continue pipeline when Hybpiper finds no paralogs"""
    inputs = [path_in + species]
    outputs = [done]
    options = {'cores': 2, 'memory': "10g", 'walltime': "0:05:00", 'account':"Coryphoideae"}

    spec = """

    touch {done}
    touch {np}

    """.format(done=done, np=no_paralogs)
    return(inputs, outputs, options, spec)

# ########################################################################################################################
# #############################################---- Intronerate ----######################################################
# ########################################################################################################################

def intronerate(species, path_in, done):
    """Intronerate the sequencec from hybpiper."""
    inputs = [path_in + species]
    outputs = [done]
    options = {'cores': 4, 'memory': "20g", 'walltime': "16:00:00", 'account':"Coryphoideae"}

    spec = """
    source activate base

    cd {path_in}

    python3 /home/owrisberg/Coryphoideae/github_code/HybPiper/intronerate.py --prefix {sp} &>> intronerate_out.txt
        
    
    touch {done}
    """.format(sp = species, done = done, path_in = path_in)

    return (inputs, outputs, options, spec)

# ########################################################################################################################
# #############################################---- Coverage ----#########################################################
# ########################################################################################################################
def coverage(species, path_in, path_out, done,all_bam,all_sorted_bam, all_sorted_bam_bai, bam, cov,fasta,fasta_amb,fasta_ann,fasta_bwt,fasta_pac,fasta_sa,trimmed_fasta,up_bam):
    """Calculating coverage of sequences."""
    inputs = [path_in+species, path_in+"done/Intronerate/"+species]
    outputs = [path_out+species+all_bam, path_out+species+all_sorted_bam, path_out+species+all_sorted_bam_bai, path_out+species+bam,
    path_out+species+cov, path_out+species+fasta, path_out+species+fasta_amb, path_out+species+fasta_ann, path_out+species+fasta_bwt,
    path_out+species+fasta_pac, path_out+species+fasta_sa, path_out+species+trimmed_fasta, path_out+species+up_bam,done] #ALL the output files
    options = {'cores': 4, 'memory': "20g", 'walltime': "08:00:00", 'account':"Coryphoideae"}

    spec = """
    source activate base
    
    cd {path_in}

    python3 /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/coverage.py {sp}
    
    touch {done}

    """.format(sp = species, done = done, path_in = path_in)

    return (inputs, outputs, options, spec)

########################################################################################################################
######################################################---- RUN ----#####################################################
########################################################################################################################
sp = ["1001","1002","1003","1004","1005","1006","1007","1008","1009","1010","1011","1012","1013","1014","1015","1016","1017","1018","1019","1020","1021","1022","1023","1024","1025","1026","1027","1028","1029","1030","1031","1032","1033","1034","1035","1036","1037","1038","1039","1040","1041","1042","1043","1044","1045","1046","1047","1048","1049","1050","1051","1052","1053","1054","1055","1056","1057","1058","1059","1060","1061","1062","1063","1064","1065","1066","1067","1068","1069","1070","1071","1072","1073","1074","1075","1076","1077","1078","1079","1080","1081","1082","1193","1194","1195","1197","1222","1224","1238","1239","1244","1245","1246","1247","1248","1249","1250","1251","1252","1256","1259","1266","1267","1268","1269","1273","1274","1275","1276","1344","1345","1346","1347","1348","1349","1350","1351","1352","1353","1354","1355","1362","1363","1364","1411","1412","1413","1414","1415","1416","1417","1418","1419","1421","1430","1431","1432","1433","1434","1435","1437","1438","1439","1440","1442","1445","1449","1459","1460","1461","1462","1463","1465","1469","1470","1475","1476","1477","1478","1488","1489","1492","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2018","2019","2028","2035","2037","2038","2039","2040","2041","2042","2043","2044","2045","2046","2047","2048","2049","2050","2055","2058","2059","2060","2061","2062","2063","2064","2065","2066","2067","2068","2069","3000","3002","3004","3006","3008","3010","3012","3014","3016","3018","3020","3022","3024","3026","3028","3030","3032","3034","3036","3038","3040","3042","3044","3046","3048","3050","3052","3054","3056","3058","3060","3062","3064","3066","3068","3070","3072","3074","3076","3078","3080","3082","3084","3086","3088","3090","3092","3094","3096","3098","3100","3102","3104","3106","3108","3110","3112","3114","3116","3118","3120","3122","3124","3126","3128","3130","3132","3134","3136","3138","3140","3142","3144","3146","3148","3150","3152","3154","3156","3158","3160","3162","3164","3166","3168","3170","3172","3174","3176","3178","3180","3182","3184","3186","3188","3190","3192","3194","3196","3198","3200","3202","3204","3206","3208","3210","3212","3214","3216","3218","3220","3222","3224","3226","3228","3230","3232","3234","3236","3238","3240","3242","3244","3246","3248","3250","3252","3254","3256","3258","3260","3262","3264","3266","3268","3270","3272","3274","3276","3278","3280","3282","3284","3286","3288","3290","3292","3294","3296","3298","3300","3302","3304","3306","3308","3310","3312","3314","3316","3318","3320","3322","3324","3326","3328","3330","3332","3334","3336","3338","3340","3342","3344","3346","3348","3350","3352","3354","3356","3358","3360","3362","3364","3366","3368","3370","3372","3374","3376","3378","3380","3382","3384","3386","3388","3390","3392","3394","3396","3398","3400","3402","3404","3406","3408","3410","3412","3414","3416","3418","3420","3422","3424","3426","3428","3430","3432","3434","3436","3438","3440","3442","3444"]

for i in range(len(sp)):
    #### Running Hybpiper
    gwf.target_from_template('Hybpiper_'+sp[i], hybpiper(species = sp[i],
                                                        p1 = "_1P.fastq",
                                                        p2 = "_2P.fastq",
                                                        un = "_UN.fastq",
                                                        path_out= "/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/",
                                                        path_in = "/home/owrisberg/Coryphoideae/work_flow/02_trimmed/",
                                                        done = "/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/done/Hybpiper/"+sp[i]))
                                                                      
    


    #### Paralogs
    if os.path.isfile("/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/"+sp[i]+"/genes_with_paralog_warnings.txt"):
        gwf.target_from_template('Paralogs_'+sp[i], paralogs(species = sp[i],
                                                            path_in = "/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/",
                                                            done = "/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/done/Paralogs/"+sp[i],
                                                            ))
    else:
        gwf.target_from_template('No_Paralogs_'+sp[i], no_paralogs(species = sp[i],
                                                                path_in = "/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/",
                                                                done = "/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/done/Paralogs/"+sp[i],
                                                                no_paralogs="/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/done/No_paralogs"+sp[i]))
     
    
    #### Getting introns
    gwf.target_from_template('Intronerate_'+sp[i], intronerate(species= sp[i],
                                                        path_in = "/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/",
                                                        done = "/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/done/Intronerate/"+sp[i]))


    #### Coverage
    gwf.target_from_template('Coverage_'+sp[i], coverage(species = sp[i],
                                                        path_in = "/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/",
                                                        all_bam = "_all.bam",
                                                        all_sorted_bam ="_all_sorted.bam",
                                                        all_sorted_bam_bai="_all_sorted.bam.bai",
                                                        bam =".bam",
                                                        cov=".cov",
                                                        fasta = ".fasta",
                                                        fasta_amb = ".fasta.amb",
                                                        fasta_ann = ".fasta.ann",
                                                        fasta_bwt = ".fasta.bwt",
                                                        fasta_pac = ".fasta.pac",
                                                        fasta_sa = ".fasta.sa",
                                                        trimmed_fasta = "_trimmed.fasta",
                                                        up_bam = "_up.bam",
                                                        path_out = "/home/owrisberg/Coryphoideae/work_flow/04_coverage/",
                                                        done = "/home/owrisberg/Coryphoideae/work_flow/04_coverage/done/Coverage/"+sp[i]))

