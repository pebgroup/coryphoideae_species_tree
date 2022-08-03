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
################################################---- Fastqc quality check raw ----#######################################
########################################################################################################################
def fastqc_raw(species,path_in ,path_out, done,):
    """Quality checking using fastqc as this should work on individual species"""
    inputs = []
    outputs = [path_out+species+"_R1_fastqc.html",path_out+species+"_R2_fastqc.html", done]
    options = {'cores': 1, 'memory': "10g", 'walltime': "00:30:00", 'account':"Coryphoideae"}


    spec = """
    source activate secapr_env

    fastqc -o {output} {path_in}{species}_R1.fastq {path_in}{species}_R2.fastq
    
    echo touching {done}
    touch {done}

    """.format(path_in = path_in,species = species, output = path_out, done = done)

    return (inputs, outputs, options, spec)

########################################################################################################################
################################################---- Fastqc quality check trimmed ----#######################################
########################################################################################################################
def fastqc_trimmed(species,path_in ,path_out, done,):
    """Quality checking using fastqc as this should work on individual species"""
    inputs = [path_in+species+"_UN.fastq", path_in+species+"_1PU.fastq", path_in+species+"_2PU.fastq","/home/owrisberg/Coryphoideae/work_flow/02_trimmed/done/"+species]
    outputs = [path_out+species+"_1PU_fastqc.html", path_out+species+"_2PU_fastqc.html",path_out+species+"_UN_fastqc.html" ,done]
    options = {'cores': 1, 'memory': "10g", 'walltime': "00:30:00", 'account':"Coryphoideae"}


    spec = """
    source activate secapr_env

    fastqc -o {output} {path_in}{species}_1PU.fastq {path_in}{species}_2PU.fastq {path_in}{species}_UN.fastq
    
    touch {done}

    """.format(path_in = path_in,species = species, output = path_out, done = done)

    return (inputs, outputs, options, spec)



# ########################################################################################################################
# ################################################---- Secapr Quality Check Raw----#######################################
# ########################################################################################################################
# def secapr_quality_check_raw(path_in, path_out, done,):
#     """Quality checking
#     Cant be set to individual species because SecapR quality check requires the input to be a folder"""
#     inputs = []
#     outputs = [path_out, done]
#     options = {'cores': 5, 'memory': "20g", 'walltime': "120:00:00", 'account':"Coryphoideae"}


#     spec = """
#     source activate secapr_env

#     secapr quality_check --input {input} --output {output}
    
#     touch {done}

#     """.format(input = path_in, output = path_out, done = done)

#     return (inputs, outputs, options, spec)


# ########################################################################################################################
# ################################################---- Trimmomatic ----###################################################
# ########################################################################################################################
def trimmomatic(species, path_in, path_out, done):
    """Trimming raw data using trimmomatic with custom adapter.
    Afterwards combines paired and unpaired reads for forward and reverse reads respectively for each species 
    to enable post-trimming secapr quality_check for comparability before and after trimming """
    inputs = []
    outputs = [path_out+species+"_UN.fastq",path_out+"secapr_postrim/"+species+"_UN.fastq", done, path_out+species+"_1P.fastq", path_out+species+"_2P.fastq"]
    options = {'cores': 16, 'memory': "10g", 'walltime': "01:00:00", 'account':"Coryphoideae"}

    spec = """
    source activate trimmomatic_env

    trimmomatic PE -threads 16 -phred33 {input}_R1.fastq {input}_R2.fastq -baseout {output}.fastq\
    ILLUMINACLIP:/home/owrisberg/miniconda3/pkgs/trimmomatic-0.39-1/share/trimmomatic-0.39-1/adapters/TruSeq3-PE-2.fa:2:30:10:1:true\
    LEADING:3\
    TRAILING:3\
    MAXINFO:40:0.8\
    MINLEN:36\
    2>> stderr_trim_loop_output.txt

    echo combining {path_out}{species}_1P.fastq and {path_out}{species}_1U.fastq into {path_out}secapr_postrim/{species}_1PU.fastq 
    cat {path_out}{species}_1P.fastq {path_out}{species}_1U.fastq > {path_out}secapr_postrim/{species}_1PU.fastq 

    echo combining {path_out}{species}_2P.fastq and {path_out}{species}_2U.fastq into {path_out}secapr_postrim/{species}_2PU.fastq 
    cat {path_out}{species}_2P.fastq {path_out}{species}_2U.fastq > {path_out}secapr_postrim/{species}_2PU.fastq

    echo combining {path_out}{species}_1U.fastq {path_out}{species}_2U.fastq > {path_out}{species}_UN.fastq
    cat {path_out}{species}_1U.fastq {path_out}{species}_2U.fastq > {path_out}{species}_UN.fastq
    cp {path_out}{species}_UN.fastq {path_out}secapr_postrim/


    echo Removing {path_out}{species}_1U.fastq
    rm {path_out}{species}_1U.fastq

    echo Removing {path_out}{species}_2U.fastq
    rm {path_out}{species}_2U.fastq

    touch {done}

    """.format(input = path_in + species, output = path_out+species, done = done, species = species, path_out = path_out)

    return (inputs, outputs, options, spec)


########################################################################################################################
################################################---- Secapr Quality Check Trimmed----#######################################
########################################################################################################################
# def secapr_quality_check_trimmed(path_in, path_out, done,trim_done):
#     """Quality checking
#     Cant be set to individual species because SecapR quality check requires the input to be a folder"""
#     inputs = [trim_done+"1001",trim_done+"1002",trim_done+"1003",trim_done+"1004",trim_done+"1005",trim_done+"1006",trim_done+"1007",trim_done+"1008",trim_done+"1009",trim_done+"1010",trim_done+"1011",trim_done+"1012",trim_done+"1013",trim_done+"1014",trim_done+"1015",trim_done+"1016",trim_done+"1017",trim_done+"1018",trim_done+"1019",trim_done+"1020",trim_done+"1021",trim_done+"1022",trim_done+"1023",trim_done+"1024",trim_done+"1025",trim_done+"1026",trim_done+"1027",trim_done+"1028",trim_done+"1029",trim_done+"1030",trim_done+"1031",trim_done+"1032",trim_done+"1033",trim_done+"1034",trim_done+"1035",trim_done+"1036",trim_done+"1037",trim_done+"1038",trim_done+"1039",trim_done+"1040",trim_done+"1041",trim_done+"1042",trim_done+"1043",trim_done+"1044",trim_done+"1045",trim_done+"1046",trim_done+"1047",trim_done+"1048",trim_done+"1049",trim_done+"1050",trim_done+"1051",trim_done+"1052",trim_done+"1053",trim_done+"1054",trim_done+"1055",trim_done+"1056",trim_done+"1057",trim_done+"1058",trim_done+"1059",trim_done+"1060",trim_done+"1061",trim_done+"1062",trim_done+"1063",trim_done+"1064",trim_done+"1065",trim_done+"1066",trim_done+"1067",trim_done+"1068",trim_done+"1069",trim_done+"1070",trim_done+"1071",trim_done+"1072",trim_done+"1073",trim_done+"1074",trim_done+"1075",trim_done+"1076",trim_done+"1077",trim_done+"1078",trim_done+"1079",trim_done+"1080",trim_done+"1081",trim_done+"1082",trim_done+"1193",trim_done+"1194",trim_done+"1195",trim_done+"1197",trim_done+"1222",trim_done+"1224",trim_done+"1238",trim_done+"1239",trim_done+"1244",trim_done+"1245",trim_done+"1246",trim_done+"1247",trim_done+"1248",trim_done+"1249",trim_done+"1250",trim_done+"1251",trim_done+"1252",trim_done+"1256",trim_done+"1259",trim_done+"1266",trim_done+"1267",trim_done+"1268",trim_done+"1269",trim_done+"1273",trim_done+"1274",trim_done+"1275",trim_done+"1276",trim_done+"1344",trim_done+"1345",trim_done+"1346",trim_done+"1347",trim_done+"1348",trim_done+"1349",trim_done+"1350",trim_done+"1351",trim_done+"1352",trim_done+"1353",trim_done+"1354",trim_done+"1355",trim_done+"1362",trim_done+"1363",trim_done+"1364",trim_done+"1411",trim_done+"1412",trim_done+"1413",trim_done+"1414",trim_done+"1415",trim_done+"1416",trim_done+"1417",trim_done+"1418",trim_done+"1419",trim_done+"1421",trim_done+"1430",trim_done+"1431",trim_done+"1432",trim_done+"1433",trim_done+"1434",trim_done+"1435",trim_done+"1437",trim_done+"1438",trim_done+"1439",trim_done+"1440",trim_done+"1442",trim_done+"1445",trim_done+"1449",trim_done+"1459",trim_done+"1460",trim_done+"1461",trim_done+"1462",trim_done+"1463",trim_done+"1465",trim_done+"1469",trim_done+"1470",trim_done+"1475",trim_done+"1476",trim_done+"1477",trim_done+"1478",trim_done+"1488",trim_done+"1489",trim_done+"1492",trim_done+"2001",trim_done+"2002",trim_done+"2003",trim_done+"2004",trim_done+"2005",trim_done+"2006",trim_done+"2007",trim_done+"2008",trim_done+"2009",trim_done+"2010",trim_done+"2011",trim_done+"2018",trim_done+"2019",trim_done+"2028",trim_done+"2035",trim_done+"2037",trim_done+"2038",trim_done+"2039",trim_done+"2040",trim_done+"2041",trim_done+"2042",trim_done+"2043",trim_done+"2044",trim_done+"2045",trim_done+"2046",trim_done+"2047",trim_done+"2048",trim_done+"2049",trim_done+"2050",trim_done+"2055",trim_done+"2058",trim_done+"2059",trim_done+"2060",trim_done+"2061",trim_done+"2062",trim_done+"2063",trim_done+"2064",trim_done+"2065",trim_done+"2066",trim_done+"2067",trim_done+"2068",trim_done+"2069",trim_done+"3000",trim_done+"3002",trim_done+"3004",trim_done+"3006",trim_done+"3008",trim_done+"3010",trim_done+"3012",trim_done+"3014",trim_done+"3016",trim_done+"3018",trim_done+"3020",trim_done+"3022",trim_done+"3024",trim_done+"3026",trim_done+"3028",trim_done+"3030",trim_done+"3032",trim_done+"3034",trim_done+"3036",trim_done+"3038",trim_done+"3040",trim_done+"3042",trim_done+"3044",trim_done+"3046",trim_done+"3048",trim_done+"3050",trim_done+"3052",trim_done+"3054",trim_done+"3056",trim_done+"3058",trim_done+"3060",trim_done+"3062",trim_done+"3064",trim_done+"3066",trim_done+"3068",trim_done+"3070",trim_done+"3072",trim_done+"3074",trim_done+"3076",trim_done+"3078",trim_done+"3080",trim_done+"3082",trim_done+"3084",trim_done+"3086",trim_done+"3088",trim_done+"3090",trim_done+"3092",trim_done+"3094",trim_done+"3096",trim_done+"3098",trim_done+"3100",trim_done+"3102",trim_done+"3104",trim_done+"3106",trim_done+"3108",trim_done+"3110",trim_done+"3112",trim_done+"3114",trim_done+"3116",trim_done+"3118",trim_done+"3120",trim_done+"3122",trim_done+"3124",trim_done+"3126",trim_done+"3128",trim_done+"3130",trim_done+"3132",trim_done+"3134",trim_done+"3136",trim_done+"3138",trim_done+"3140",trim_done+"3142",trim_done+"3144",trim_done+"3146",trim_done+"3148",trim_done+"3150",trim_done+"3152",trim_done+"3154",trim_done+"3156",trim_done+"3158",trim_done+"3160",trim_done+"3162",trim_done+"3164",trim_done+"3166",trim_done+"3168",trim_done+"3170",trim_done+"3172",trim_done+"3174",trim_done+"3176",trim_done+"3178",trim_done+"3180",trim_done+"3182",trim_done+"3184",trim_done+"3186",trim_done+"3188",trim_done+"3190",trim_done+"3192",trim_done+"3194",trim_done+"3196",trim_done+"3198",trim_done+"3200",trim_done+"3202",trim_done+"3204",trim_done+"3206",trim_done+"3208",trim_done+"3210",trim_done+"3212",trim_done+"3214",trim_done+"3216",trim_done+"3218",trim_done+"3220",trim_done+"3222",trim_done+"3224",trim_done+"3226",trim_done+"3228",trim_done+"3230",trim_done+"3232",trim_done+"3234",trim_done+"3236",trim_done+"3238",trim_done+"3240",trim_done+"3242",trim_done+"3244",trim_done+"3246",trim_done+"3248",trim_done+"3250",trim_done+"3252",trim_done+"3254",trim_done+"3256",trim_done+"3258",trim_done+"3260",trim_done+"3262",trim_done+"3264",trim_done+"3266",trim_done+"3268",trim_done+"3270",trim_done+"3272",trim_done+"3274",trim_done+"3276",trim_done+"3278",trim_done+"3280",trim_done+"3282",trim_done+"3284",trim_done+"3286",trim_done+"3288",trim_done+"3290",trim_done+"3292",trim_done+"3294",trim_done+"3296",trim_done+"3298",trim_done+"3300",trim_done+"3302",trim_done+"3304",trim_done+"3306",trim_done+"3308",trim_done+"3310",trim_done+"3312",trim_done+"3314",trim_done+"3316",trim_done+"3318",trim_done+"3320",trim_done+"3322",trim_done+"3324",trim_done+"3326",trim_done+"3328",trim_done+"3330",trim_done+"3332",trim_done+"3334",trim_done+"3336",trim_done+"3338",trim_done+"3340",trim_done+"3342",trim_done+"3344",trim_done+"3346",trim_done+"3348",trim_done+"3350",trim_done+"3352",trim_done+"3354",trim_done+"3356",trim_done+"3358",trim_done+"3360",trim_done+"3362",trim_done+"3364",trim_done+"3366",trim_done+"3368",trim_done+"3370",trim_done+"3372",trim_done+"3374",trim_done+"3376",trim_done+"3378",trim_done+"3380",trim_done+"3382",trim_done+"3384",trim_done+"3386",trim_done+"3388",trim_done+"3390",trim_done+"3392",trim_done+"3394",trim_done+"3396",trim_done+"3398",trim_done+"3400",trim_done+"3402",trim_done+"3404",trim_done+"3406",trim_done+"3408",trim_done+"3410",trim_done+"3412",trim_done+"3414",trim_done+"3416",trim_done+"3418",trim_done+"3420",trim_done+"3422",trim_done+"3424",trim_done+"3426",trim_done+"3428",trim_done+"3430",trim_done+"3432",trim_done+"3434",trim_done+"3436",trim_done+"3438",trim_done+"3440",trim_done+"3442",trim_done+"3444",trim_done+"4000",trim_done+"4001",trim_done+"4002",trim_done+"4003",trim_done+"4004",trim_done+"4005",trim_done+"4006",trim_done+"4007",trim_done+"4008",trim_done+"4009",trim_done+"4010",trim_done+"4011",trim_done+"4012",trim_done+"4013",trim_done+"4014",trim_done+"4015",trim_done+"4016",trim_done+"4017",trim_done+"4018",trim_done+"4019",trim_done+"4020",trim_done+"4021",trim_done+"4022",trim_done+"4023",trim_done+"4024",trim_done+"4025",trim_done+"4026",trim_done+"4027",trim_done+"4028",trim_done+"4029",trim_done+"4030",trim_done+"4031",trim_done+"4032",trim_done+"4033",trim_done+"4034",trim_done+"4035",trim_done+"4036",trim_done+"4037",trim_done+"4038",trim_done+"4039",trim_done+"4040",trim_done+"4052",trim_done+"4053",trim_done+"4054",trim_done+"4055",trim_done+"4056",trim_done+"4057",trim_done+"4058",trim_done+"4059",trim_done+"4060",trim_done+"4061",trim_done+"4062",trim_done+"4063",trim_done+"4064",trim_done+"4065",trim_done+"4066",trim_done+"4067",trim_done+"4068",trim_done+"4069",trim_done+"4070",trim_done+"4071",trim_done+"4072",trim_done+"4073",trim_done+"4074",trim_done+"4075",trim_done+"4076",trim_done+"4077",trim_done+"4078",trim_done+"4079",trim_done+"4080",trim_done+"4081",trim_done+"4082",trim_done+"4083",trim_done+"4084",trim_done+"4085",trim_done+"4086",trim_done+"4087",trim_done+"4088",trim_done+"4089"]
#     outputs = [path_out, done]
#     options = {'cores': 5, 'memory': "20g", 'walltime': "120:00:00", 'account':"Coryphoideae"}


#     spec = """
#     source activate secapr_env

#     cd /home/owrisberg/Coryphoideae/work_flow/02_trimmed/secapr_postrim

#     secapr quality_check --input {input} --output {output}
    
#     touch {done}

#     """.format(input = path_in, output = path_out, done = done)

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
    echo activating base environment

    #This line should enable the activation of specific conda environments
    source /home/owrisberg/miniconda3/etc/profile.d/conda.sh

    conda activate base

    echo cding to slurm job id
    cd /scratch/$SLURM_JOBID
        
    /home/owrisberg/Coryphoideae/github_code/HybPiper/reads_first.py --cpu 1 --readfiles {p1} {p2} --unpaired {un} -b /home/owrisberg/Coryphoideae/target_sequence/PhyloPalms_loci_renamed_794-176_HEYcorrected.fasta --prefix {species} --bwa

    mv {species} /home/owrisberg/Coryphoideae/work_flow/03_hybpiper/

    touch {done}

    """.format(species=species, p1 = path_in + species + p1,p2 = path_in + species + p2, un = path_in + species + un , out = path_out, done = done)


    return (inputs, outputs, options, spec)

########################################################################################################################
#############################################---- Paralogs ----#########################################################
########################################################################################################################

def paralogs(species,path_in, done, no_paralogs, in_done):
    """Find Paralog genes and write them in the file called paralog.txt"""
    inputs = [path_in + species, in_done]
    outputs = [done]
    options = {'cores': 2, 'memory': "10g", 'walltime': "0:30:00", 'account':"Coryphoideae"}

    spec = """
    source activate base
    
    if test -f /home/owrisberg/Coryphoideae/work_flow/03_hybpiper/{sp}/genes_with_paralog_warnings.txt; then
        echo "/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/{sp}/genes_with_paralog_warnings.txt exists" 
        cd {path_in}
        python /home/owrisberg/Coryphoideae/github_code/HybPiper/paralog_investigator.py {sp} 2>> paralog.txt
    else
        echo "the genes_with_paralog_warnings.txt does not exist and we run the no parallels part"
        touch {np}
    fi
    
    touch {done}

    """.format(sp = species, done = done, path_in = path_in, np = no_paralogs)
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
    inputs = [path_in + species, path_in+"done/Hybpiper/"+species]
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
def coverage(species, path_in, path_out, done,all_bam,all_sorted_bam, all_sorted_bam_bai, bam, cov,fasta,fasta_amb,fasta_ann,fasta_bwt,fasta_pac,fasta_sa,trimmed_fasta,up_bam,dir_in,dir_out):
    """Calculating coverage of sequences."""
    inputs = [path_in+species, path_in+"done/Intronerate/"+species]
    outputs = [path_out+species+all_bam, path_out+species+all_sorted_bam, path_out+species+all_sorted_bam_bai, path_out+species+bam,
    path_out+species+cov, path_out+species+fasta, path_out+species+fasta_amb, path_out+species+fasta_ann, path_out+species+fasta_bwt,
    path_out+species+fasta_pac, path_out+species+fasta_sa, path_out+species+trimmed_fasta, path_out+species+up_bam,done] #ALL the output files
    options = {'cores': 4, 'memory': "20g", 'walltime': "08:00:00", 'account':"Coryphoideae"}

    spec = """
    source activate base
    
    cd {path_in}

    python3 /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/coverage.py {sp} {dir_in} {dir_out}
    
    touch {done}

    """.format(sp = species, done = done, path_in = path_in, dir_in = dir_in, dir_out = dir_out)

    return (inputs, outputs, options, spec)

########################################################################################################################
######################################################---- RUN ----#####################################################
########################################################################################################################
sp = ["1001","1002","1003","1004","1005","1006","1007","1008","1009","1010","1011","1012","1013","1014","1015","1016","1017","1018","1019","1020","1021","1022","1023","1024","1025","1026","1027","1028","1029","1030","1031","1032","1033","1034","1035","1036","1037","1038","1039","1040","1041","1042","1043","1044","1045","1046","1047","1048","1049","1050","1051","1052","1053","1054","1055","1056","1057","1058","1059","1060","1061","1062","1063","1064","1065","1066","1067","1068","1069","1070","1071","1072","1073","1074","1075","1076","1077","1078","1079","1080","1099","1081","1082","1193","1194","1195","1196","1197","1222","1224","1237","1238","1239","1244","1245","1246","1247","1248","1249","1250","1251","1252","1255","1256","1257","1258","1259","1266","1267","1268","1269","1273","1274","1275","1276","1344","1345","1346","1347","1348","1349","1350","1351","1352","1353","1354","1355","1362","1363","1364","1411","1412","1413","1414","1415","1416","1417","1418","1419","1421","1430","1431","1432","1433","1434","1435","1437","1438","1439","1440","1442","1443","1444","1445","1447","1449","1459","1460","1461","1462","1463","1465","1469","1470","1475","1476","1477","1478","1488","1489","1492","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2018","2019","2028","2035","2037","2038","2039","2040","2041","2042","2043","2044","2045","2046","2047","2048","2049","2050","2055","2058","2059","2060","2061","2062","2063","2064","2065","2066","2067","2068","2069","3000","3002","3004","3006","3008","3010","3012","3014","3016","3018","3020","3022","3024","3026","3028","3030","3032","3034","3036","3038","3040","3042","3044","3046","3048","3050","3052","3054","3056","3058","3060","3062","3064","3066","3068","3070","3072","3074","3076","3078","3080","3082","3084","3086","3088","3090","3092","3094","3096","3098","3100","3102","3104","3106","3108","3110","3112","3114","3116","3118","3120","3122","3124","3126","3128","3130","3132","3134","3136","3138","3140","3142","3144","3146","3148","3150","3152","3154","3156","3158","3160","3162","3164","3166","3168","3170","3172","3174","3176","3178","3180","3182","3184","3186","3188","3190","3192","3194","3196","3198","3200","3202","3204","3206","3208","3210","3212","3214","3216","3218","3220","3222","3224","3226","3228","3230","3232","3234","3236","3238","3240","3242","3244","3246","3248","3250","3252","3254","3256","3258","3260","3262","3264","3266","3268","3270","3272","3274","3276","3278","3280","3282","3284","3286","3288","3290","3292","3294","3296","3298","3300","3302","3304","3306","3308","3310","3312","3314","3316","3318","3320","3322","3324","3326","3328","3330","3332","3334","3336","3338","3340","3342","3344","3346","3348","3350","3352","3354","3356","3358","3360","3362","3364","3366","3368","3370","3372","3374","3376","3378","3380","3382","3384","3386","3388","3390","3392","3394","3396","3398","3400","3402","3404","3406","3408","3410","3412","3414","3416","3418","3420","3422","3424","3426","3428","3430","3432","3434","3436","3438","3440","3442","3444","4000","4001","4002","4003","4004","4005","4006","4007","4008","4009","4010","4011","4012","4013","4014","4015","4016","4017","4018","4019","4020","4021","4022","4023","4024","4025","4026","4027","4028","4029","4030","4031","4032","4033","4034","4035","4036","4037","4038","4039","4040","4052","4053","4054","4055","4056","4057","4058","4059","4060","4061","4062","4063","4064","4065","4066","4067","4068","4069","4070","4071","4072","4073","4074","4075","4076","4077","4078","4079","4080","4081","4082","4083","4084","4085","4086","4087","4088","4089"]

# #### Running Secapr quality check on raw data
# gwf.target_from_template('Quality_check_Raw', secapr_quality_check_raw(path_out= "/home/owrisberg/Coryphoideae/work_flow/00_secapr/0_data/",
#                                                         path_in = "/home/owrisberg/Coryphoideae/work_flow/01_data/",
#                                                         done = "/home/owrisberg/Coryphoideae/work_flow/00_secapr/done/secapr_raw"))



# #### Running Secapr quality check on trimmed data
# gwf.target_from_template('Quality_check_trimmed', secapr_quality_check_trimmed(path_out= "/home/owrisberg/Coryphoideae/work_flow/00_secapr/1_trimmed/",
#                                                         path_in = "/home/owrisberg/Coryphoideae/work_flow/02_trimmed/",
#                                                         done = "/home/owrisberg/Coryphoideae/work_flow/00_secapr/done/secapr_trimmed",
#                                                         trim_done = "/home/owrisberg/Coryphoideae/work_flow/02_trimmed/done/"))



for i in range(len(sp)):
    #### Running fastqc on raw data
    gwf.target_from_template('fastqc_raw_'+sp[i], fastqc_raw(species = sp[i],
                                                        path_in= "/home/owrisberg/Coryphoideae/work_flow/01_data/",
                                                        path_out = "/home/owrisberg/Coryphoideae/work_flow/00_secapr/0_data/",
                                                        done = "/home/owrisberg/Coryphoideae/work_flow/00_secapr/done/raw_data/"+sp[i]))


    #### Running Trimmomatic
    gwf.target_from_template('trimmomatic_'+sp[i], trimmomatic(species = sp[i],
                                                        path_in= "/home/owrisberg/Coryphoideae/work_flow/01_data/",
                                                        path_out = "/home/owrisberg/Coryphoideae/work_flow/02_trimmed/",
                                                        done = "/home/owrisberg/Coryphoideae/work_flow/00_secapr/done/trimmed_data"+sp[i]))

    #### Running fastqc on the trimmed data
    gwf.target_from_template('fastqc_trimmed_'+sp[i], fastqc_trimmed(species = sp[i],
                                                        path_in= "/home/owrisberg/Coryphoideae/work_flow/02_trimmed/secapr_postrim/",
                                                        path_out = "/home/owrisberg/Coryphoideae/work_flow/00_secapr/1_trimmed/",
                                                        done = "/home/owrisberg/Coryphoideae/work_flow/00_secapr/done/trimmed_data/"+sp[i]))                                                   

    #### Running Hybpiper
    gwf.target_from_template('Hybpiper_'+sp[i], hybpiper(species = sp[i],
                                                        p1 = "_1P.fastq",
                                                        p2 = "_2P.fastq",
                                                        un = "_UN.fastq",
                                                        path_out= "/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/",
                                                        path_in = "/home/owrisberg/Coryphoideae/work_flow/02_trimmed/",
                                                        done = "/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/done/Hybpiper/"+sp[i],))
                                                                      

    #### Paralogs
    
    gwf.target_from_template('Paralogs_'+sp[i], paralogs(species = sp[i],
                                                        path_in = "/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/",
                                                        done = "/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/done/Paralogs/"+sp[i],
                                                        no_paralogs="/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/done/No_paralogs/"+sp[i],
                                                        in_done="/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/done/Hybpiper/"+sp[i]))
    # else:
    #     gwf.target_from_template('No_Paralogs_'+sp[i], no_paralogs(species = sp[i],
    #                                                             path_in = "/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/",
    #                                                             done = "/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/done/Paralogs/"+sp[i],
    #                                                             no_paralogs="/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/done/No_paralogs/"+sp[i]))
     
    
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
                                                        done = "/home/owrisberg/Coryphoideae/work_flow/04_coverage/done/Coverage/"+sp[i],
                                                        dir_in ="/home/owrisberg/Coryphoideae/work_flow/02_trimmed/", #Folder with raw reads
                                                        dir_out ="/home/owrisberg/Coryphoideae/work_flow/04_coverage/")) # folder with coverage

