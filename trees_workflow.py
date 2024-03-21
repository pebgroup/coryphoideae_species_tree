'''
------------------------------------------------------------------------------------------------------------------------
***GUIDE***:
This workflow is used to build the trees from the manually aligned gene sequences.
It should automatically check if the genes in the manual aligned gene folder are newer than the ones used to create the
last phylogeny, and if they are it should run the whole pipe again.

------------------------------------------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------------------------------------
Author: Oscar Wrisberg
Date: 7/12/2021
------------------------------------------------------------------------------------------------------------------------
'''

from os import O_SYNC
from gwf import Workflow, AnonymousTarget
import os.path
# import math
# import glob

gwf = Workflow()

##########################################################################################################################
# ########################---- Copying alignments and creating partitions ----############################################
# ########################################################################################################################

def partitioner(path_in,path_out, gene, done):
    """Copying alignments from the manual alignment folder to the treebuilding folder and creating partition files"""
    inputs = ["/home/owrisberg/Coryphoideae/work_flow/09_mapping/"+gene+"_output_tapper_mapped.fasta"]
    outputs = [path_out+gene+"_part.txt",path_out+gene+"_clean.fasta", done]
    options = {'cores': 1, 'memory': "5g", 'walltime': "00:30:00", 'account':"Coryphoideae"}

    spec = """
    source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
    conda activate treebuilder_env
	#Going to folder with data
	cd {path_in}
    
	#The partitioner should produce 2 files for each gene
	#one file called gene_aligned_part.txt which is the partitioning file
	#another called gene_aligned_clean.fasta which are just the sequences without the exons

	python3 /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/partitioner.py --smoother 10 --gene {gene} --file_ending _output_tapper_mapped.fasta
	
	#Moving files to the correct folders
	mv {gene}_clean.fasta {path_out}
	mv {gene}_part.txt {path_out} 

	touch {done}


    """.format(path_in = path_in, gene = gene, done = done, path_out = path_out)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# ########################################################################################################################
# #############################################---- IQ-tree ----#############################################################
# ########################################################################################################################

def iq_tree(path_in, gene,path_out ):
    """Using Iq-tree to produce trees for each gene with a partition file to use individual substitution rates for each gene"""
    inputs = [path_in+gene+"_part.txt", path_in+gene+"_clean.fasta"]
    outputs = [path_out+gene+".txt.tre"]
    options = {'cores': 20, 'memory': "20g", 'walltime': "24:00:00", 'account':"Coryphoideae"}

    spec = """
	source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
	conda activate treebuilder_env

	cd {path_in}
    
    echo "Running IQ-tree for {gene} at:"
    date

	#Actual IQtree tree search. 
	iqtree2 -s {gene}_clean.fasta -p {gene}_part.txt -T AUTO -ntmax 20 -m MFP -B 1000 -redo 



	mv {gene}*.treefile {path_out}{gene}.txt.tre
	


	""".format(path_in = path_in, gene = gene, path_out=path_out)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# ########################################################################################################################
# #####################################---- Renaming & Rerooting ----#####################################################
# ########################################################################################################################
def rename_reroot(path_in, gene):
    """Using genenameremover to remove the gene names from all the tip labels
	And using rerooter.py to root each individual gene tree based on the available outgroup"""
    inputs = [path_in+gene+".txt.tre"] # changed _part.txt.tre to .txt.tre
    outputs = [path_in+gene+"_rooted.tre"]
    options = {'cores': 5, 'memory': "10g", 'walltime': "00:30:00", 'account':"Coryphoideae"}

    spec = """
	source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
	conda activate treebuilder_env

	cd {path_in}

	echo Removing {gene} from tip labels
	python3 /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/genenameremover.py --gene {gene} --treefile {gene}.txt.tre

	#Removing _R_ from sequences which have been reversed
	sed -i -e 's/_R_//g' {gene}.txt.tre

	echo Rerooting each genetree based on the outgroup	
	python3 /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/rooter.py --gene {gene} --treefile {gene}.txt.tre 


	""".format(path_in = path_in, gene = gene)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# ########################################################################################################################
# ############################---- Newick Contracting and gathering of genetrees ----#####################################
# ########################################################################################################################

def newick_contracting(path_in,path_out ):
    """Gathering all the gene-trees in a single file and removing branches with low bootstrap support
	I have had to remove the following genes for varius reasons such as outgroup species missing and low gene recovery and bad optrimal values"""
    inputs = [ path_in+"EGU105032175_rooted.tre",path_in+"EGU105032229_rooted.tre",path_in+"EGU105032337_rooted.tre",path_in+"EGU105032379_rooted.tre",path_in+"EGU105033063_rooted.tre",path_in+"EGU105033626_rooted.tre",path_in+"EGU105034121_rooted.tre",path_in+"EGU105034616_rooted.tre",path_in+"EGU105034893_rooted.tre",path_in+"EGU105034993_rooted.tre",path_in+"EGU105035046_rooted.tre",path_in+"EGU105035196_rooted.tre",path_in+"EGU105035203_rooted.tre",path_in+"EGU105035462_rooted.tre",path_in+"EGU105035555_rooted.tre",path_in+"EGU105035989_rooted.tre",path_in+"EGU105036031_rooted.tre",path_in+"EGU105036385_rooted.tre",path_in+"EGU105036774_rooted.tre",path_in+"EGU105037749_rooted.tre",path_in+"EGU105037800_rooted.tre",path_in+"EGU105037890_rooted.tre",path_in+"EGU105037902_rooted.tre",path_in+"EGU105037930_rooted.tre",path_in+"EGU105037938_rooted.tre",path_in+"EGU105038008_rooted.tre",path_in+"EGU105038036_rooted.tre",path_in+"EGU105038098_rooted.tre",path_in+"EGU105038099_rooted.tre",path_in+"EGU105038100_rooted.tre",path_in+"EGU105038110_rooted.tre",path_in+"EGU105038114_rooted.tre",path_in+"EGU105038118_rooted.tre",path_in+"EGU105038123_rooted.tre",path_in+"EGU105038179_rooted.tre",path_in+"EGU105038201_rooted.tre",path_in+"EGU105038228_rooted.tre",path_in+"EGU105038234_rooted.tre",path_in+"EGU105038245_rooted.tre",path_in+"EGU105038252_rooted.tre",path_in+"EGU105038310_rooted.tre",path_in+"EGU105038382_rooted.tre",path_in+"EGU105038400_rooted.tre",path_in+"EGU105038419_rooted.tre",path_in+"EGU105038431_rooted.tre",path_in+"EGU105038499_rooted.tre",path_in+"EGU105038513_rooted.tre",path_in+"EGU105038571_rooted.tre",path_in+"EGU105038580_rooted.tre",path_in+"EGU105038603_rooted.tre",path_in+"EGU105038631_rooted.tre",path_in+"EGU105038680_rooted.tre",path_in+"EGU105038693_rooted.tre",path_in+"EGU105038720_rooted.tre",path_in+"EGU105038747_rooted.tre",path_in+"EGU105038794_rooted.tre",path_in+"EGU105038832_rooted.tre",path_in+"EGU105038882_rooted.tre",path_in+"EGU105038986_rooted.tre",path_in+"EGU105038988_rooted.tre",path_in+"EGU105039013_rooted.tre",path_in+"EGU105039062_rooted.tre",path_in+"EGU105039067_rooted.tre",path_in+"EGU105039082_rooted.tre",path_in+"EGU105039099_rooted.tre",path_in+"EGU105039101_rooted.tre",path_in+"EGU105039107_rooted.tre",path_in+"EGU105039121_rooted.tre",path_in+"EGU105039164_rooted.tre",path_in+"EGU105039178_rooted.tre",path_in+"EGU105039221_rooted.tre",path_in+"EGU105039236_rooted.tre",path_in+"EGU105039255_rooted.tre",path_in+"EGU105039282_rooted.tre",path_in+"EGU105039298_rooted.tre",path_in+"EGU105039313_rooted.tre",path_in+"EGU105039403_rooted.tre",path_in+"EGU105039431_rooted.tre",path_in+"EGU105039449_rooted.tre",path_in+"EGU105039460_rooted.tre",path_in+"EGU105039480_rooted.tre",path_in+"EGU105039494_rooted.tre",path_in+"EGU105039501_rooted.tre",path_in+"EGU105039512_rooted.tre",path_in+"EGU105039542_rooted.tre",path_in+"EGU105039587_rooted.tre",path_in+"EGU105039595_rooted.tre",path_in+"EGU105039609_rooted.tre",path_in+"EGU105039660_rooted.tre",path_in+"EGU105039685_rooted.tre",path_in+"EGU105039690_rooted.tre",path_in+"EGU105039699_rooted.tre",path_in+"EGU105039763_rooted.tre",path_in+"EGU105039783_rooted.tre",path_in+"EGU105039809_rooted.tre",path_in+"EGU105039822_rooted.tre",path_in+"EGU105039925_rooted.tre",path_in+"EGU105039947_rooted.tre",path_in+"EGU105039957_rooted.tre",path_in+"EGU105039962_rooted.tre",path_in+"EGU105040073_rooted.tre",path_in+"EGU105040088_rooted.tre",path_in+"EGU105040099_rooted.tre",path_in+"EGU105040114_rooted.tre",path_in+"EGU105040115_rooted.tre",path_in+"EGU105040125_rooted.tre",path_in+"EGU105040139_rooted.tre",path_in+"EGU105040185_rooted.tre",path_in+"EGU105040186_rooted.tre",path_in+"EGU105040189_rooted.tre",path_in+"EGU105040206_rooted.tre",path_in+"EGU105040207_rooted.tre",path_in+"EGU105040242_rooted.tre",path_in+"EGU105040281_rooted.tre",path_in+"EGU105040302_rooted.tre",path_in+"EGU105040308_rooted.tre",path_in+"EGU105040359_rooted.tre",path_in+"EGU105040368_rooted.tre",path_in+"EGU105040426_rooted.tre",path_in+"EGU105040452_rooted.tre",path_in+"EGU105040462_rooted.tre",path_in+"EGU105040530_rooted.tre",path_in+"EGU105040583_rooted.tre",path_in+"EGU105040667_rooted.tre",path_in+"EGU105040675_rooted.tre",path_in+"EGU105040684_rooted.tre",path_in+"EGU105040690_rooted.tre",path_in+"EGU105040700_rooted.tre",path_in+"EGU105040756_rooted.tre",path_in+"EGU105040758_rooted.tre",path_in+"EGU105040813_rooted.tre",path_in+"EGU105040837_rooted.tre",path_in+"EGU105040842_rooted.tre",path_in+"EGU105040850_rooted.tre",path_in+"EGU105040851_rooted.tre",path_in+"EGU105040863_rooted.tre",path_in+"EGU105040887_rooted.tre",path_in+"EGU105040914_rooted.tre",path_in+"EGU105040918_rooted.tre",path_in+"EGU105040922_rooted.tre",path_in+"EGU105040957_rooted.tre",path_in+"EGU105040970_rooted.tre",path_in+"EGU105041055_rooted.tre",path_in+"EGU105041100_rooted.tre",path_in+"EGU105041117_rooted.tre",path_in+"EGU105041125_rooted.tre",path_in+"EGU105041127_rooted.tre",path_in+"EGU105041133_rooted.tre",path_in+"EGU105041179_rooted.tre",path_in+"EGU105041182_rooted.tre",path_in+"EGU105041189_rooted.tre",path_in+"EGU105041217_rooted.tre",path_in+"EGU105041246_rooted.tre",path_in+"EGU105041283_rooted.tre",path_in+"EGU105041337_rooted.tre",path_in+"EGU105041353_rooted.tre",path_in+"EGU105041650_rooted.tre",path_in+"EGU105041657_rooted.tre",path_in+"EGU105041665_rooted.tre",path_in+"EGU105041680_rooted.tre",path_in+"EGU105041687_rooted.tre",path_in+"EGU105041710_rooted.tre",path_in+"EGU105041807_rooted.tre",path_in+"EGU105041816_rooted.tre",path_in+"EGU105041872_rooted.tre",path_in+"EGU105041902_rooted.tre",path_in+"EGU105041903_rooted.tre",path_in+"EGU105041929_rooted.tre",path_in+"EGU105041933_rooted.tre",path_in+"EGU105041982_rooted.tre",path_in+"EGU105042090_rooted.tre",path_in+"EGU105042113_rooted.tre",path_in+"EGU105042128_rooted.tre",path_in+"EGU105042147_rooted.tre",path_in+"EGU105042168_rooted.tre",path_in+"EGU105042205_rooted.tre",path_in+"EGU105042290_rooted.tre",path_in+"EGU105042307_rooted.tre",path_in+"EGU105042323_rooted.tre",path_in+"EGU105042329_rooted.tre",path_in+"EGU105042368_rooted.tre",path_in+"EGU105042422_rooted.tre",path_in+"EGU105042525_rooted.tre",path_in+"EGU105042558_rooted.tre",path_in+"EGU105042560_rooted.tre",path_in+"EGU105042584_rooted.tre",path_in+"EGU105042633_rooted.tre",path_in+"EGU105042644_rooted.tre",path_in+"EGU105042651_rooted.tre",path_in+"EGU105042664_rooted.tre",path_in+"EGU105042722_rooted.tre",path_in+"EGU105042781_rooted.tre",path_in+"EGU105042808_rooted.tre",path_in+"EGU105042820_rooted.tre",path_in+"EGU105042873_rooted.tre",path_in+"EGU105042965_rooted.tre",path_in+"EGU105043011_rooted.tre",path_in+"EGU105043037_rooted.tre",path_in+"EGU105043042_rooted.tre",path_in+"EGU105043061_rooted.tre",path_in+"EGU105043069_rooted.tre",path_in+"EGU105043119_rooted.tre",path_in+"EGU105043155_rooted.tre",path_in+"EGU105043160_rooted.tre",path_in+"EGU105043164_rooted.tre",path_in+"EGU105043193_rooted.tre",path_in+"EGU105043320_rooted.tre",path_in+"EGU105043338_rooted.tre",path_in+"EGU105043374_rooted.tre",path_in+"EGU105043419_rooted.tre",path_in+"EGU105043430_rooted.tre",path_in+"EGU105043469_rooted.tre",path_in+"EGU105043485_rooted.tre",path_in+"EGU105043499_rooted.tre",path_in+"EGU105043601_rooted.tre",path_in+"EGU105043633_rooted.tre",path_in+"EGU105043666_rooted.tre",path_in+"EGU105043685_rooted.tre",path_in+"EGU105043686_rooted.tre",path_in+"EGU105043730_rooted.tre",path_in+"EGU105043786_rooted.tre",path_in+"EGU105043816_rooted.tre",path_in+"EGU105043827_rooted.tre",path_in+"EGU105043926_rooted.tre",path_in+"EGU105043975_rooted.tre",path_in+"EGU105044063_rooted.tre",path_in+"EGU105044120_rooted.tre",path_in+"EGU105044133_rooted.tre",path_in+"EGU105044174_rooted.tre",path_in+"EGU105044182_rooted.tre",path_in+"EGU105044183_rooted.tre",path_in+"EGU105044203_rooted.tre",path_in+"EGU105044252_rooted.tre",path_in+"EGU105044281_rooted.tre",path_in+"EGU105044307_rooted.tre",path_in+"EGU105044309_rooted.tre",path_in+"EGU105044350_rooted.tre",path_in+"EGU105044378_rooted.tre",path_in+"EGU105044400_rooted.tre",path_in+"EGU105044407_rooted.tre",path_in+"EGU105044445_rooted.tre",path_in+"EGU105044446_rooted.tre",path_in+"EGU105044481_rooted.tre",path_in+"EGU105044588_rooted.tre",path_in+"EGU105044613_rooted.tre",path_in+"EGU105044614_rooted.tre",path_in+"EGU105044668_rooted.tre",path_in+"EGU105044676_rooted.tre",path_in+"EGU105044710_rooted.tre",path_in+"EGU105044758_rooted.tre",path_in+"EGU105044844_rooted.tre",path_in+"EGU105044846_rooted.tre",path_in+"EGU105044854_rooted.tre",path_in+"EGU105044885_rooted.tre",path_in+"EGU105044893_rooted.tre",path_in+"EGU105044896_rooted.tre",path_in+"EGU105044978_rooted.tre",path_in+"EGU105044982_rooted.tre",path_in+"EGU105044983_rooted.tre",path_in+"EGU105044984_rooted.tre",path_in+"EGU105045005_rooted.tre",path_in+"EGU105045043_rooted.tre",path_in+"EGU105045070_rooted.tre",path_in+"EGU105045078_rooted.tre",path_in+"EGU105045094_rooted.tre",path_in+"EGU105045099_rooted.tre",path_in+"EGU105045102_rooted.tre",path_in+"EGU105045137_rooted.tre",path_in+"EGU105045148_rooted.tre",path_in+"EGU105045232_rooted.tre",path_in+"EGU105045248_rooted.tre",path_in+"EGU105045254_rooted.tre",path_in+"EGU105045282_rooted.tre",path_in+"EGU105045310_rooted.tre",path_in+"EGU105045358_rooted.tre",path_in+"EGU105045367_rooted.tre",path_in+"EGU105045424_rooted.tre",path_in+"EGU105045464_rooted.tre",path_in+"EGU105045467_rooted.tre",path_in+"EGU105045489_rooted.tre",path_in+"EGU105045507_rooted.tre",path_in+"EGU105045509_rooted.tre",path_in+"EGU105045514_rooted.tre",path_in+"EGU105045520_rooted.tre",path_in+"EGU105045529_rooted.tre",path_in+"EGU105045544_rooted.tre",path_in+"EGU105045640_rooted.tre",path_in+"EGU105045658_rooted.tre",path_in+"EGU105045703_rooted.tre",path_in+"EGU105045726_rooted.tre",path_in+"EGU105045732_rooted.tre",path_in+"EGU105045760_rooted.tre",path_in+"EGU105045782_rooted.tre",path_in+"EGU105045788_rooted.tre",path_in+"EGU105045820_rooted.tre",path_in+"EGU105045827_rooted.tre",path_in+"EGU105045828_rooted.tre",path_in+"EGU105045835_rooted.tre",path_in+"EGU105045898_rooted.tre",path_in+"EGU105045932_rooted.tre",path_in+"EGU105045946_rooted.tre",path_in+"EGU105046030_rooted.tre",path_in+"EGU105046050_rooted.tre",path_in+"EGU105046056_rooted.tre",path_in+"EGU105046099_rooted.tre",path_in+"EGU105046103_rooted.tre",path_in+"EGU105046147_rooted.tre",path_in+"EGU105046168_rooted.tre",path_in+"EGU105046245_rooted.tre",path_in+"EGU105046297_rooted.tre",path_in+"EGU105046360_rooted.tre",path_in+
	"EGU105046387_rooted.tre",path_in+"EGU105046393_rooted.tre",path_in+"EGU105046401_rooted.tre",path_in+"EGU105046449_rooted.tre",path_in+"EGU105046454_rooted.tre",path_in+"EGU105046456_rooted.tre",path_in+"EGU105046503_rooted.tre",path_in+"EGU105046518_rooted.tre",path_in+"EGU105046530_rooted.tre",path_in+"EGU105046549_rooted.tre",path_in+"EGU105046559_rooted.tre",path_in+"EGU105046562_rooted.tre",path_in+"EGU105046574_rooted.tre",path_in+"EGU105046630_rooted.tre",path_in+"EGU105046632_rooted.tre",path_in+"EGU105046696_rooted.tre",path_in+"EGU105046735_rooted.tre",path_in+"EGU105046766_rooted.tre",path_in+"EGU105046786_rooted.tre",path_in+"EGU105046827_rooted.tre",path_in+"EGU105046875_rooted.tre",path_in+"EGU105046918_rooted.tre",path_in+"EGU105047024_rooted.tre",path_in+"EGU105047029_rooted.tre",path_in+"EGU105047253_rooted.tre",path_in+"EGU105047288_rooted.tre",path_in+"EGU105047293_rooted.tre",path_in+"EGU105047342_rooted.tre",path_in+"EGU105047357_rooted.tre",path_in+"EGU105047362_rooted.tre",path_in+"EGU105047379_rooted.tre",path_in+"EGU105047385_rooted.tre",path_in+"EGU105047395_rooted.tre",path_in+"EGU105047433_rooted.tre",path_in+"EGU105047434_rooted.tre",path_in+"EGU105047446_rooted.tre",path_in+"EGU105047519_rooted.tre",path_in+"EGU105047533_rooted.tre",path_in+"EGU105047546_rooted.tre",path_in+"EGU105047553_rooted.tre",path_in+"EGU105047578_rooted.tre",path_in+"EGU105047585_rooted.tre",path_in+"EGU105047597_rooted.tre",path_in+"EGU105047621_rooted.tre",path_in+"EGU105047644_rooted.tre",path_in+"EGU105047662_rooted.tre",path_in+"EGU105047689_rooted.tre",path_in+"EGU105047751_rooted.tre",path_in+"EGU105047777_rooted.tre",path_in+"EGU105047790_rooted.tre",path_in+"EGU105047907_rooted.tre",path_in+"EGU105047916_rooted.tre",path_in+"EGU105047922_rooted.tre",path_in+"EGU105047940_rooted.tre",path_in+"EGU105047945_rooted.tre",path_in+"EGU105047970_rooted.tre",path_in+"EGU105048009_rooted.tre",path_in+"EGU105048015_rooted.tre",path_in+"EGU105048028_rooted.tre",path_in+"EGU105048054_rooted.tre",path_in+"EGU105048056_rooted.tre",path_in+"EGU105048129_rooted.tre",path_in+"EGU105048130_rooted.tre",path_in+"EGU105048137_rooted.tre",path_in+"EGU105048159_rooted.tre",path_in+"EGU105048182_rooted.tre",path_in+"EGU105048199_rooted.tre",path_in+"EGU105048300_rooted.tre",path_in+"EGU105048357_rooted.tre",path_in+"EGU105048410_rooted.tre",path_in+"EGU105048474_rooted.tre",path_in+"EGU105048476_rooted.tre",path_in+"EGU105048479_rooted.tre",path_in+"EGU105048484_rooted.tre",path_in+"EGU105048486_rooted.tre",path_in+"EGU105048493_rooted.tre",path_in+"EGU105048527_rooted.tre",path_in+"EGU105048541_rooted.tre",path_in+"EGU105048581_rooted.tre",path_in+"EGU105048612_rooted.tre",path_in+"EGU105048694_rooted.tre",path_in+"EGU105048725_rooted.tre",path_in+"EGU105048751_rooted.tre",path_in+"EGU105048796_rooted.tre",path_in+"EGU105048839_rooted.tre",path_in+"EGU105048867_rooted.tre",path_in+"EGU105048886_rooted.tre",path_in+"EGU105048898_rooted.tre",path_in+"EGU105048909_rooted.tre",path_in+"EGU105048915_rooted.tre",path_in+"EGU105048926_rooted.tre",path_in+"EGU105048961_rooted.tre",path_in+"EGU105048968_rooted.tre",path_in+"EGU105049007_rooted.tre",path_in+"EGU105049016_rooted.tre",path_in+"EGU105049020_rooted.tre",path_in+"EGU105049025_rooted.tre",path_in+"EGU105049052_rooted.tre",path_in+"EGU105049097_rooted.tre",path_in+"EGU105049274_rooted.tre",path_in+"EGU105049312_rooted.tre",path_in+"EGU105049318_rooted.tre",path_in+"EGU105049360_rooted.tre",path_in+"EGU105049426_rooted.tre",path_in+"EGU105049539_rooted.tre",path_in+"EGU105049543_rooted.tre",path_in+"EGU105049583_rooted.tre",path_in+"EGU105049690_rooted.tre",path_in+"EGU105049729_rooted.tre",path_in+"EGU105049737_rooted.tre",path_in+"EGU105049761_rooted.tre",path_in+"EGU105049827_rooted.tre",path_in+"EGU105049882_rooted.tre",path_in+"EGU105049902_rooted.tre",path_in+"EGU105049903_rooted.tre",path_in+"EGU105049934_rooted.tre",path_in+"EGU105049947_rooted.tre",path_in+"EGU105050012_rooted.tre",path_in+"EGU105050023_rooted.tre",path_in+"EGU105050036_rooted.tre",path_in+"EGU105050058_rooted.tre",path_in+"EGU105050114_rooted.tre",path_in+"EGU105050126_rooted.tre",path_in+"EGU105050202_rooted.tre",path_in+"EGU105050207_rooted.tre",path_in+"EGU105050328_rooted.tre",path_in+"EGU105050344_rooted.tre",path_in+"EGU105050362_rooted.tre",path_in+"EGU105050366_rooted.tre",path_in+"EGU105050383_rooted.tre",path_in+"EGU105050387_rooted.tre",path_in+"EGU105050404_rooted.tre",path_in+"EGU105050432_rooted.tre",path_in+"EGU105050450_rooted.tre",path_in+"EGU105050521_rooted.tre",path_in+"EGU105050532_rooted.tre",path_in+"EGU105050644_rooted.tre",path_in+"EGU105050670_rooted.tre",path_in+"EGU105050680_rooted.tre",path_in+"EGU105050681_rooted.tre",path_in+"EGU105050682_rooted.tre",path_in+"EGU105050831_rooted.tre",path_in+"EGU105050841_rooted.tre",path_in+"EGU105050853_rooted.tre",path_in+"EGU105050854_rooted.tre",path_in+"EGU105050961_rooted.tre",path_in+"EGU105050970_rooted.tre",path_in+"EGU105050972_rooted.tre",path_in+"EGU105051087_rooted.tre",path_in+"EGU105051146_rooted.tre",path_in+"EGU105051156_rooted.tre",path_in+"EGU105051188_rooted.tre",path_in+"EGU105051345_rooted.tre",path_in+"EGU105051362_rooted.tre",path_in+"EGU105051366_rooted.tre",path_in+"EGU105051373_rooted.tre",path_in+"EGU105051391_rooted.tre",path_in+"EGU105051395_rooted.tre",path_in+"EGU105051403_rooted.tre",path_in+"EGU105051481_rooted.tre",path_in+"EGU105051499_rooted.tre",path_in+"EGU105051503_rooted.tre",path_in+"EGU105051560_rooted.tre",path_in+"EGU105051564_rooted.tre",path_in+"EGU105051582_rooted.tre",path_in+"EGU105051614_rooted.tre",path_in+"EGU105051677_rooted.tre",path_in+"EGU105051704_rooted.tre",path_in+"EGU105051726_rooted.tre",path_in+"EGU105051740_rooted.tre",path_in+"EGU105051748_rooted.tre",path_in+"EGU105051764_rooted.tre",path_in+"EGU105051795_rooted.tre",path_in+"EGU105051802_rooted.tre",path_in+"EGU105051821_rooted.tre",path_in+"EGU105051823_rooted.tre",path_in+"EGU105051832_rooted.tre",path_in+"EGU105051847_rooted.tre",path_in+"EGU105051857_rooted.tre",path_in+"EGU105051860_rooted.tre",path_in+"EGU105051870_rooted.tre",path_in+"EGU105051891_rooted.tre",path_in+"EGU105051924_rooted.tre",path_in+"EGU105051953_rooted.tre",path_in+"EGU105051985_rooted.tre",path_in+"EGU105052035_rooted.tre",path_in+"EGU105052070_rooted.tre",path_in+"EGU105052170_rooted.tre",path_in+"EGU105052178_rooted.tre",path_in+"EGU105052304_rooted.tre",path_in+"EGU105052307_rooted.tre",path_in+"EGU105052346_rooted.tre",path_in+"EGU105052351_rooted.tre",path_in+"EGU105052386_rooted.tre",path_in+"EGU105052389_rooted.tre",path_in+"EGU105052394_rooted.tre",path_in+"EGU105052428_rooted.tre",path_in+"EGU105052446_rooted.tre",path_in+"EGU105052476_rooted.tre",path_in+"EGU105052483_rooted.tre",path_in+"EGU105052492_rooted.tre",path_in+"EGU105052495_rooted.tre",path_in+"EGU105052527_rooted.tre",path_in+"EGU105052529_rooted.tre",path_in+"EGU105052538_rooted.tre",path_in+"EGU105052552_rooted.tre",path_in+"EGU105052573_rooted.tre",path_in+"EGU105052580_rooted.tre",path_in+"EGU105052623_rooted.tre",path_in+"EGU105052650_rooted.tre",path_in+"EGU105052694_rooted.tre",path_in+"EGU105052704_rooted.tre",path_in+"EGU105052739_rooted.tre",path_in+"EGU105052743_rooted.tre",path_in+"EGU105052750_rooted.tre",path_in+"EGU105052771_rooted.tre",path_in+"EGU105052804_rooted.tre",path_in+"EGU105052818_rooted.tre",path_in+"EGU105052849_rooted.tre",path_in+"EGU105052855_rooted.tre",path_in+"EGU105052865_rooted.tre",path_in+"EGU105052888_rooted.tre",path_in+"EGU105052944_rooted.tre",path_in+"EGU105052947_rooted.tre",path_in+"EGU105052956_rooted.tre",path_in+"EGU105053006_rooted.tre",path_in+"EGU105053055_rooted.tre",path_in+"EGU105053059_rooted.tre",path_in+"EGU105053079_rooted.tre",path_in+"EGU105053105_rooted.tre",path_in+"EGU105053124_rooted.tre",path_in+"EGU105053130_rooted.tre",path_in+"EGU105053136_rooted.tre",path_in+"EGU105053172_rooted.tre",path_in+"EGU105053204_rooted.tre",path_in+"EGU105053227_rooted.tre",path_in+"EGU105053263_rooted.tre",path_in+"EGU105053403_rooted.tre",path_in+"EGU105053422_rooted.tre",path_in+"EGU105053426_rooted.tre",path_in+"EGU105053457_rooted.tre",path_in+"EGU105053465_rooted.tre",path_in+"EGU105053468_rooted.tre",path_in+"EGU105053482_rooted.tre",path_in+"EGU105053549_rooted.tre",path_in+"EGU105053642_rooted.tre",path_in+"EGU105053654_rooted.tre",path_in+"EGU105053735_rooted.tre",path_in+"EGU105053747_rooted.tre",path_in+"EGU105053770_rooted.tre",path_in+"EGU105053835_rooted.tre",path_in+"EGU105053848_rooted.tre",path_in+"EGU105053866_rooted.tre",path_in+"EGU105053889_rooted.tre",path_in+"EGU105053901_rooted.tre",path_in+"EGU105053932_rooted.tre",path_in+"EGU105053961_rooted.tre",path_in+"EGU105053969_rooted.tre",path_in+"EGU105053974_rooted.tre",path_in+"EGU105053980_rooted.tre",path_in+"EGU105054002_rooted.tre",path_in+"EGU105054124_rooted.tre",path_in+"EGU105054130_rooted.tre",path_in+"EGU105054153_rooted.tre",path_in+"EGU105054204_rooted.tre",path_in+"EGU105054280_rooted.tre",path_in+"EGU105054293_rooted.tre",path_in+"EGU105054405_rooted.tre",path_in+"EGU105054435_rooted.tre",path_in+"EGU105054440_rooted.tre",path_in+"EGU105054455_rooted.tre",path_in+"EGU105054457_rooted.tre",path_in+"EGU105054469_rooted.tre",path_in+"EGU105054478_rooted.tre",path_in+"EGU105054486_rooted.tre",path_in+"EGU105054498_rooted.tre",path_in+"EGU105054529_rooted.tre",path_in+"EGU105054534_rooted.tre",path_in+"EGU105054595_rooted.tre",path_in+"EGU105054649_rooted.tre",path_in+"EGU105054653_rooted.tre",path_in+"EGU105054668_rooted.tre",path_in+"EGU105054723_rooted.tre",path_in+"EGU105054765_rooted.tre",path_in+"EGU105054786_rooted.tre",path_in+"EGU105054827_rooted.tre",path_in+"EGU105054845_rooted.tre",path_in+"EGU105054864_rooted.tre",path_in+"EGU105054891_rooted.tre",path_in+"EGU105054896_rooted.tre",path_in+"EGU105054898_rooted.tre",path_in+"EGU105054924_rooted.tre",path_in+"EGU105054930_rooted.tre",path_in+"EGU105054936_rooted.tre",path_in+"EGU105054948_rooted.tre",path_in+"EGU105054972_rooted.tre",path_in+"EGU105055008_rooted.tre",path_in+"EGU105055015_rooted.tre",path_in+"EGU105055023_rooted.tre",path_in+"EGU105055024_rooted.tre",path_in+"EGU105055030_rooted.tre",path_in+"EGU105055047_rooted.tre",path_in+"EGU105055052_rooted.tre",path_in+"EGU105055065_rooted.tre",path_in+"EGU105055072_rooted.tre",path_in+"EGU105055075_rooted.tre",path_in+"EGU105055077_rooted.tre",path_in+"EGU105055090_rooted.tre",path_in+"EGU105055093_rooted.tre",path_in+"EGU105055098_rooted.tre",path_in+"EGU105055114_rooted.tre",path_in+"EGU105055115_rooted.tre",path_in+
	"EGU105055130_rooted.tre",path_in+"EGU105055144_rooted.tre",path_in+"EGU105055157_rooted.tre",path_in+"EGU105055201_rooted.tre",path_in+"EGU105055283_rooted.tre",path_in+"EGU105055433_rooted.tre",path_in+"EGU105055438_rooted.tre",path_in+"EGU105055490_rooted.tre",path_in+"EGU105055499_rooted.tre",path_in+"EGU105055507_rooted.tre",path_in+"EGU105055550_rooted.tre",path_in+"EGU105055569_rooted.tre",path_in+"EGU105055621_rooted.tre",path_in+"EGU105055634_rooted.tre",path_in+"EGU105055664_rooted.tre",path_in+"EGU105055709_rooted.tre",path_in+"EGU105055755_rooted.tre",path_in+"EGU105055761_rooted.tre",path_in+"EGU105055771_rooted.tre",path_in+"EGU105055800_rooted.tre",path_in+"EGU105055862_rooted.tre",path_in+"EGU105055873_rooted.tre",path_in+"EGU105055883_rooted.tre",path_in+"EGU105055889_rooted.tre",path_in+"EGU105055908_rooted.tre",path_in+"EGU105055912_rooted.tre",path_in+"EGU105055913_rooted.tre",path_in+"EGU105056032_rooted.tre",path_in+"EGU105056091_rooted.tre",path_in+"EGU105056151_rooted.tre",path_in+"EGU105056269_rooted.tre",path_in+"EGU105056287_rooted.tre",path_in+"EGU105056289_rooted.tre",path_in+"EGU105056313_rooted.tre",path_in+"EGU105056323_rooted.tre",path_in+"EGU105056365_rooted.tre",path_in+"EGU105056382_rooted.tre",path_in+"EGU105056393_rooted.tre",path_in+"EGU105056460_rooted.tre",path_in+"EGU105056468_rooted.tre",path_in+"EGU105056469_rooted.tre",path_in+"EGU105056496_rooted.tre",path_in+"EGU105056530_rooted.tre",path_in+"EGU105056534_rooted.tre",path_in+"EGU105056539_rooted.tre",path_in+"EGU105056654_rooted.tre",path_in+"EGU105056662_rooted.tre",path_in+"EGU105056684_rooted.tre",path_in+"EGU105056688_rooted.tre",path_in+"EGU105056714_rooted.tre",path_in+"EGU105056726_rooted.tre",path_in+"EGU105056817_rooted.tre",path_in+"EGU105056848_rooted.tre",path_in+"EGU105056881_rooted.tre",path_in+"EGU105056943_rooted.tre",path_in+"EGU105056960_rooted.tre",path_in+"EGU105056998_rooted.tre",path_in+"EGU105057013_rooted.tre",path_in+"EGU105057015_rooted.tre",path_in+"EGU105057019_rooted.tre",path_in+"EGU105057074_rooted.tre",path_in+"EGU105057090_rooted.tre",path_in+"EGU105057110_rooted.tre",path_in+"EGU105057130_rooted.tre",path_in+"EGU105057194_rooted.tre",path_in+"EGU105057235_rooted.tre",path_in+"EGU105057256_rooted.tre",path_in+"EGU105057335_rooted.tre",path_in+"EGU105057357_rooted.tre",path_in+"EGU105057553_rooted.tre",path_in+"EGU105057579_rooted.tre",path_in+"EGU105057634_rooted.tre",path_in+"EGU105057666_rooted.tre",path_in+"EGU105057669_rooted.tre",path_in+"EGU105057721_rooted.tre",path_in+"EGU105057742_rooted.tre",path_in+"EGU105057795_rooted.tre",path_in+"EGU105057841_rooted.tre",path_in+"EGU105057912_rooted.tre",path_in+"EGU105057919_rooted.tre",path_in+"EGU105057941_rooted.tre",path_in+"EGU105058078_rooted.tre",path_in+"EGU105058081_rooted.tre",path_in+"EGU105058083_rooted.tre",path_in+"EGU105058094_rooted.tre",path_in+"EGU105058107_rooted.tre",path_in+"EGU105058131_rooted.tre",path_in+"EGU105058170_rooted.tre",path_in+"EGU105058175_rooted.tre",path_in+"EGU105058180_rooted.tre",path_in+"EGU105058202_rooted.tre",path_in+"EGU105058237_rooted.tre",path_in+"EGU105058241_rooted.tre",path_in+"EGU105058245_rooted.tre",path_in+"EGU105058326_rooted.tre",path_in+"EGU105058366_rooted.tre",path_in+"EGU105058377_rooted.tre",path_in+"EGU105058418_rooted.tre",path_in+"EGU105058469_rooted.tre",path_in+"EGU105058499_rooted.tre",path_in+"EGU105058547_rooted.tre",path_in+"EGU105058556_rooted.tre",path_in+"EGU105058567_rooted.tre",path_in+"EGU105058576_rooted.tre",path_in+"EGU105058582_rooted.tre",path_in+"EGU105058592_rooted.tre",path_in+"EGU105058598_rooted.tre",path_in+"EGU105058614_rooted.tre",path_in+"EGU105058633_rooted.tre",path_in+"EGU105058682_rooted.tre",path_in+"EGU105058683_rooted.tre",path_in+"EGU105058687_rooted.tre",path_in+"EGU105058702_rooted.tre",path_in+"EGU105058723_rooted.tre",path_in+"EGU105058731_rooted.tre",path_in+"EGU105058781_rooted.tre",path_in+"EGU105058798_rooted.tre",path_in+"EGU105058802_rooted.tre",path_in+"EGU105058808_rooted.tre",path_in+"EGU105058863_rooted.tre",path_in+"EGU105058889_rooted.tre",path_in+"EGU105058890_rooted.tre",path_in+"EGU105058894_rooted.tre",path_in+"EGU105058904_rooted.tre",path_in+"EGU105058918_rooted.tre",path_in+"EGU105058989_rooted.tre",path_in+"EGU105058990_rooted.tre",path_in+"EGU105059003_rooted.tre",path_in+"EGU105059008_rooted.tre",path_in+"EGU105059023_rooted.tre",path_in+"EGU105059035_rooted.tre",path_in+"EGU105059042_rooted.tre",path_in+"EGU105059054_rooted.tre",path_in+"EGU105059108_rooted.tre",path_in+"EGU105059112_rooted.tre",path_in+"EGU105059113_rooted.tre",path_in+"EGU105059126_rooted.tre",path_in+"EGU105059131_rooted.tre",path_in+"EGU105059138_rooted.tre",path_in+"EGU105059176_rooted.tre",path_in+"EGU105059186_rooted.tre",path_in+"EGU105059193_rooted.tre",path_in+"EGU105059276_rooted.tre",path_in+"EGU105059342_rooted.tre",path_in+"EGU105059366_rooted.tre",path_in+"EGU105059367_rooted.tre",path_in+"EGU105059381_rooted.tre",path_in+"EGU105059441_rooted.tre",path_in+"EGU105059453_rooted.tre",path_in+"EGU105059458_rooted.tre",path_in+"EGU105059479_rooted.tre",path_in+"EGU105059480_rooted.tre",path_in+"EGU105059490_rooted.tre",path_in+"EGU105059570_rooted.tre",path_in+"EGU105059573_rooted.tre",path_in+"EGU105059575_rooted.tre",path_in+"EGU105059587_rooted.tre",path_in+"EGU105059612_rooted.tre",path_in+"EGU105059624_rooted.tre",path_in+"EGU105059636_rooted.tre",path_in+"EGU105059639_rooted.tre",path_in+"EGU105059671_rooted.tre",path_in+"EGU105059853_rooted.tre",path_in+"EGU105059900_rooted.tre",path_in+"EGU105060095_rooted.tre",path_in+"EGU105060589_rooted.tre",path_in+"EGU105061025_rooted.tre",path_in+"EGU105061385_rooted.tre",path_in+"EGU105061427_rooted.tre",path_in+"HEY1007_rooted.tre",path_in+"HEY1013_rooted.tre",path_in+"HEY1017_rooted.tre",path_in+"HEY1020_rooted.tre",path_in+"HEY1025_rooted.tre",path_in+"HEY1035_rooted.tre",path_in+"HEY1050_rooted.tre",path_in+"HEY1052_rooted.tre",path_in+"HEY1064_rooted.tre",path_in+"HEY110_rooted.tre",path_in+"HEY1119_rooted.tre",path_in+"HEY1168_rooted.tre",path_in+"HEY1171_rooted.tre",path_in+"HEY1197_rooted.tre",path_in+"HEY1201_rooted.tre",path_in+"HEY120_rooted.tre",path_in+"HEY122_rooted.tre",path_in+"HEY125_rooted.tre",path_in+"HEY12_rooted.tre",path_in+"HEY136_rooted.tre",path_in+"HEY139_rooted.tre",path_in+"HEY1484_rooted.tre",path_in+"HEY148_rooted.tre",path_in+"HEY1494_rooted.tre",path_in+"HEY14_rooted.tre",path_in+"HEY150_rooted.tre",path_in+"HEY1615_rooted.tre",path_in+"HEY164_rooted.tre",path_in+"HEY17_rooted.tre",path_in+"HEY1801_rooted.tre",path_in+"HEY180_rooted.tre",path_in+"HEY1815_rooted.tre",path_in+"HEY182_rooted.tre",path_in+"HEY1842_rooted.tre",path_in+"HEY1854_rooted.tre",path_in+"HEY1877_rooted.tre",path_in+"HEY1901_rooted.tre",path_in+"HEY191_rooted.tre",path_in+"HEY194_rooted.tre",path_in+"HEY197_rooted.tre",path_in+"HEY1986_rooted.tre",path_in+"HEY201_rooted.tre",path_in+"HEY204e_rooted.tre",path_in+"HEY204s_rooted.tre",path_in+"HEY2056_rooted.tre",path_in+"HEY207_rooted.tre",path_in+"HEY215_rooted.tre",path_in+"HEY2164_rooted.tre",path_in+"HEY218_rooted.tre",path_in+"HEY21_rooted.tre",path_in+"HEY2238_rooted.tre",path_in+"HEY225_rooted.tre",path_in+"HEY226_rooted.tre",path_in+"HEY2291_rooted.tre",path_in+"HEY231_rooted.tre",path_in+"HEY2339_rooted.tre",path_in+"HEY2363_rooted.tre",path_in+"HEY2370_rooted.tre",path_in+"HEY2377_rooted.tre",path_in+"HEY237_rooted.tre",path_in+"HEY2388_rooted.tre",path_in+"HEY240_rooted.tre",path_in+"HEY2459_rooted.tre",path_in+"HEY245_rooted.tre",path_in+"HEY24_rooted.tre",path_in+"HEY250_rooted.tre",path_in+"HEY252e_rooted.tre",path_in+"HEY252p_rooted.tre",path_in+"HEY252s_rooted.tre",path_in+"HEY2550_rooted.tre",path_in+"HEY2561_rooted.tre",path_in+"HEY257_rooted.tre",path_in+"HEY267_rooted.tre",path_in+"HEY269_rooted.tre",path_in+"HEY277_rooted.tre",path_in+"HEY280_rooted.tre",path_in+"HEY281_rooted.tre",path_in+"HEY282_rooted.tre",path_in+"HEY290_rooted.tre",path_in+"HEY293_rooted.tre",path_in+"HEY296_rooted.tre",path_in+"HEY299_rooted.tre",path_in+"HEY305_rooted.tre",path_in+"HEY308_rooted.tre",path_in+"HEY310_rooted.tre",path_in+"HEY31_rooted.tre",path_in+"HEY323_rooted.tre",path_in+"HEY326_rooted.tre",path_in+"HEY32e_rooted.tre",path_in+"HEY32s_rooted.tre",path_in+"HEY332_rooted.tre",path_in+"HEY334_rooted.tre",path_in+"HEY340_rooted.tre",path_in+"HEY357_rooted.tre",path_in+"HEY360_rooted.tre",path_in+"HEY362_rooted.tre",path_in+"HEY363_rooted.tre",path_in+"HEY369_rooted.tre",path_in+"HEY378e_rooted.tre",path_in+"HEY378s_rooted.tre",path_in+"HEY38_rooted.tre",path_in+"HEY391_rooted.tre",path_in+"HEY392_rooted.tre",path_in+"HEY415_rooted.tre",path_in+"HEY417_rooted.tre",path_in+"HEY421_rooted.tre",path_in+"HEY429_rooted.tre",path_in+"HEY449_rooted.tre",path_in+"HEY464_rooted.tre",path_in+"HEY484_rooted.tre",path_in+"HEY490_rooted.tre",path_in+"HEY497_rooted.tre",path_in+"HEY4_rooted.tre",path_in+"HEY508_rooted.tre",path_in+"HEY514_rooted.tre",path_in+"HEY51_rooted.tre",path_in+"HEY52_rooted.tre",path_in+"HEY556_rooted.tre",path_in+"HEY563_rooted.tre",path_in+"HEY576_rooted.tre",path_in+"HEY581_rooted.tre",path_in+"HEY587_rooted.tre",path_in+"HEY604_rooted.tre",path_in+"HEY609_rooted.tre",path_in+"HEY61_rooted.tre",path_in+"HEY629_rooted.tre",path_in+"HEY630_rooted.tre",path_in+"HEY637_rooted.tre",path_in+"HEY673_rooted.tre",path_in+"HEY680_rooted.tre",path_in+"HEY703_rooted.tre",path_in+"HEY717_rooted.tre",path_in+"HEY727_rooted.tre",path_in+"HEY728_rooted.tre",path_in+"HEY732_rooted.tre",path_in+"HEY736_rooted.tre",path_in+"HEY740_rooted.tre",path_in+"HEY743_rooted.tre",path_in+"HEY757_rooted.tre",path_in+"HEY758_rooted.tre",path_in+"HEY762_rooted.tre",path_in+"HEY763_rooted.tre",path_in+"HEY785_rooted.tre",path_in+"HEY790_rooted.tre",path_in+"HEY793_rooted.tre",path_in+"HEY7_rooted.tre",path_in+"HEY807_rooted.tre",path_in+"HEY808_rooted.tre",path_in+"HEY822_rooted.tre",path_in+"HEY825_rooted.tre",path_in+"HEY82_rooted.tre",path_in+"HEY83_rooted.tre",path_in+"HEY84_rooted.tre",path_in+"HEY855_rooted.tre",path_in+"HEY856_rooted.tre",path_in+"HEY863_rooted.tre",path_in+"HEY872_rooted.tre",path_in+"HEY874_rooted.tre",path_in+"HEY883e_rooted.tre",path_in+"HEY883n_rooted.tre",path_in+"HEY886_rooted.tre",path_in+"HEY88_rooted.tre",path_in+"HEY897_rooted.tre",path_in+"HEY89_rooted.tre",path_in+"HEY938_rooted.tre",path_in+"HEY948_rooted.tre",path_in+"HEY94_rooted.tre",path_in+"HEY950_rooted.tre",path_in+"HEY958_rooted.tre",path_in+"HEY964_rooted.tre",path_in+"HEY977_rooted.tre",path_in+"HEY982_rooted.tre",path_in+
	"HEY985_rooted.tre",path_in+"HEY989_rooted.tre"]
    outputs = [path_out+"genetrees.tre"]
    options = {'cores': 2, 'memory': "10g", 'walltime': "00:10:00", 'account':"Coryphoideae"}

    spec = """
	source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
	conda activate treebuilder_env

	echo "this is the path in: {path_in}"
	echo "this is the path out: {path_out}"

	cd {path_in}


	for f in *_rooted.tre
	do 
		nw_ed $f 'i & (b<30)' o >> {path_out}genetrees.tre #Moves trees used in treebuilding

	done

	""".format(path_in = path_in, path_out=path_out)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)



# ########################################################################################################################
# ############################---- Newick Contracting and gathering of genetrees ----#####################################
# ###############################---- Using only genes which are not paralogs ----########################################
# ########################################################################################################################

def newick_contracting_orthologs(path_in,path_out ):
    """Gathering all the gene-trees in a single file and removing branches with low bootstrap support
	I have had to remove the following genes for varius reasons such as outgroup species missing and low gene recovery and bad optrimal values
	
	I should have Included HEY168 but I have no recovery of this gene in any of my samples."""

    inputs = [path_in+"ortholog_genes.txt",path_in+"EGU105035046_rooted.tre",path_in+"EGU105035203_rooted.tre",path_in+"EGU105035555_rooted.tre",path_in+"EGU105035989_rooted.tre",path_in+"EGU105036774_rooted.tre",path_in+"EGU105038201_rooted.tre",path_in+"EGU105038431_rooted.tre",path_in+"EGU105038513_rooted.tre",path_in+"EGU105038747_rooted.tre",path_in+"EGU105038832_rooted.tre",path_in+"EGU105039082_rooted.tre",path_in+"EGU105039099_rooted.tre",path_in+"EGU105039164_rooted.tre",path_in+"EGU105039255_rooted.tre",path_in+"EGU105039449_rooted.tre",path_in+"EGU105039494_rooted.tre",path_in+"EGU105039501_rooted.tre",path_in+"EGU105039512_rooted.tre",path_in+"EGU105039542_rooted.tre",path_in+"EGU105039660_rooted.tre",path_in+"EGU105039690_rooted.tre",path_in+"EGU105039783_rooted.tre",path_in+"EGU105040185_rooted.tre",path_in+"EGU105040242_rooted.tre",path_in+"EGU105040462_rooted.tre",path_in+"EGU105040842_rooted.tre",path_in+"EGU105040851_rooted.tre",path_in+"EGU105040914_rooted.tre",path_in+"EGU105040918_rooted.tre",path_in+"EGU105040970_rooted.tre",path_in+"EGU105041100_rooted.tre",path_in+"EGU105041179_rooted.tre",path_in+"EGU105041217_rooted.tre",path_in+"EGU105041337_rooted.tre",path_in+"EGU105041650_rooted.tre",path_in+"EGU105041903_rooted.tre",path_in+"EGU105041933_rooted.tre",path_in+"EGU105042090_rooted.tre",path_in+"EGU105042633_rooted.tre",path_in+"EGU105042781_rooted.tre",path_in+"EGU105043037_rooted.tre",path_in+"EGU105043469_rooted.tre",path_in+"EGU105043485_rooted.tre",path_in+"EGU105043633_rooted.tre",path_in+"EGU105043686_rooted.tre",path_in+"EGU105043786_rooted.tre",path_in+"EGU105044407_rooted.tre",path_in+"EGU105044668_rooted.tre",path_in+"EGU105044710_rooted.tre",path_in+"EGU105045078_rooted.tre",path_in+"EGU105045099_rooted.tre",path_in+"EGU105045282_rooted.tre",path_in+"EGU105045367_rooted.tre",path_in+"EGU105045467_rooted.tre",path_in+"EGU105046245_rooted.tre",path_in+"EGU105046456_rooted.tre",path_in+"EGU105047379_rooted.tre",path_in+"EGU105047434_rooted.tre",path_in+"EGU105047446_rooted.tre",path_in+"EGU105047553_rooted.tre",path_in+"EGU105047621_rooted.tre",path_in+"EGU105047790_rooted.tre",path_in+"EGU105047940_rooted.tre",path_in+"EGU105048476_rooted.tre",path_in+"EGU105048541_rooted.tre",path_in+"EGU105048694_rooted.tre",path_in+"EGU105048909_rooted.tre",path_in+"EGU105048961_rooted.tre",path_in+"EGU105049052_rooted.tre",path_in+"EGU105049312_rooted.tre",path_in+"EGU105049360_rooted.tre",path_in+"EGU105049539_rooted.tre",path_in+"EGU105050126_rooted.tre",path_in+"EGU105050344_rooted.tre",path_in+"EGU105050366_rooted.tre",path_in+"EGU105050383_rooted.tre",path_in+"EGU105050521_rooted.tre",path_in+"EGU105050682_rooted.tre",path_in+"EGU105050853_rooted.tre",path_in+"EGU105051156_rooted.tre",path_in+"EGU105051188_rooted.tre",path_in+"EGU105051362_rooted.tre",path_in+"EGU105051403_rooted.tre",path_in+"EGU105051560_rooted.tre",path_in+"EGU105051726_rooted.tre",path_in+"EGU105051764_rooted.tre",path_in+"EGU105051795_rooted.tre",path_in+"EGU105051847_rooted.tre",path_in+"EGU105051860_rooted.tre",path_in+"EGU105051870_rooted.tre",path_in+"EGU105051985_rooted.tre",path_in+"EGU105052307_rooted.tre",path_in+"EGU105052476_rooted.tre",path_in+"EGU105052492_rooted.tre",path_in+"EGU105052580_rooted.tre",path_in+"EGU105052804_rooted.tre",path_in+"EGU105052818_rooted.tre",path_in+"EGU105052855_rooted.tre",path_in+"EGU105052888_rooted.tre",path_in+"EGU105052956_rooted.tre",path_in+"EGU105053006_rooted.tre",path_in+"EGU105053055_rooted.tre",path_in+"EGU105053079_rooted.tre",path_in+"EGU105053422_rooted.tre",path_in+"EGU105053482_rooted.tre",path_in+"EGU105053549_rooted.tre",path_in+"EGU105053848_rooted.tre",path_in+"EGU105053866_rooted.tre",path_in+"EGU105053889_rooted.tre",path_in+"EGU105053901_rooted.tre",path_in+"EGU105053980_rooted.tre",path_in+"EGU105054153_rooted.tre",path_in+"EGU105054405_rooted.tre",path_in+"EGU105054595_rooted.tre",path_in+"EGU105054786_rooted.tre",path_in+"EGU105054930_rooted.tre",path_in+"EGU105054948_rooted.tre",path_in+"EGU105055065_rooted.tre",path_in+"EGU105055072_rooted.tre",path_in+"EGU105055075_rooted.tre",path_in+"EGU105055114_rooted.tre",path_in+"EGU105055115_rooted.tre",path_in+"EGU105055499_rooted.tre",path_in+"EGU105055664_rooted.tre",path_in+"EGU105055800_rooted.tre",path_in+"EGU105055873_rooted.tre",path_in+"EGU105056091_rooted.tre",path_in+"EGU105056289_rooted.tre",path_in+"EGU105056469_rooted.tre",path_in+"EGU105056654_rooted.tre",path_in+"EGU105057074_rooted.tre",path_in+"EGU105057634_rooted.tre",path_in+"EGU105057666_rooted.tre",path_in+"EGU105057721_rooted.tre",path_in+"EGU105058081_rooted.tre",path_in+"EGU105058131_rooted.tre",path_in+"EGU105058170_rooted.tre",path_in+"EGU105058180_rooted.tre",path_in+"EGU105058237_rooted.tre",path_in+"EGU105058469_rooted.tre",path_in+"EGU105058702_rooted.tre",path_in+"EGU105058798_rooted.tre",path_in+"EGU105058889_rooted.tre",path_in+"EGU105059023_rooted.tre",path_in+"EGU105059126_rooted.tre",path_in+"EGU105059138_rooted.tre",path_in+"EGU105059186_rooted.tre",path_in+"EGU105059366_rooted.tre",path_in+"EGU105059381_rooted.tre",path_in+"EGU105059480_rooted.tre",path_in+"EGU105059624_rooted.tre",path_in+"EGU105059639_rooted.tre",path_in+"EGU105059853_rooted.tre",path_in+"EGU105061427_rooted.tre",path_in+"HEY1007_rooted.tre",path_in+"HEY1017_rooted.tre",path_in+"HEY1020_rooted.tre",path_in+"HEY1035_rooted.tre",path_in+"HEY1052_rooted.tre",path_in+"HEY1064_rooted.tre",path_in+"HEY1119_rooted.tre",path_in+"HEY1197_rooted.tre",path_in+"HEY1201_rooted.tre",path_in+"HEY120_rooted.tre",path_in+"HEY122_rooted.tre",path_in+"HEY12_rooted.tre",path_in+"HEY150_rooted.tre",path_in+"HEY1615_rooted.tre",path_in+"HEY1815_rooted.tre",path_in+"HEY182_rooted.tre",path_in+"HEY1854_rooted.tre",path_in+"HEY194_rooted.tre",path_in+"HEY197_rooted.tre",path_in+"HEY201_rooted.tre",path_in+"HEY204e_rooted.tre",path_in+"HEY2056_rooted.tre",path_in+"HEY2238_rooted.tre",path_in+"HEY226_rooted.tre",path_in+"HEY231_rooted.tre",path_in+"HEY250_rooted.tre",path_in+"HEY252e_rooted.tre",path_in+"HEY252p_rooted.tre",path_in+"HEY252s_rooted.tre",path_in+"HEY2550_rooted.tre",path_in+"HEY2561_rooted.tre",path_in+"HEY257_rooted.tre",path_in+"HEY269_rooted.tre",path_in+"HEY277_rooted.tre",path_in+"HEY281_rooted.tre",path_in+"HEY282_rooted.tre",path_in+"HEY290_rooted.tre",path_in+"HEY293_rooted.tre",path_in+"HEY299_rooted.tre",path_in+"HEY305_rooted.tre",path_in+"HEY31_rooted.tre",path_in+"HEY326_rooted.tre",path_in+"HEY32e_rooted.tre",path_in+"HEY32s_rooted.tre",path_in+"HEY340_rooted.tre",path_in+"HEY357_rooted.tre",path_in+"HEY362_rooted.tre",path_in+"HEY363_rooted.tre",path_in+"HEY369_rooted.tre",path_in+"HEY378s_rooted.tre",path_in+"HEY38_rooted.tre",path_in+"HEY392_rooted.tre",path_in+"HEY51_rooted.tre",path_in+"HEY576_rooted.tre",path_in+"HEY587_rooted.tre",path_in+"HEY604_rooted.tre",path_in+"HEY61_rooted.tre",path_in+"HEY629_rooted.tre",path_in+"HEY630_rooted.tre",path_in+"HEY637_rooted.tre",path_in+"HEY703_rooted.tre",path_in+"HEY740_rooted.tre",path_in+"HEY758_rooted.tre",path_in+"HEY790_rooted.tre",path_in+"HEY7_rooted.tre",path_in+"HEY807_rooted.tre",path_in+"HEY83_rooted.tre",path_in+"HEY855_rooted.tre",path_in+"HEY863_rooted.tre",path_in+"HEY872_rooted.tre",path_in+"HEY886_rooted.tre",path_in+"HEY88_rooted.tre",path_in+"HEY948_rooted.tre",path_in+"HEY94_rooted.tre",path_in+"HEY950_rooted.tre",path_in+"HEY977_rooted.tre",path_in+"HEY985_rooted.tre"]

    outputs = [path_out+"genetrees_orthologs.tre"]
    options = {'cores': 2, 'memory': "10g", 'walltime': "00:10:00", 'account':"Coryphoideae"}

    spec = """
	source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
	conda activate treebuilder_env

	echo "this is the path in: {path_in}"
	echo "this is the path out: {path_out}"

	cd {path_in}

	while IFS= read -r line; do
		echo "$line"
		nw_ed $line 'i & (b<30)' o >> {path_out}genetrees_orthologs.tre #Moves trees used in treebuilding
	done < ortholog_genes.txt


	""".format(path_in = path_in, path_out=path_out)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)





# ########################################################################################################################
# #####################################---- Astral Tree Search ----#####################################################
# ########################################################################################################################
def astral(path_in, gene_tree_file,output):
    """Using Astral to construct a species tree based on the genetrees"""
    inputs = [path_in+gene_tree_file]
    outputs = [path_in+output]
    options = {'cores': 20, 'memory': "40g", 'walltime': "1:00:00", 'account':"Coryphoideae"}

    spec = """
	source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
	conda activate treebuilder_env

	cd {path_in}

	java -jar /home/owrisberg/Coryphoideae/github_code/ASTRAL/astral.5.7.7.jar -i {gene_tree_file} -o {output}


	""".format(path_in = path_in, gene_tree_file = gene_tree_file, output=output)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# ########################################################################################################################
# ############################################---- Renaming ----#####################################################
# ########################################################################################################################
def renaming(path_in, tree_in,gene_tree_file, tree_out):
    """Renaming the tips in the phylogeny based on the names_for_tips.csv"""
    inputs = [path_in+tree_in]
    outputs = [path_in+tree_out]
    options = {'cores': 1, 'memory': "10g", 'walltime': "00:10:00", 'account':"Coryphoideae"}

    spec = """
	source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
	conda activate treebuilder_env

	cd {path_in}

	#Renaming tips in tree
	python3 /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/renamer.py /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/names_for_tips.csv {tree_in} {tree_out} --bs 1

	""".format(path_in = path_in, tree_in=tree_in, tree_out=tree_out, gene_tree_file = gene_tree_file)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# ########################################################################################################################
# ############################################---- Quartetscores ----#####################################################
# ########################################################################################################################
def quartet_scores(path_in):
    """Using Astral to construct a species tree based on the genetrees"""
    inputs = [path_in+"genetrees.tre", path_in+"astral_tree_renamed.tre"]
    outputs = [path_in+"astral_tree_QS_renamed.tre"]
    options = {'cores': 20, 'memory': "40g", 'walltime': "10:30:00", 'account':"Coryphoideae"}

    spec = """
	source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
	conda activate treebuilder_env

	cd {path_in}

	#Running quartet scores
    /home/owrisberg/Coryphoideae/github_code/QuartetScores -o astral_tree_QS.tre -e genetrees.tre -r astral_tree.tre -v

	#Correcting labels from Quartetscores
	sed astral_tree_QS.tre -i'.old' -e s/[0-9]\.*[0-9]*\(:[0-9]\.*[0-9]*\)\[qp-ic:-*[0-9]\.[0-9]*;lq-ic:-*[0-9]\.[0-9]*;eqp-ic:\(-*[0-9]\.[0-9]*\)\]/\2\1/g`
	sed astral_tree_QS.tre -i'.old' -e 's/\[eqp-ic:-*[0-9]\.*[0-9]*\]//g`


	""".format(path_in = path_in)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# ########################################################################################################################
# #####################################---- Astral Tree annotation ----#####################################################
# ########################################################################################################################
def astral_annotation(path_in, gene_tree_file, species_tree_file, outfile):
    """Using Astral to construct a species tree based on the genetrees"""
    inputs = [path_in+gene_tree_file, path_in+species_tree_file]
    outputs = [path_in+outfile]
    options = {'cores': 20, 'memory': "40g", 'walltime': "04:00:00", 'account':"Coryphoideae"}

    spec = """
	source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
	conda activate treebuilder_env

	cd {path_in}

	java -jar /home/owrisberg/Coryphoideae/github_code/ASTRAL/astral.5.7.7.jar -q {species_tree_file} -i {gene_tree_file} -t 2 -o {outfile}


	""".format(path_in = path_in, gene_tree_file = gene_tree_file, species_tree_file = species_tree_file, outfile=outfile)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# ########################################################################################################################
# #####################################---- SortaDate ----#####################################################
# ########################################################################################################################
def sorta_date(path_in,path_out,astral_tree, done):
    """Using SortaDate to produce a CSV file which can be used to evaluate the use of different genes in dating the trees"""
    inputs = [astral_tree]
    outputs = [path_out+"var",path_out+"bp",path_out+"comb", path_out+"gg", done]
    options = {'cores': 3, 'memory': "10g", 'walltime': "00:10:00", 'account':"Coryphoideae"}

    spec = """
	#Activating conda base environment 
	source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
	conda activate base

	#Get the root-to-tip variance with
	python {SortaDate}get_var_length.py {path_in} --flend _rooted.tre --outf {path_out}var --outg 1079,1080,1081,1082

	#Get the bipartition support with
	python {SortaDate}get_bp_genetrees.py {path_in} {astral_tree} --flend _rooted.tre --outf {path_out}bp

	#Combine the results from these two runs with
	python {SortaDate}combine_results.py {path_out}var {path_out}bp --outf {path_out}comb

	#Sort and get the list of the good genes with
	python {SortaDate}get_good_genes.py {path_out}comb --max 1000 --order 3,1,2 --outf {path_out}gg

	touch {done}

	""".format(path_in = path_in,path_out = path_out, SortaDate = "/home/owrisberg/Coryphoideae/github_code/SortaDate/src/", astral_tree = astral_tree, done = done)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)





########################################################################################################################
######################################################---- RUN ----#####################################################
########################################################################################################################


genes = ["EGU105032175","EGU105032229","EGU105032337","EGU105032379","EGU105033063","EGU105033626","EGU105034121","EGU105034616","EGU105034893","EGU105034993","EGU105035046","EGU105035196","EGU105035203","EGU105035462","EGU105035555","EGU105035989","EGU105036031","EGU105036385","EGU105036774","EGU105037749","EGU105037800","EGU105037890","EGU105037902","EGU105037930","EGU105037938","EGU105038008","EGU105038036","EGU105038098","EGU105038099","EGU105038100","EGU105038110","EGU105038114","EGU105038118","EGU105038123","EGU105038179","EGU105038201","EGU105038228","EGU105038234","EGU105038245","EGU105038252","EGU105038310","EGU105038382","EGU105038400","EGU105038419","EGU105038431","EGU105038499","EGU105038513","EGU105038571","EGU105038580","EGU105038603","EGU105038631","EGU105038680","EGU105038693","EGU105038720","EGU105038747","EGU105038794","EGU105038832","EGU105038882","EGU105038986","EGU105038988","EGU105039013","EGU105039062","EGU105039067","EGU105039082","EGU105039099","EGU105039101","EGU105039107","EGU105039121","EGU105039164","EGU105039178","EGU105039221","EGU105039236","EGU105039255","EGU105039282","EGU105039298","EGU105039313","EGU105039403","EGU105039431","EGU105039449","EGU105039460","EGU105039480","EGU105039494","EGU105039501","EGU105039512","EGU105039542","EGU105039587","EGU105039595","EGU105039609","EGU105039660","EGU105039685","EGU105039690","EGU105039699","EGU105039763","EGU105039783","EGU105039809","EGU105039822","EGU105039925","EGU105039947","EGU105039957","EGU105039962","EGU105040073","EGU105040088","EGU105040099","EGU105040114","EGU105040115","EGU105040125","EGU105040139","EGU105040185","EGU105040186","EGU105040189","EGU105040206","EGU105040207","EGU105040242","EGU105040281","EGU105040302","EGU105040308","EGU105040359","EGU105040368","EGU105040426","EGU105040452","EGU105040462","EGU105040530","EGU105040583","EGU105040667","EGU105040675","EGU105040684","EGU105040690","EGU105040700","EGU105040756","EGU105040758","EGU105040813","EGU105040837","EGU105040842","EGU105040850","EGU105040851","EGU105040863","EGU105040887","EGU105040914","EGU105040918","EGU105040922","EGU105040957","EGU105040970","EGU105041055","EGU105041100","EGU105041117","EGU105041125","EGU105041127","EGU105041133","EGU105041179","EGU105041182","EGU105041189","EGU105041217","EGU105041246","EGU105041283","EGU105041337","EGU105041353","EGU105041650","EGU105041657","EGU105041665","EGU105041680","EGU105041687","EGU105041710","EGU105041807","EGU105041816","EGU105041872","EGU105041902","EGU105041903","EGU105041929","EGU105041933","EGU105041982","EGU105042090","EGU105042113","EGU105042128","EGU105042147","EGU105042168","EGU105042205","EGU105042290","EGU105042307","EGU105042323","EGU105042329","EGU105042368","EGU105042422","EGU105042525","EGU105042558","EGU105042560","EGU105042584","EGU105042633","EGU105042644","EGU105042651","EGU105042664","EGU105042722","EGU105042781","EGU105042808","EGU105042820","EGU105042873","EGU105042965","EGU105043011","EGU105043037","EGU105043042","EGU105043061","EGU105043069","EGU105043119","EGU105043155","EGU105043160","EGU105043164","EGU105043193","EGU105043320","EGU105043338","EGU105043374","EGU105043419","EGU105043430","EGU105043469","EGU105043485","EGU105043499","EGU105043601","EGU105043633","EGU105043666","EGU105043685","EGU105043686","EGU105043730","EGU105043786","EGU105043816","EGU105043827","EGU105043926","EGU105043975","EGU105044063","EGU105044120","EGU105044133","EGU105044174","EGU105044182","EGU105044183","EGU105044203","EGU105044252","EGU105044281","EGU105044307","EGU105044309","EGU105044350","EGU105044378","EGU105044400","EGU105044407","EGU105044445","EGU105044446","EGU105044481","EGU105044588","EGU105044613","EGU105044614","EGU105044668","EGU105044676","EGU105044710","EGU105044758","EGU105044844","EGU105044846","EGU105044854","EGU105044885","EGU105044893","EGU105044896","EGU105044978","EGU105044982","EGU105044983","EGU105044984","EGU105045005","EGU105045043","EGU105045070","EGU105045078","EGU105045094","EGU105045099","EGU105045102","EGU105045137","EGU105045148","EGU105045232","EGU105045248","EGU105045254","EGU105045282","EGU105045310","EGU105045358","EGU105045367","EGU105045424","EGU105045464","EGU105045467","EGU105045489","EGU105045507","EGU105045509","EGU105045514","EGU105045520","EGU105045529","EGU105045544","EGU105045640","EGU105045658","EGU105045703","EGU105045726","EGU105045732","EGU105045760","EGU105045782","EGU105045788","EGU105045820","EGU105045827","EGU105045828","EGU105045835","EGU105045898","EGU105045932","EGU105045946","EGU105046030","EGU105046050","EGU105046056","EGU105046099","EGU105046103","EGU105046147","EGU105046168","EGU105046245","EGU105046297","EGU105046360","EGU105046387","EGU105046393","EGU105046401","EGU105046449","EGU105046454","EGU105046456","EGU105046503","EGU105046518","EGU105046530","EGU105046549","EGU105046559","EGU105046562","EGU105046574","EGU105046630","EGU105046632","EGU105046696","EGU105046735","EGU105046766","EGU105046786","EGU105046827","EGU105046875","EGU105046918","EGU105047024","EGU105047029","EGU105047253","EGU105047288","EGU105047293","EGU105047342","EGU105047357","EGU105047362","EGU105047379","EGU105047385","EGU105047395","EGU105047433","EGU105047434","EGU105047446","EGU105047519","EGU105047533","EGU105047546","EGU105047553","EGU105047578","EGU105047585","EGU105047597","EGU105047621","EGU105047644","EGU105047662","EGU105047689","EGU105047751","EGU105047777","EGU105047790","EGU105047907","EGU105047916","EGU105047922","EGU105047940","EGU105047945","EGU105047970","EGU105048009","EGU105048015","EGU105048028","EGU105048054","EGU105048056","EGU105048129","EGU105048130","EGU105048137","EGU105048159","EGU105048182","EGU105048199","EGU105048300","EGU105048357","EGU105048410","EGU105048474","EGU105048476","EGU105048479","EGU105048484","EGU105048486","EGU105048493","EGU105048527","EGU105048541","EGU105048581","EGU105048612","EGU105048694","EGU105048725","EGU105048751","EGU105048796","EGU105048839","EGU105048867","EGU105048886","EGU105048898","EGU105048909","EGU105048915","EGU105048926","EGU105048961","EGU105048968","EGU105049007","EGU105049016","EGU105049020","EGU105049025","EGU105049052","EGU105049097","EGU105049274","EGU105049312","EGU105049318","EGU105049360","EGU105049426","EGU105049539","EGU105049543","EGU105049583","EGU105049690","EGU105049729","EGU105049737","EGU105049761","EGU105049827","EGU105049882","EGU105049902","EGU105049903","EGU105049934","EGU105049947","EGU105050012","EGU105050023","EGU105050036","EGU105050058","EGU105050114","EGU105050126","EGU105050202","EGU105050207","EGU105050328","EGU105050344","EGU105050362","EGU105050366","EGU105050383","EGU105050387","EGU105050404","EGU105050432","EGU105050450","EGU105050521","EGU105050532","EGU105050644","EGU105050670","EGU105050680","EGU105050681","EGU105050682","EGU105050831","EGU105050841","EGU105050853","EGU105050854","EGU105050961","EGU105050970","EGU105050972","EGU105051087","EGU105051146","EGU105051156","EGU105051188","EGU105051345","EGU105051362","EGU105051366","EGU105051373","EGU105051391","EGU105051395","EGU105051403","EGU105051481","EGU105051499","EGU105051503","EGU105051560","EGU105051564","EGU105051582","EGU105051614","EGU105051677","EGU105051704","EGU105051726","EGU105051740","EGU105051748","EGU105051764","EGU105051795","EGU105051802","EGU105051821","EGU105051823","EGU105051832","EGU105051847","EGU105051857","EGU105051860","EGU105051870","EGU105051891","EGU105051924","EGU105051953","EGU105051985","EGU105052035","EGU105052070","EGU105052170","EGU105052178","EGU105052304","EGU105052307","EGU105052346","EGU105052351","EGU105052386","EGU105052389","EGU105052394","EGU105052428","EGU105052446","EGU105052476","EGU105052483","EGU105052492","EGU105052495","EGU105052527","EGU105052529","EGU105052538","EGU105052552","EGU105052573","EGU105052580","EGU105052623","EGU105052650","EGU105052694","EGU105052704","EGU105052739","EGU105052743","EGU105052750","EGU105052771","EGU105052804","EGU105052818","EGU105052849","EGU105052855","EGU105052865","EGU105052888","EGU105052944","EGU105052947","EGU105052956","EGU105053006","EGU105053055","EGU105053059","EGU105053079","EGU105053105","EGU105053124","EGU105053130","EGU105053136","EGU105053172","EGU105053204","EGU105053227","EGU105053263","EGU105053403","EGU105053422","EGU105053426","EGU105053457","EGU105053465","EGU105053468","EGU105053482","EGU105053549","EGU105053642","EGU105053654","EGU105053735","EGU105053747","EGU105053770","EGU105053835","EGU105053848","EGU105053866","EGU105053889","EGU105053901","EGU105053932","EGU105053961","EGU105053969","EGU105053974","EGU105053980","EGU105054002","EGU105054124","EGU105054130","EGU105054153","EGU105054204","EGU105054280","EGU105054293","EGU105054405","EGU105054435","EGU105054440","EGU105054455","EGU105054457","EGU105054469","EGU105054478","EGU105054486","EGU105054498","EGU105054529","EGU105054534","EGU105054595","EGU105054649","EGU105054653","EGU105054668","EGU105054723","EGU105054765","EGU105054786","EGU105054827","EGU105054845","EGU105054864","EGU105054891","EGU105054896","EGU105054898","EGU105054924","EGU105054930","EGU105054936","EGU105054948","EGU105054972","EGU105055008","EGU105055015","EGU105055023","EGU105055024","EGU105055030","EGU105055047","EGU105055052","EGU105055065","EGU105055072","EGU105055075","EGU105055077","EGU105055090","EGU105055093","EGU105055098","EGU105055114","EGU105055115","EGU105055130","EGU105055144","EGU105055157","EGU105055201","EGU105055283","EGU105055433","EGU105055438","EGU105055490","EGU105055499","EGU105055507","EGU105055550","EGU105055569","EGU105055621","EGU105055634","EGU105055664","EGU105055709","EGU105055755","EGU105055761","EGU105055771","EGU105055800","EGU105055862","EGU105055873","EGU105055883","EGU105055889","EGU105055908","EGU105055912","EGU105055913","EGU105056032","EGU105056091","EGU105056151","EGU105056269","EGU105056287","EGU105056289","EGU105056313","EGU105056323","EGU105056365","EGU105056382","EGU105056393","EGU105056460","EGU105056468","EGU105056469","EGU105056496","EGU105056530","EGU105056534","EGU105056539","EGU105056654","EGU105056662","EGU105056684","EGU105056688","EGU105056714","EGU105056726","EGU105056817","EGU105056848","EGU105056881","EGU105056943","EGU105056960","EGU105056998","EGU105057013","EGU105057015","EGU105057019","EGU105057074","EGU105057090","EGU105057110","EGU105057130","EGU105057194","EGU105057235","EGU105057256","EGU105057335","EGU105057357","EGU105057553","EGU105057579","EGU105057634","EGU105057666","EGU105057669","EGU105057721","EGU105057742","EGU105057795","EGU105057841","EGU105057912","EGU105057919","EGU105057941","EGU105058078","EGU105058081","EGU105058083","EGU105058094","EGU105058107","EGU105058131","EGU105058170","EGU105058175","EGU105058180","EGU105058202","EGU105058237","EGU105058241","EGU105058245","EGU105058326","EGU105058366","EGU105058377","EGU105058418","EGU105058469","EGU105058499","EGU105058547","EGU105058556","EGU105058567","EGU105058576","EGU105058582","EGU105058592","EGU105058598","EGU105058614","EGU105058633","EGU105058682","EGU105058683","EGU105058687","EGU105058702","EGU105058723","EGU105058731","EGU105058781","EGU105058798","EGU105058802","EGU105058808","EGU105058863","EGU105058889","EGU105058890","EGU105058894","EGU105058904","EGU105058918","EGU105058989","EGU105058990","EGU105059003","EGU105059008","EGU105059023","EGU105059035","EGU105059042","EGU105059054","EGU105059108","EGU105059112","EGU105059113","EGU105059126","EGU105059131","EGU105059138","EGU105059176","EGU105059186","EGU105059193","EGU105059276","EGU105059342","EGU105059366","EGU105059367","EGU105059381","EGU105059441","EGU105059453","EGU105059458","EGU105059479","EGU105059480","EGU105059490","EGU105059570","EGU105059573","EGU105059575","EGU105059587","EGU105059612","EGU105059624","EGU105059636","EGU105059639","EGU105059671","EGU105059853","EGU105059900","EGU105060095","EGU105060589","EGU105061025","EGU105061385","EGU105061427","HEY1007","HEY1013","HEY1017","HEY1020","HEY1025","HEY1035","HEY1050","HEY1052","HEY1064","HEY110","HEY1119","HEY1168","HEY1171","HEY1197","HEY1201","HEY120","HEY122","HEY125","HEY12","HEY136","HEY139","HEY1484","HEY148","HEY1494","HEY14","HEY150","HEY1615","HEY164","HEY168","HEY17","HEY1801","HEY180","HEY1815","HEY182","HEY1842","HEY1854","HEY1877","HEY1901","HEY191","HEY194","HEY197","HEY1986","HEY201","HEY204e","HEY204s","HEY2056","HEY207","HEY215","HEY2164","HEY218","HEY21","HEY2238","HEY225","HEY226","HEY2291","HEY231","HEY2339","HEY2363","HEY2370","HEY2377","HEY237","HEY2388","HEY240","HEY2459","HEY245","HEY24","HEY250","HEY252e","HEY252p","HEY252s","HEY2550","HEY2561","HEY257","HEY267","HEY269","HEY277","HEY280","HEY281","HEY282","HEY290","HEY293","HEY296","HEY299","HEY305","HEY308","HEY310","HEY31","HEY323","HEY326","HEY32e","HEY32s","HEY332","HEY334","HEY340","HEY357","HEY360","HEY362","HEY363","HEY369","HEY378e","HEY378s","HEY38","HEY391","HEY392","HEY415","HEY417","HEY421","HEY429","HEY449","HEY464","HEY484","HEY490","HEY497","HEY4","HEY508","HEY514","HEY51","HEY52","HEY556","HEY563","HEY576","HEY581","HEY587","HEY604","HEY609","HEY61","HEY629","HEY630","HEY637","HEY673","HEY680","HEY703","HEY717","HEY727","HEY728","HEY732","HEY736","HEY740","HEY743","HEY757","HEY758","HEY762","HEY763","HEY785","HEY790","HEY793","HEY7","HEY807","HEY808","HEY822","HEY825","HEY82","HEY83","HEY84","HEY855","HEY856","HEY863","HEY872","HEY874","HEY883e","HEY883n","HEY886","HEY88","HEY897","HEY89","HEY938","HEY948","HEY94","HEY950","HEY958","HEY964","HEY977","HEY982","HEY985","HEY989"]


bad_genes = ["EGU105059594","EGU105059996"]


#Main workflow for trees
for i in range(len(genes)):

    ### Creating the partition files for each gene
    gwf.target_from_template('Partition_'+genes[i], partitioner(gene = genes[i],
                                                        path_in = "/home/owrisberg/Coryphoideae/work_flow/09_mapping/",
                                                        path_out = "/home/owrisberg/Coryphoideae/work_flow/10_tree_building/01_genetrees/partitions_and_clean_fastas/",
                                                        done = "/home/owrisberg/Coryphoideae/work_flow/10_tree_building/01_genetrees/done/partitioner/"+genes[i]))

	#Running IQ_tree
    gwf.target_from_template('IQtree_'+genes[i], iq_tree(gene = genes[i],
                                                        path_in = "/home/owrisberg/Coryphoideae/work_flow/10_tree_building/01_genetrees/partitions_and_clean_fastas/",
                                                        path_out = "/home/owrisberg/Coryphoideae/work_flow/10_tree_building/01_genetrees/"))
	
	#Running Rename Reroot
    gwf.target_from_template('RR_'+genes[i], rename_reroot(gene = genes[i],
                                                        path_in = "/home/owrisberg/Coryphoideae/work_flow/10_tree_building/01_genetrees/"))

# Gathering Genetrees into single file and contracting low support branches
gwf.target_from_template('Newick_Contracting_all', newick_contracting(path_in = "/home/owrisberg/Coryphoideae/work_flow/10_tree_building/01_genetrees/",
                                                        path_out = "/home/owrisberg/Coryphoideae/work_flow/10_tree_building/02_speciestree/"))

# Running Astral on the Genetrees
gwf.target_from_template('Astral', astral(path_in = "/home/owrisberg/Coryphoideae/work_flow/10_tree_building/02_speciestree/",
                                                        gene_tree_file="genetrees.tre",
														output="astral_tree.tre"))

# Renaming the tips
gwf.target_from_template('Renaming', renaming(path_in = "/home/owrisberg/Coryphoideae/work_flow/10_tree_building/02_speciestree/",
                                                        tree_in="astral_tree.tre",
														gene_tree_file="genetrees.tre",
														tree_out="astral_tree_renamed.tre"))

# Running Quartet scores
#gwf.target_from_template('Quartet_Scores', quartet_scores(path_in = "/home/owrisberg/Coryphoideae/work_flow/10_tree_building/02_speciestree/"))

# Running Astral_annotation on the Genetrees
gwf.target_from_template('Astral_annotation', astral_annotation(path_in = "/home/owrisberg/Coryphoideae/work_flow/10_tree_building/02_speciestree/",
                                                        gene_tree_file="genetrees.tre",
														species_tree_file="astral_tree.tre",
														outfile="astral_tree_annotated.tre"))

# Running SortaDate on the Astral tree using the genetrees
gwf.target_from_template('Sorta_date', sorta_date(path_in = "/home/owrisberg/Coryphoideae/work_flow/10_tree_building/01_genetrees/",
                                                        path_out ="/home/owrisberg/Coryphoideae/work_flow/11_dating_the_tree/00_sortadate/all_genes/",
														astral_tree="/home/owrisberg/Coryphoideae/work_flow/10_tree_building/02_speciestree/astral_tree.tre",
														done = "/home/owrisberg/Coryphoideae/work_flow/11_dating_the_tree/00_sortadate/all_genes/done_sorta_date"))



############################################################################################
###############################------Orthologs------########################################
############################################################################################

#Running simmilar analysis on only the genes which Sidonie deemed to be orthologs.
gwf.target_from_template('Newick_Contracting_orthologs', newick_contracting_orthologs(path_in = "/home/owrisberg/Coryphoideae/work_flow/10_tree_building/01_genetrees/",
                                                        path_out = "/home/owrisberg/Coryphoideae/work_flow/10_tree_building/02_speciestree/"))


# Running Astral on the Genetrees
gwf.target_from_template('Astral_orthologs', astral(path_in = "/home/owrisberg/Coryphoideae/work_flow/10_tree_building/02_speciestree/",
                                                        gene_tree_file="genetrees_orthologs.tre",
														output="astral_tree_orthologs.tre"))

# Renaming the tips
gwf.target_from_template('Renaming_orthologs', renaming(path_in = "/home/owrisberg/Coryphoideae/work_flow/10_tree_building/02_speciestree/",
                                                        tree_in="astral_tree_orthologs.tre",
														gene_tree_file="genetrees_orthologs.tre",
														tree_out="astral_tree_orthologs_renamed.tre"))

# Running Quartet scores
#gwf.target_from_template('Quartet_Scores', quartet_scores(path_in = "/home/owrisberg/Coryphoideae/work_flow/10_tree_building/02_speciestree/"))

# Running Astral_annotation on the Genetrees
gwf.target_from_template('Astral_annotation_orthologs', astral_annotation(path_in = "/home/owrisberg/Coryphoideae/work_flow/10_tree_building/02_speciestree/",
                                                        gene_tree_file="genetrees_orthologs.tre",
														species_tree_file="astral_tree_orthologs.tre",
														outfile="astral_tree_orthologs_annotated.tre"))

# Running SortaDate on the Astral tree using the genetrees
gwf.target_from_template('Sorta_date_orthologs', sorta_date(path_in = "/home/owrisberg/Coryphoideae/work_flow/10_tree_building/01_genetrees/copy_of_ortholog_gene_trees/",
                                                        path_out ="/home/owrisberg/Coryphoideae/work_flow/11_dating_the_tree/00_sortadate/orthologs/",
														astral_tree="/home/owrisberg/Coryphoideae/work_flow/10_tree_building/02_speciestree/astral_tree_orthologs.tre",
														done = "/home/owrisberg/Coryphoideae/work_flow/11_dating_the_tree/00_sortadate/orthologs/done_sorta_date_orthologs"))