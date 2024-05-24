library(ape)
library(readr)

#data_dir <- "/Users/au265104/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/gis07.st.client.au.dk – SFTP/data_vol/wolf/Dypsis"
data_dir <- "/home/au543206/Documents/Coryphoideae/Figures_in_r/data"

# read name translation table (SECAPR No. to tip name)
rename <- read.table("/home/au543206/GenomeDK/Coryphoideae/github_code/coryphoideae_species_tree/names_for_tips.csv", sep=",", colClasses = "character")

# read figurename translation table (SECAPR No. to figure name)
figurename <- read.table("/home/au543206/GenomeDK/Coryphoideae/github_code/coryphoideae_species_tree/names_for_tips.csv", sep=",", colClasses = "character")
#figurename <- figurename[1:518,]
figurename_idx <- figurename$V2
names(figurename_idx) <- figurename$V1
figurename_idx

# read species trees
astral_tree <- read.tree(paste(data_dir, "/astral_tree_renamed.tre", sep=""))
astral_tree_orthologs <- read.tree(paste(data_dir, "/astral_tree_orthologs_renamed.tre", sep=""))
astral_tree_orthologs <- drop.tip(astral_tree_orthologs, "Maxburretiagracilis")


#Removing troublesome labels that we dont need
#Tree using all genes
tree_string <- read_file(paste(data_dir,"/astral_tree_annotated.tre", sep = ""))
tree_string <- gsub('\\[q1.*?QC=[0-9]*;EN=', '', tree_string)
tree_string <- gsub("\\]", "", tree_string)
tree_string <- gsub("'","", tree_string)
write_file(tree_string,paste0(data_dir,"/astral_tree_annotated_EN.tre"), append = FALSE)
astral_tree_EN <- ape::read.tree(paste(data_dir, "/astral_tree_annotated_EN.tre", sep=""))


#Tree using only orthologues genes
tree_string_orthologs <- read_file(paste(data_dir,"/astral_tree_orthologs_annotated.tre", sep = ""))
tree_string_orthologs <- gsub('\\[q1.*?QC=[0-9]*;EN=', '', tree_string_orthologs)
tree_string_orthologs <- gsub("\\]", "", tree_string_orthologs)
tree_string_orthologs <- gsub("'","", tree_string_orthologs)
write_file(tree_string_orthologs,paste0(data_dir,"/astral_tree_orthologs_annotated_EN.tre"), append = FALSE)
astral_tree_orthologs_EN <- ape::read.tree(paste(data_dir, "/astral_tree_orthologs_annotated_EN.tre", sep=""))
astral_tree_orthologs_EN <- drop.tip(astral_tree_orthologs_EN, "4039")

# Create a list of all the genetrees
gtnames = list.files(paste("/home/au543206/GenomeDK/Coryphoideae/work_flow/10_tree_building/01_genetrees", sep=""), pattern="*_rooted.tre", full.names=FALSE)


gts <- list()
for(i in 1:length(gtnames)){
  print(paste("/home/au543206/GenomeDK/Coryphoideae/work_flow/10_tree_building/01_genetrees/", gtnames[i], sep=""))
  gts[[i]] <- read.tree(paste("/home/au543206/GenomeDK/Coryphoideae/work_flow/10_tree_building/01_genetrees/", gtnames[i], sep=""))
}
rm(i)

# Create a list of all the ortholog genetrees
gtnames_orthologs <- list(
  "EGU105035046_rooted.tre",
  "EGU105035203_rooted.tre",
  "EGU105035555_rooted.tre",
  "EGU105035989_rooted.tre",
  "EGU105036774_rooted.tre",
  "EGU105038201_rooted.tre",
  "EGU105038431_rooted.tre",
  "EGU105038513_rooted.tre",
  "EGU105038747_rooted.tre",
  "EGU105038832_rooted.tre",
  "EGU105039082_rooted.tre",
  "EGU105039099_rooted.tre",
  "EGU105039164_rooted.tre",
  "EGU105039255_rooted.tre",
  "EGU105039449_rooted.tre",
  "EGU105039494_rooted.tre",
  "EGU105039501_rooted.tre",
  "EGU105039512_rooted.tre",
  "EGU105039542_rooted.tre",
  "EGU105039660_rooted.tre",
  "EGU105039690_rooted.tre",
  "EGU105039783_rooted.tre",
  "EGU105040185_rooted.tre",
  "EGU105040242_rooted.tre",
  "EGU105040462_rooted.tre",
  "EGU105040842_rooted.tre",
  "EGU105040851_rooted.tre",
  "EGU105040914_rooted.tre",
  "EGU105040918_rooted.tre",
  "EGU105040970_rooted.tre",
  "EGU105041100_rooted.tre",
  "EGU105041179_rooted.tre",
  "EGU105041217_rooted.tre",
  "EGU105041337_rooted.tre",
  "EGU105041650_rooted.tre",
  "EGU105041903_rooted.tre",
  "EGU105041933_rooted.tre",
  "EGU105042090_rooted.tre",
  "EGU105042633_rooted.tre",
  "EGU105042781_rooted.tre",
  "EGU105043037_rooted.tre",
  "EGU105043469_rooted.tre",
  "EGU105043485_rooted.tre",
  "EGU105043633_rooted.tre",
  "EGU105043686_rooted.tre",
  "EGU105043786_rooted.tre",
  "EGU105044407_rooted.tre",
  "EGU105044668_rooted.tre",
  "EGU105044710_rooted.tre",
  "EGU105045078_rooted.tre",
  "EGU105045099_rooted.tre",
  "EGU105045282_rooted.tre",
  "EGU105045367_rooted.tre",
  "EGU105045467_rooted.tre",
  "EGU105046245_rooted.tre",
  "EGU105046456_rooted.tre",
  "EGU105047379_rooted.tre",
  "EGU105047434_rooted.tre",
  "EGU105047446_rooted.tre",
  "EGU105047553_rooted.tre",
  "EGU105047621_rooted.tre",
  "EGU105047790_rooted.tre",
  "EGU105047940_rooted.tre",
  "EGU105048476_rooted.tre",
  "EGU105048541_rooted.tre",
  "EGU105048694_rooted.tre",
  "EGU105048909_rooted.tre",
  "EGU105048961_rooted.tre",
  "EGU105049052_rooted.tre",
  "EGU105049312_rooted.tre",
  "EGU105049360_rooted.tre",
  "EGU105049539_rooted.tre",
  "EGU105050126_rooted.tre",
  "EGU105050344_rooted.tre",
  "EGU105050366_rooted.tre",
  "EGU105050383_rooted.tre",
  "EGU105050521_rooted.tre",
  "EGU105050682_rooted.tre",
  "EGU105050853_rooted.tre",
  "EGU105051156_rooted.tre",
  "EGU105051188_rooted.tre",
  "EGU105051362_rooted.tre",
  "EGU105051403_rooted.tre",
  "EGU105051560_rooted.tre",
  "EGU105051726_rooted.tre",
  "EGU105051764_rooted.tre",
  "EGU105051795_rooted.tre",
  "EGU105051847_rooted.tre",
  "EGU105051860_rooted.tre",
  "EGU105051870_rooted.tre",
  "EGU105051985_rooted.tre",
  "EGU105052307_rooted.tre",
  "EGU105052476_rooted.tre",
  "EGU105052492_rooted.tre",
  "EGU105052580_rooted.tre",
  "EGU105052804_rooted.tre",
  "EGU105052818_rooted.tre",
  "EGU105052855_rooted.tre",
  "EGU105052888_rooted.tre",
  "EGU105052956_rooted.tre",
  "EGU105053006_rooted.tre",
  "EGU105053055_rooted.tre",
  "EGU105053079_rooted.tre",
  "EGU105053422_rooted.tre",
  "EGU105053482_rooted.tre",
  "EGU105053549_rooted.tre",
  "EGU105053848_rooted.tre",
  "EGU105053866_rooted.tre",
  "EGU105053889_rooted.tre",
  "EGU105053901_rooted.tre",
  "EGU105053980_rooted.tre",
  "EGU105054153_rooted.tre",
  "EGU105054405_rooted.tre",
  "EGU105054595_rooted.tre",
  "EGU105054786_rooted.tre",
  "EGU105054930_rooted.tre",
  "EGU105054948_rooted.tre",
  "EGU105055065_rooted.tre",
  "EGU105055072_rooted.tre",
  "EGU105055075_rooted.tre",
  "EGU105055114_rooted.tre",
  "EGU105055115_rooted.tre",
  "EGU105055499_rooted.tre",
  "EGU105055664_rooted.tre",
  "EGU105055800_rooted.tre",
  "EGU105055873_rooted.tre",
  "EGU105056091_rooted.tre",
  "EGU105056289_rooted.tre",
  "EGU105056469_rooted.tre",
  "EGU105056654_rooted.tre",
  "EGU105057074_rooted.tre",
  "EGU105057634_rooted.tre",
  "EGU105057666_rooted.tre",
  "EGU105057721_rooted.tre",
  "EGU105058081_rooted.tre",
  "EGU105058131_rooted.tre",
  "EGU105058170_rooted.tre",
  "EGU105058180_rooted.tre",
  "EGU105058237_rooted.tre",
  "EGU105058469_rooted.tre",
  "EGU105058702_rooted.tre",
  "EGU105058798_rooted.tre",
  "EGU105058889_rooted.tre",
  "EGU105059023_rooted.tre",
  "EGU105059126_rooted.tre",
  "EGU105059138_rooted.tre",
  "EGU105059186_rooted.tre",
  "EGU105059366_rooted.tre",
  "EGU105059381_rooted.tre",
  "EGU105059480_rooted.tre",
  "EGU105059624_rooted.tre",
  "EGU105059639_rooted.tre",
  "EGU105059853_rooted.tre",
  "EGU105061427_rooted.tre",
  "HEY1007_rooted.tre",
  "HEY1017_rooted.tre",
  "HEY1020_rooted.tre",
  "HEY1035_rooted.tre",
  "HEY1052_rooted.tre",
  "HEY1064_rooted.tre",
  "HEY1119_rooted.tre",
  "HEY1197_rooted.tre",
  "HEY1201_rooted.tre",
  "HEY120_rooted.tre",
  "HEY122_rooted.tre",
  "HEY12_rooted.tre",
  "HEY150_rooted.tre",
  "HEY1615_rooted.tre",
  "HEY1815_rooted.tre",
  "HEY182_rooted.tre",
  "HEY1854_rooted.tre",
  "HEY194_rooted.tre",
  "HEY197_rooted.tre",
  "HEY201_rooted.tre",
  "HEY204e_rooted.tre",
  "HEY2056_rooted.tre",
  "HEY2238_rooted.tre",
  "HEY226_rooted.tre",
  "HEY231_rooted.tre",
  "HEY250_rooted.tre",
  "HEY252e_rooted.tre",
  "HEY252p_rooted.tre",
  "HEY252s_rooted.tre",
  "HEY2550_rooted.tre",
  "HEY2561_rooted.tre",
  "HEY257_rooted.tre",
  "HEY269_rooted.tre",
  "HEY277_rooted.tre",
  "HEY281_rooted.tre",
  "HEY282_rooted.tre",
  "HEY290_rooted.tre",
  "HEY293_rooted.tre",
  "HEY299_rooted.tre",
  "HEY305_rooted.tre",
  "HEY31_rooted.tre",
  "HEY326_rooted.tre",
  "HEY32e_rooted.tre",
  "HEY32s_rooted.tre",
  "HEY340_rooted.tre",
  "HEY357_rooted.tre",
  "HEY362_rooted.tre",
  "HEY363_rooted.tre",
  "HEY369_rooted.tre",
  "HEY378s_rooted.tre",
  "HEY38_rooted.tre",
  "HEY392_rooted.tre",
  "HEY51_rooted.tre",
  "HEY576_rooted.tre",
  "HEY587_rooted.tre",
  "HEY604_rooted.tre",
  "HEY61_rooted.tre",
  "HEY629_rooted.tre",
  "HEY630_rooted.tre",
  "HEY637_rooted.tre",
  "HEY703_rooted.tre",
  "HEY740_rooted.tre",
  "HEY758_rooted.tre",
  "HEY790_rooted.tre",
  "HEY7_rooted.tre",
  "HEY807_rooted.tre",
  "HEY83_rooted.tre",
  "HEY855_rooted.tre",
  "HEY863_rooted.tre",
  "HEY872_rooted.tre",
  "HEY886_rooted.tre",
  "HEY88_rooted.tre",
  "HEY948_rooted.tre",
  "HEY94_rooted.tre",
  "HEY950_rooted.tre",
  "HEY977_rooted.tre",
  "HEY985_rooted.tre"
)

gts_orthologs <- list()
for(i in 1:length(gtnames_orthologs)){
  print(paste("/home/au543206/GenomeDK/Coryphoideae/work_flow/10_tree_building/01_genetrees/", gtnames_orthologs[i], sep=""))
  gts_orthologs[[i]] <- read.tree(paste("/home/au543206/GenomeDK/Coryphoideae/work_flow/10_tree_building/01_genetrees/", gtnames_orthologs[i], sep=""))
}
rm(i)


