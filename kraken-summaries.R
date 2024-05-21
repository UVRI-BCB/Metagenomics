rm(list = ls())

# load packages
require(ggplot2)
require(reshape)
require(dplyr)
require(pheatmap)

krakenDNAprep <- read.delim("kraken-combined-clean.txt",row.names = 1)
krakenDNAprep$Undetermined_S0.report<-NULL
krakenDNAprep$PC_S18.report<-NULL

ktmp <- krakenDNAprep

df1 <- data.frame(abundance=rowSums(krakenDNAprep))
df1$Taxa <- rownames(df1)
df1 <- df1[df1$abundance>50,]
krakenDNAprep$Taxa <- rownames(krakenDNAprep)
krakenDNAprep <- krakenDNAprep[krakenDNAprep$Taxa%in%df1$Taxa,]

pathogenic <- read.csv("genera.txt",header = F)
head(pathogenic)  

kraken.combined.clean <- krakenDNAprep
gs<-sub(".*g__","",kraken.combined.clean$Taxa)
gs
kraken.combined.clean$genus <- sub("\\|.*", "", gs)
kraken.combined.clean <- kraken.combined.clean[kraken.combined.clean$genus%in%pathogenic$V1,]
kraken.combined.clean$Species <- sub(".*s__","",kraken.combined.clean$Taxa)
rownames(kraken.combined.clean)<-kraken.combined.clean$Species

k_bact <- kraken.combined.clean

kraken.combined.clean$Taxa<-NULL
kraken.combined.clean$genus<-NULL
kraken.combined.clean$Species<-NULL

df3 <- data.frame(counts=rowSums(kraken.combined.clean))
df3$species <- rownames(df3)
df3<-df3[df3$counts>10000,]
df3<-df3[order(df3$counts, decreasing = TRUE),]
df3$species<-factor(df3$species,levels=df3$species, labels = df3$species)
df3$Phyla<-phylum_krak$phylum[match(
  df3$species, phylum_krak$Species
)]

p0<-ggplot(df3, aes(x=species, y = log2(counts), fill=Phyla)) + #fill="#4682B4"
  geom_bar(stat = "identity")+
  theme(axis.title.x = element_text(color="black", size=15))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15))+
  theme(axis.text.y = element_text(color="black", size=12))+
  theme(axis.text.x = element_text(angle=0)
        ,panel.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(size=12),#
        panel.border = element_rect(fill=NA),
        title = element_text(size = 12)#, face = "bold")
  )+
  scale_y_continuous(labels = function(x) format(x, digits=1))+
  scale_fill_brewer(palette = "Dark2", direction = 1)+
  coord_flip()+theme(legend.position = "right")
p0$labels$fill<-"Phyla"
p0$labels$x<-'Bacteria'
p0$labels$y<-'Log2 transformed read counts'
p0
ggsave("~/Desktop/species-phyla-counts.pdf", plot=p0, device = "pdf", width = 10, height = 12)


kraken.combined.clean <- kraken.combined.clean[rownames(kraken.combined.clean)%in%df3$species,]

kraken.combined.clean<-kraken.combined.clean[sort(colnames(kraken.combined.clean))]

colnames(kraken.combined.clean) <- sub("_.*","",colnames(kraken.combined.clean))
colnames(kraken.combined.clean)<-paste0("ENV-", 1:dim(kraken.combined.clean)[2])
pj<-pheatmap(as.matrix(log2(kraken.combined.clean+1)), 
             cluster_cols = TRUE,
             cluster_rows = TRUE, 
             show_colnames = TRUE,
             treeheight_row = 0, 
             treeheight_col = 0
)
ggsave("~/Desktop/WW-Bacteria-10000-reads-and-above.pdf", plot=pj, device = "pdf", width = 10, height = 15)

genus_krak <- k_bact
genus_krak$Taxa<-NULL
genus_krak$Species<-NULL
genus_krak$Undetermined_S0.report<-NULL

# Aggregate by 'group' column while keeping all other columns
genus_krak <- aggregate(. ~ genus, data = genus_krak, FUN = sum)
rownames(genus_krak) <- genus_krak$genus
genus_krak$genus<-NULL
genus_df <- data.frame(counts=rowSums(genus_krak))
genus_df$genus <- rownames(genus_df)
genus_df<-genus_df[order(genus_df$counts, decreasing = FALSE),]
genus_df$genus<-factor(genus_df$genus,levels=genus_df$genus, labels = genus_df$genus)
genus_df$phylum <- phylum_krak$phylum[match(genus_df$genus, phylum_krak$genus)]

p0<-ggplot(genus_df, aes(x=genus, y = log2(counts), fill=phylum)) + #fill="#4682B4"
  geom_bar(stat = "identity")+
  theme(axis.title.x = element_text(color="black", size=15))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15))+
  theme(axis.text.y = element_text(color="black", size=12))+
  theme(axis.text.x = element_text(angle=0)
        ,panel.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(size=12),#
        panel.border = element_rect(fill=NA),
        title = element_text(size = 12)#, face = "bold")
  )+
  scale_y_continuous(labels = function(x) format(x, digits=1))+
  scale_fill_brewer(palette = "Dark2", direction = 1)+
  coord_flip()+theme(legend.position = "right")
p0$labels$fill<-"Phyla"
p0$labels$x<-'Genera'
p0$labels$y<-'Log2 transformed read counts'
p0
ggsave("~/Desktop/genera-phyla-counts.pdf", plot=p0, device = "pdf", width = 10, height = 10)

colnames(genus_krak)<-paste0("ENV-", 1:dim(genus_krak)[2])
pj<-pheatmap(as.matrix(log2(genus_krak+1)), 
             cluster_cols = TRUE,
             cluster_rows = TRUE, 
             show_colnames = TRUE,
             treeheight_row = 0, 
             treeheight_col = 0
)
ggsave("~/Desktop/WW-Bacteria-genus.pdf", plot=pj, device = "pdf", width = 10, height = 10)

#family
family_genus_krak <- k_bact
family_genus_krak$family <- stringr::str_extract(family_genus_krak$Taxa, '(?<=f__)[^|]+')

# Aggregate by 'group' column while keeping all other columns
family_genus_krak$Species<-NULL
family_genus_krak$genus<-NULL
family_genus_krak$Taxa<-NULL

family_genus_krak <- aggregate(. ~ family, data = family_genus_krak, FUN = sum)
rownames(family_genus_krak) <- family_genus_krak$family
family_genus_krak$family<-NULL

colnames(family_genus_krak)<-paste0("ENV-", 1:dim(family_genus_krak)[2])
pj<-pheatmap(as.matrix(log2(family_genus_krak+1)), 
             cluster_cols = TRUE,
             cluster_rows = TRUE, 
             show_colnames = TRUE,
             treeheight_row = 0, 
             treeheight_col = 0
)
ggsave("~/Desktop/WW-Bacteria-family.pdf", plot=pj, device = "pdf", width = 10, height = 10)

#phylum
phylum_krak <- k_bact
phylum_krak$phylum <- stringr::str_extract(phylum_krak$Taxa, '(?<=p__)[^|]+')

# Aggregate by 'group' column while keeping all other columns
phylum_krak$Species<-NULL
phylum_krak$genus<-NULL
phylum_krak$Taxa<-NULL

phylum_krak <- aggregate(. ~ phylum, data = phylum_krak, FUN = sum)
rownames(phylum_krak) <- phylum_krak$phylum
phylum_krak$phylum<-NULL

colnames(phylum_krak)<-paste0("ENV-", 1:dim(phylum_krak)[2])
pj<-pheatmap(as.matrix(log2(phylum_krak+1)), 
             cluster_cols = TRUE,
             cluster_rows = TRUE, 
             show_colnames = TRUE,
             treeheight_row = 0, 
             treeheight_col = 0
)
ggsave("~/Desktop/WW-Bacteria-phyla.pdf", plot=pj, device = "pdf", width = 10, height = 10)


# viruses
kraken.combined.clean<-ktmp
#kraken.combined.clean$Undetermined_S0.report<-NULL
colnames(kraken.combined.clean)<-sub(".report","",colnames(kraken.combined.clean))

###########
colnames(kraken.combined.clean)<-gsub("_",".",colnames(kraken.combined.clean))
no.of.viralhits<-data.frame(colSums(kraken.combined.clean))
no.of.viralhits$sample<-rownames(no.of.viralhits)
#write.csv(no.of.viralhits, "number-of-viralhits.csv")

############
kraken.combined.clean$Tax<-rownames(kraken.combined.clean)
kraken.combined.clean<-kraken.combined.clean[grepl('k__Viruses',kraken.combined.clean$Tax),]
kraken.combined.clean<-kraken.combined.clean[grepl('f__',kraken.combined.clean$Tax),]
kraken.combined.clean$Family<-sub("\\|g__.*","",sub("\\|s__.*","",sub(".*f__","",kraken.combined.clean$Tax)))
kraken.combined.clean$Species <- sub(".*s__","",kraken.combined.clean$Tax)
rownames(kraken.combined.clean)<-kraken.combined.clean$Species
########## Start summarizing by virus family ##########

k<-kraken.combined.clean

virus_fam <- k
virus_fam$Species<-NULL  
virus_fam$Tax<-NULL

virus_fam_krak <- aggregate(. ~ Family, data = virus_fam, FUN = sum)
rownames(virus_fam_krak) <- virus_fam_krak$Family
virus_fam_krak$Family<-NULL

View(data.frame(rowSums(virus_fam_krak)))
############# Heatmap #################################
kraken.combined.clean$species2<-sub("_.*","",kraken.combined.clean$Species)
#kraken.combined.clean<-kraken.combined.clean[kraken.combined.clean$species2!="Tequintavirus",]
#k2<-kraken.combined.clean["A22.10.002_S10"]
#k2$Species<-rownames(k2)

kraken.combined.clean$Virus<-NULL
kraken.combined.clean$Tax<-NULL
kraken.combined.clean$species2<-NULL
kraken.combined.clean$Family<-NULL
tmp_df <- reshape2::melt(kraken.combined.clean, id.vars = "Species")
#write.csv(tmp_df, file = "classified-05-reads-and-above.csv", quote = F, row.names = F)
kraken.combined.clean$Species<-NULL

virus_counts <- rowSums(kraken.combined.clean)
virus_counts <- data.frame(reads=virus_counts, virus=names(virus_counts))
#virus_counts<-virus_counts[virus_counts$reads>=10, ]
virus_counts<-virus_counts[order(virus_counts$reads, decreasing = T),]
kraken.combined.clean<-kraken.combined.clean[rownames(kraken.combined.clean)%in%virus_counts$virus,]
rownames(kraken.combined.clean)<-gsub("_"," ",rownames(kraken.combined.clean))

tf<-rowSums(kraken.combined.clean)
tf<-data.frame(reads=tf,species=names(tf))
tf<-tf[order(tf$reads, decreasing = F),]
tf$species<-factor(tf$species, levels = tf$species, labels = tf$species)

top50<-tf[tf$reads>=20,]

p0<-ggplot(top50, aes(x=species, y = log2(reads))) + #fill="#4682B4"
  geom_bar(stat = "identity", fill="#6082B6")+ylab("Log2 (number of reads)")+xlab("")+
  #geom_text(aes(label=paste0(reads+100)))+
  theme(axis.title.x = element_text(color="black", size=15,face = "bold"))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15,face = "bold"))+
  theme(axis.text.y = element_text(color="black", size=12))+
  #scale_fill_manual("Mutational frequency",values=rev(c("#DAA520","#4682B4"))) +
  theme(axis.text.x = element_text(angle=0)
        ,panel.background = element_blank(),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(size=15),#
        panel.border = element_rect(fill=NA),
        title = element_text(size = 12, face = "bold")
  )+
  scale_y_continuous(labels = function(x) format(x, digits=1))+
  scale_fill_brewer(palette = "Dark2", direction = 1)+
  coord_flip()+theme(legend.position = "none")
p0
ggsave("~/Desktop/viruses-counts.pdf", plot=p0, device = "pdf", width = 10, height = 12)

colnames(kraken.combined.clean)<-paste0("ENV-", 1:dim(kraken.combined.clean)[2])
k2<-kraken.combined.clean[rownames(kraken.combined.clean)%in%top50$species,]

pj<-pheatmap(as.matrix(log2(k2+1)), 
             cluster_cols = TRUE,
             cluster_rows = TRUE, 
             show_colnames = TRUE,
             treeheight_row = 0, 
             treeheight_col = 0
)
#ggsave("Virus-hits-10-reads-and-above.pdf", plot=pj, device = "pdf", width = 10, height = 10)
ggsave("~/Desktop/WW-Viruses-10-reads-and-above.pdf", plot=pj, device = "pdf", width = 10, height = 12)


# ARCHEA

