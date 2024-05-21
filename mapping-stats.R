rm(list = ls())

require(ggplot2)

samstats <- read.csv("mapping-stats.csv")

df <- samstats
df <- df[df$sample!="Undetermined_S0:",]
df$sample <- paste0("ARUG-BA-",1:length(df$sample))
df <- df[df$coverage>0.4,]
df<-df[order(df$coverage, decreasing = FALSE),]
df$sample<-factor(df$sample, levels = df$sample, labels = df$sample)

p<- ggplot(data = df, aes(x=sample,y=coverage))+ylab("Percentage genome coverage")+xlab("Samples")+
  geom_bar(stat = "identity", fill="#6082B6")+
  #geom_text(aes(label=paste0(reads+100)))+
  theme(axis.title.x = element_text(color="black", size=15))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15))+
  theme(axis.text.y = element_text(color="black", size=12))+
  #scale_fill_manual("Mutational frequency",values=rev(c("#DAA520","#4682B4"))) +
  theme(axis.text.x = element_text(angle=0)
        ,panel.background = element_blank(),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(size=15),#
        panel.border = element_rect(fill=NA),
        title = element_text(size = 12)
  )+
  scale_y_continuous(labels = function(x) format(x, digits=1))+
  scale_fill_brewer(palette = "Dark2", direction = 1)+
  coord_flip()+theme(legend.position = "none")
p

#pick samples with top coverages
top_anthrax_samples <- samstats$sample[samstats$coverage>10]

depths <- read.delim("~/Desktop/arua-tmp/depths.txt", header=FALSE)
depths$V4 <- gsub('-anth', '', depths$V4)
depths<-depths[depths$V3>=10,]
depths2 <- depths[depths$V4%in%top_anthrax_samples, ]

depths1 <- depths2[depths2$V4==top_anthrax_samples[1],]
depths3 <- depths2[depths2$V4==top_anthrax_samples[3],]
depths4 <- depths2[depths2$V4==top_anthrax_samples[4],]

p4<-ggplot2::ggplot(data = depths1, aes(x=V2, y=log2(V3)))+ 
  #facet_wrap(~V4, nrow = 4, scales = 'free_y')+
  theme(axis.title.x = element_text(color="black", size=15))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15))+
  theme(axis.text.y = element_text(color="black", size=12))+
  #scale_fill_manual("Mutational frequency",values=rev(c("#DAA520","#4682B4"))) +
  theme(axis.text.x = element_text(angle=0)
        ,panel.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(size=12),#
        panel.border = element_rect(fill=NA),
        title = element_text(size = 12)
  )+#scale_x_continuous(labels = scales::scientific_format())+
  geom_col(color='#3C5488B2',width = 0.2)+ylab('Log2 (Depth of coverage)')+xlab('Position along the reference genome')
#p4
ggsave('~/Desktop/A24-03-014-041_S4-Bacillus-Coverage-plot.pdf', plot = p4, device = 'pdf', width = 12, height = 4)


p5<-ggplot2::ggplot(data = depths3, aes(x=V2, y=log2(V3)))+ 
  #facet_wrap(~V4, nrow = 4, scales = 'free_y')+
  theme(axis.title.x = element_text(color="black", size=15))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15))+
  theme(axis.text.y = element_text(color="black", size=12))+
  #scale_fill_manual("Mutational frequency",values=rev(c("#DAA520","#4682B4"))) +
  theme(axis.text.x = element_text(angle=0)
        ,panel.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(size=12),#
        panel.border = element_rect(fill=NA),
        title = element_text(size = 12)
  )+#scale_x_continuous(labels = scales::scientific_format())+
  geom_col(color='#3C5488B2',width = 0.2)+ylab('Log2 (Depth of coverage)')+xlab('Position along the reference genome')
#p5
ggsave('~/Desktop/A24-03-020-057_S6-Bacillus-Coverage-plot.pdf', plot = p5, device = 'pdf', width = 12, height = 4)

p6<-ggplot2::ggplot(data = depths4, aes(x=V2, y=log2(V3)))+ 
  #facet_wrap(~V4, nrow = 4, scales = 'free_y')+
  theme(axis.title.x = element_text(color="black", size=15))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15))+
  theme(axis.text.y = element_text(color="black", size=12))+
  #scale_fill_manual("Mutational frequency",values=rev(c("#DAA520","#4682B4"))) +
  theme(axis.text.x = element_text(angle=0)
        ,panel.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(size=12),#
        panel.border = element_rect(fill=NA),
        title = element_text(size = 12)
  )+#scale_x_continuous(labels = scales::scientific_format())+
  geom_col(color='#3C5488B2',width = 0.2)+ylab('Log2 (Depth of coverage)')+xlab('Position along the reference genome')
#p6
ggsave('~/Desktop/A24-03-007-027_S7-Bacillus-Coverage-plot.pdf', plot = p6, device = 'pdf', width = 12, height = 4)

arua3.depths <- read.delim("~/Desktop/arua-tmp/arua3-depths.txt", header=FALSE)
arua3<-arua3.depths[arua3.depths$V1=='A24-03-066-215-Blood_S19',]
arua3<-arua3[arua3$V4>=10,]
p7<-ggplot2::ggplot(data = arua3, aes(x=V3, y=log2(V4)))+ 
  #facet_wrap(~V4, nrow = 4, scales = 'free_y')+
  theme(axis.title.x = element_text(color="black", size=15))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15))+
  theme(axis.text.y = element_text(color="black", size=12))+
  #scale_fill_manual("Mutational frequency",values=rev(c("#DAA520","#4682B4"))) +
  theme(axis.text.x = element_text(angle=0)
        ,panel.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(size=12),#
        panel.border = element_rect(fill=NA),
        title = element_text(size = 12)
  )+#scale_x_continuous(labels = scales::scientific_format())+
  geom_col(color='#3C5488B2',width = 0.2)+ylab('Log2 (Depth of coverage)')+xlab('Position along the reference genome')
#p7
ggsave('~/Desktop/A24-03-066-215-Blood_S19-Bacillus-Coverage-plot.pdf', plot = p7, device = 'pdf', width = 12, height = 4)

p20 <- cowplot::plot_grid(p4, p5, p6, p7, nrow=2)
p21 <- cowplot::plot_grid(p, p20, ncol = 2)
ggsave("~/Desktop/coverage-plots-2.pdf", p21, device = 'pdf', width = 20, height = 10)

blast <- read.delim("~/Desktop/arua-tmp/final-blast.txt", header=FALSE)
dim(blast)
blast<-blast[blast$V4>=1000,]
dim(blast)
blast <- blast[blast$V3==100,]
dim(blast)

blast <- blast[blast$V7==1,]
dim(blast)

length(unique(blast$V13))
length(unique(blast$V2))

blastn_short <- blast[c('V13','V1','V2')]
dim(blastn_short)
blastn_short <- blastn_short[!duplicated(blastn_short),]
dim(blastn_short)

annot <- read.delim("~/Desktop/arua-tmp/anthrax-strains.txt", header=FALSE)

blastn_short$Strain <- annot$V2[match(blastn_short$V2, trimws(annot$V1))]

df4 <- data.frame(table(blastn_short$Strain))
unique(df4$Var1)
df4 <- df4[df4$Freq>20,]
unique(df4$Var2)


rm(list = ls())
df <- read.csv("~/Desktop/arua-tmp/mapping-stats.csv")
df <- df[df$sample!="Undetermined_S0:",]
df$sample <- paste0("ARUG-BA-",1:length(df$sample))
df <- df[df$coverage>0.4,]
df<-df[order(df$coverage, decreasing = FALSE),]
df$sample<-factor(df$sample, levels = df$sample, labels = df$sample)

p<- ggplot(data = df, aes(x=sample,y=coverage))+ylab("Percentage genome coverage")+xlab("Samples")+
  geom_bar(stat = "identity", fill="#6082B6")+
  #geom_text(aes(label=paste0(reads+100)))+
  theme(axis.title.x = element_text(color="black", size=15))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15))+
  theme(axis.text.y = element_text(color="black", size=12))+
  #scale_fill_manual("Mutational frequency",values=rev(c("#DAA520","#4682B4"))) +
  theme(axis.text.x = element_text(angle=0)
        ,panel.background = element_blank(),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(size=15),#
        panel.border = element_rect(fill=NA),
        title = element_text(size = 12)
  )+
  scale_y_continuous(labels = function(x) format(x, digits=1))+
  scale_fill_brewer(palette = "Dark2", direction = 1)+
  coord_flip()+theme(legend.position = "none")
p

df <- df[df$coverage>0.4,]
df<-df[order(df$coverage, decreasing = FALSE),]
df$sample<-factor(df$sample, levels = df$sample, labels = df$sample)

p<- ggplot(data = df, aes(x=sample,y=avgdepth))+ylab("Average depth")+xlab("Samples")+
  geom_bar(stat = "identity", fill="#6082B6")+
  #geom_text(aes(label=paste0(reads+100)))+
  theme(axis.title.x = element_text(color="black", size=15))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15))+
  theme(axis.text.y = element_text(color="black", size=12))+
  #scale_fill_manual("Mutational frequency",values=rev(c("#DAA520","#4682B4"))) +
  theme(axis.text.x = element_text(angle=0)
        ,panel.background = element_blank(),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(size=15),#
        panel.border = element_rect(fill=NA),
        title = element_text(size = 12)
  )+
  scale_y_continuous(labels = function(x) format(x, digits=1))+
  scale_fill_brewer(palette = "Dark2", direction = 1)+
  coord_flip()+theme(legend.position = "none")
p

