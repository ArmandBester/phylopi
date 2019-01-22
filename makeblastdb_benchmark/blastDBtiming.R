
library(ggplot2)


system("scp pi@196.254.115.64:/home/pi/blastDataBaseTiming/makeBlastTimes.txt .")


timeFile <- read.csv("makeBlastTimes.txt", header = F)

colnames(timeFile) <- c("n", "s")


lin <- lm(s ~ n, data = timeFile)

summary(lin)

p <- ggplot(timeFile, aes(x = n, y = s))+
  geom_point(size = 0.5)+
  geom_smooth(method = 'lm', color = "black", linetype = "dashed")+
  theme_bw()+
  xlab("Number of sequences")+
  ylab("Time (s)")
p

ggsave(filename = "makeblastdb.tiff", plot = p, dpi = 300)
