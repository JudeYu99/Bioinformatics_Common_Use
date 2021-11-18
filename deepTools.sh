### deepTools Analysis 
# @ Prepare .bw file and .bed file(narrowPeak file in this example)
# @ Author: Yu Zhu
# @ Email: 1830416012@stu.suda.edu.cn
####################################################################################

computeMatrix reference-point -S SRR035186_trimed_sort.bw -R ../sxjns/macs2.narrowPeak.filt_blacklist --referencePoint center -a 2000 -b 2000 -out matrix_black_filt.tab.gz

plotHeatmap -m matrix_black_filt.tab.gz -out hm_black_filt.png --heatmapHeight 15 --refPointLabel TSS


computeMatrix reference-point -S SRR035187_trimed_sort.bw -R ../sxjns/macs2.narrowPeak.filt_blacklist --referencePoint center -a 2000 -b 2000 -out matrix_black_filt_2.tab.gz

plotHeatmap -m matrix_black_filt_2.tab.gz -out hm_black_filt_2.png --heatmapHeight 15 --refPointLabel TSS


computeMatrix reference-point -S SRR015350_sort.bw -R ../sxjns/macs2.narrowPeak.filt_blacklist --referencePoint center -a 2000 -b 2000 -out matrix_black_filt_1.tab.gz

plotHeatmap -m matrix_black_filt_1.tab.gz -out hm_matrix_black_filt_1.png --heatmapHeight 15 --refPointLabel TSS


computeMatrix reference-point -S SRR015349_sort.bw -R ../sxjns/macs2.narrowPeak.filt_blacklist --referencePoint center -a 2000 -b 2000 -out matrix_black_filt_2.tab.gz

plotHeatmap -m matrix_black_filt_2.tab.gz -out hm_matrix_black_filt_2.png --heatmapHeight 15 --refPointLabel TSS

# @ Cite: Ramírez, Fidel, Devon P. Ryan, Björn Grüning, Vivek Bhardwaj, Fabian Kilpert, Andreas S. Richter, Steffen Heyne, Friederike Dündar, and Thomas Manke. deepTools2: A next Generation Web Server for Deep-Sequencing Data Analysis. Nucleic Acids Research (2016).