# **Bam_Spatial_Velocity**

bam ---> raw gem --> gem --> Annadata --> ScanPy&Scvelo

:smiley:

## 1.需要准备的文件

1. gtf文件：不多解释，放了demo
2. bam文件：需要有GE:Z ; XF:Z ; CB:Z 标签（放了demo，可以查看STAR比对的参数）
3. poison range文件：一列多行，每行由x_y组成（例如：223333_14512，放了demo），代表挑出来的合法坐标

## 2.依赖

- Linux or MacOS

- pysam>=0.18.0   解析bam文件
- tqdm>=4.62.3    进度条
- numpy>=1.20.3   老朋友，不解释
- pandas>=1.3.5   老朋友，不解释
- anndata>=0.7.8  用来生成anndata
- scipy>=1.7.3   用来操作矩阵
- scvelo>=0.2.4  RNA速率分析
- scanpy>=1.8.2  单细胞分析整合包

> 注意：由于时间有限，并未测试低版本的python package，按需安装即可
>
> 由于未创建虚拟环境，并未freeze成requirement

## 3.使用帮助

使用顺序：

1. python3 bam2rawgem.py -h
2. python3 gem_format.py -h
3. python3 gem2adata.py -h
4. analysis.py (自行分析)

使用帮助：自行查看help