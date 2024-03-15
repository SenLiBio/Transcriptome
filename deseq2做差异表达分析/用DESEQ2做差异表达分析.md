# 用DESEQ2做差异表达基因分析

# 0. 用到的软件

| 软件 | 版本 |
| --- | --- |
| fastp | 0.20.1 |
| hisat | 2.2.1 |
| samtools | 1.12 |
| stringtie | v2.1.6 |
| gffread | 0.12.7 |
| gffcompare | 0.12.6 |

# 1. RNA测序的序列过滤和比对

```bash
fastp -i Prefix.R1.fastq.gz -I Prefix.R2.fastq.gz -o Prefix.R1.clean.fastq.gz -O Prefix.R2.clean.fastq.gz -A -q 20 -u 10 -n 5 -w 4 -j summary/Prefix.json -h summary/Prefix.html
hisat2-build Ref.fa Rer
hisat2 --phred33 -p 36 Ref -1 Prefix.R1.clean.fastq.gz -2 Prefix.R2.clean.fastq.gz -S Prefix.sam --summary-file Prefix.txt
# -p, threads
samtools sort -@ 36 -o Prefix.sorted.bam Prefix.sam
# -@, threads, sometimes it dosen't work well.
```

---

# 2. stringtie做转录本组装和定量

<aside>
💡 可以用gffread将gff转为gtf
`conda install -c bioconda gffread
gffread Bgy.gff -T -o bgy.gtf`

</aside>

## 2.1 对每个个体的转录本进行组装

<aside>
💡 如果不做2.2的内容，就加上-e，如果继续做，就不加。此外不加-e产生的结果不能用于后续prepDE.py的分析

</aside>

```bash
stringtie -e -p 16 -G Ref.gtf -l Ref -o Prefix.gtf Prefix.sorted.bam
# -p, threads
# -G, reference annotation to use for guiding the assembly process
# -l, name prefix for output transcripts
# -o, output path/file name for the assembled transcripts GTF
# -e, this option directs StringTie to operate in expression estimation mode; this limits the processing of read alignments to estimating the coverage of the transcripts given with the -G option (hence this option requires -G).
```

## 2.2 检测转录本检测的精确性和完整性（选做）

<aside>
💡 这一步主要是在无参的时候做，或者说要预测一些可变剪切做后续分析，由于stringtie主要是根据reads比对到的位置来预测转录本，所以有可能同一个基因会产生多个可变剪切，这种情况下有可能反而对后续分析有影响。做基因组的分析的时候可以不做2.2的内容，一个有参考价值的[博客](https://www.jianshu.com/p/40c1c0606db6?utm_medium=timeline)和[github issue](https://github.com/gpertea/stringtie/issues/192)。

</aside>

### 2.2.1 合并转录本

```bash
stringtie --merge -p 16 -G Ref.gtf -o Ref.merged.gtf -l Ref Ref.mergelist.txt
# -p, threads
# -G, reference annotation to include in the merging (GTF/GFF3)
# -o, output file name for the merged transcripts GTF
# -l, name prefix for output transcripts
```

```bash
# Ref.mergelist.txt
# Each line contains the name of one of the gtfs to merge.
SRR15464825.gtf
SRR15464826.gtf
SRR15464827.gtf
SRR15464828.gtf
SRR15464830.gtf
SRR15464831.gtf
SRR15497174.gtf
```

### 2.2.2 比较参考转录本和合并的转录本

通过这一步来检验转录本组装的精确度和准确度。具体的检验结果说明可以参考[简书的这篇说明](https://www.jianshu.com/p/08bb35f8e7c2)和[官网的documentary](http://ccb.jhu.edu/software/stringtie/gffcompare.shtml)。

```bash
gffcompare -r Ref.gtf -G -o Ref Ref.merged.gtf
# -r, reference annotation file (GTF/GFF)
# -G, 比较输入的gtf中所有的转录本，即使它们有可能是冗余的
# -o, prefix of outputs
```

### 2.2.3 重新组装转录本

如果说发现合并的转录本相较于原本的参考转录本（Ref.gtf）文件有较多新发现的转录本，则可以用合并的转录本作为reference重新组装转录本，否则直接用之前组装的结果即可。

```bash
stringtie -e -p 16 -G Ref.merged.gtf -l Ref -A Prefix.tab -o Prefix.new.gtf Prefix.bam

```

---

# 3. 准备deseq2所需的文件

stringtie中提供了脚本，可以直接使用，首先要把所有的gtf按照名字放在对应的文件夹。

![截屏2021-12-16 下午2.48.13.png](%E7%94%A8DESEQ2%E5%81%9A%E5%B7%AE%E5%BC%82%2029da4/%E6%88%AA%E5%B1%8F2021-12-16_%E4%B8%8B%E5%8D%882.48.13.png)

```bash
prepDE.py -i ./newgtfs/ -g Ref.gene.matrix -t Ref.transcript.matrix
# -i, a folder containing all sample sub-directories, or a text file with sample               ID and path to its GTF file on each line
# -g, where to output the gene count matrix
# -t, where to output the transcript count matrix (如果没有做2.2，这两个文件其实是一样的）
```

---

# 3. deseq2分析

```r
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library("DESeq2")

# 新版的r里面已经用biomanager替代了bioclite
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
library("DESeq2")

countData <- as.matrix(read.csv("Ref.gene.matrix ", row.names="gene_id"))
# Ref.gene.matrix是prepDE.py生成的文件

colData <- read.csv("PhenoData.txt", sep="\t", row.names=1)
# 记得检查分隔符

all(rownames(colData) %in% colnames(countData))
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))
# 检查coldata和countdata的sample id是否一一对应（TURE)

dds <- DESeqDataSetFromMatrix(countData = countData,colData = colData, design = ~ CHOOSE_FEATURE)
# 记得看design的表达式，以colData里面的列的名字做组合。

dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]
resOrdered=as.data.frame(resOrdered)
write.table(resOrdered,"DEGoutput.txt",sep="\t")
 
```

---

# 4. 对WGD产生的paralog做表达量差异分析

## 4.1 提取每个基因的TPM值

<aside>
💡 做不做2.2这一步的脚本是不一样的，这里的是不做2.2的版本，如果做了2.2，不同的transcript可能对应同一个geneid，这一点要注意。

</aside>

```bash
perl extractTPM.pl >Ref.TPM.matrix
```

- extractTPM.pl
    
    ```perl
    use strict;
    use warnings;
    use autodie;
    use 5.010;
    
    my %TPM;
    my $NumOfFile=-1;
    my @individual_names;
    my @GeneIds;
    my @gtfs=glob ("*gtf");
    foreach my $gtf (@gtfs) {
            $NumOfFile++;
            (my $individual_name=$gtf)=~s/\.gtf//;
            push @individual_names,$individual_name;
            open GTF,"<","$gtf";
            while (<GTF>) {
                next if /^#/;
                next if (/exon/);
    	        	$_=~m/gene_id "(?<GeneId>.*?)";.*TPM "(?<TPM>.*)";/;
                my $GeneId=$+{GeneId};
    						my $TPM=sprintf "%.2f",$+{TPM};
    						$TPM{$GeneId}[$NumOfFile]=$TPM;
                push @GeneIds,$GeneId;
            }
            close GTF;
    }
    
    my $individual_names=join(',',@individual_names);
    say "gene_id,$individual_names";
    my %bak;
    map {$bak{$_}++} @GeneIds;
    @GeneIds=keys %bak;
    foreach my $GeneId (@GeneIds) {
            my $exp=join(',',$GeneId,@{$TPM{$GeneId}});
            say $exp;
    }
    ```
    

## 4.2 做Wilcoxon成对秩和检验

### 4.2.1 检验

```perl
perl wilcoxon.pair.pl [TPMMatrixFile] [KsMedianFile] [T0File] [CollinearityFile]
```

- wilcoxon.pair.pl
    
    ```perl
    ###     本脚本用来做wilcoxon成对秩和检验
    ###     T0值通过查表根据样本大小获得
    ###     记得修改39行ks中位数的区间
    ###     记得修改第56行输出文件的名字
    
    use 5.012;
    use strict;
    use warnings;
    use autodie;
    
    die "Usage: perl $0 [TPMMatrixFile] [KsMedianFile] [T0File] [CollinearityFile]\n" if @ARGV<4;
    
    my $TPMMatrixFile=shift;
    my $KsMedianFile=shift;
    my $T0File=shift;
    my $CollinearFile=shift;
    
    # TPM data
    open TPM,"<","$TPMMatrixFile";
    my %Exp;
    while (<TPM>) {
            next if /^gene/;
    		    s/evm\.model\.Bruguiera_gymnorrhiza_PacBio_hic_scaffold_/bgy\./g;
            s/\s+$//;
            my @eles=split/,/;
            my $GeneId=shift @eles;
            $Exp{$GeneId}=[@eles];
    }
    close TPM;
    
    # 筛选近期wgd产生的block
    open MEDIAN,"<","$KsMedianFile";
    my @RecentWGDBlocks;
    while (<MEDIAN>) {
                    s/\s+$//;
                    my ($block,$ks)=split/\s+/;
                    $block=~m/(\d+)$/;
                    $block=$1;
                    push @RecentWGDBlocks,$block if ($ks<=0.8 && $ks>=0.3);
    }
    close MEDIAN;
    
    # T0的阈值
    open TVALUE,"<","$T0File";
    my %TValue;
    while (<TVALUE>) {
    	s/\s+$//;
    	my ($n,$t0)=split/\s+/;
    	$TValue{$n}=$t0;
    }
    close TVALUE;
    
    # Extract WGD retention gene pairs.
    open COL,"<","$CollinearFile";
    open OUT,">","TPM.wilcoxon.txt";
    say OUT "# lt: lower than\n# gt: greater than\n#ns: no significance";
    L: while (<COL>) {
            next if /^#/;
            $_=~m/^\s?(\d+)\-/;
            my $block=$1;
            next if (($block~~@RecentWGDBlocks)==0);
            $_=~m/(bgy.*?)\s+(bgy.*?)\s+/;
            my $Gene1=$1;
            my $Gene2=$2;
            next L if (!exists $Exp{$Gene1});
            next L if (!exists $Exp{$Gene2});
            WilcoxonPairedTest ($Gene1,$Gene2);
    }
    close COL;
    close OUT;
    
    sub WilcoxonPairedTest {
            my $Gene1=shift;
            my $Gene2=shift;
            my @Gene1Exp=@{$Exp{$Gene1}};
            my @Gene2Exp=@{$Exp{$Gene2}};
            my $Gene1Exp=join(',',@{$Exp{$Gene1}});
            my $Gene2Exp=join(',',@{$Exp{$Gene2}});        
    print "$Gene1\n$Exp{$Gene1}[0]\n$Gene1Exp[0]\n$Gene2\n$Exp{$Gene2}[0]\n$Gene2Exp[0]\n";
            my (@Differences,@ABSDifferences);
            my ($Rank,@RankOrder);
    		my $PositiveRankSum=0;
    		my $NegativeRankSum=0;
    # 求表达量的差值
            for (1..@Gene1Exp) {
                    my $Gene1Exp=shift @Gene1Exp;
                    my $Gene2Exp=shift @Gene2Exp;
                    my $Difference=$Gene1Exp-$Gene2Exp;
                    push @Differences,$Difference;
            }
    # 排序并记录绝对值
            @Differences=sort {abs($a)<=>abs($b)} @Differences;
            @ABSDifferences=map {abs ($_)} @Differences;
            my %num;
            map {$num{$_}++} @ABSDifferences;
    # 差值为0的话在排序的时候要删除
            delete $num{'0'};
    		my @NonzeroValues=keys %num;
    		my $NonzeroValueNum=@NonzeroValues;
    		next L if $NonzeroValueNum<6;
            $Rank=0;
            @RankOrder;
            @ABSDifferences=keys %num;
    # 根据差值的绝对值求秩
            foreach my $ABSDifference (@ABSDifferences) {
                    if ($num{$ABSDifference}==1) {
                            $Rank+=1;
                            push @RankOrder,$Rank;
                    } else {
                            my $RankLow=$Rank+1;
                            my $RankHigh=$Rank+$num{$ABSDifference};
                            my $AverageRank=($RankLow+$RankHigh)/2;
                            for (1..$num{$ABSDifference}) {push @RankOrder,$AverageRank};
                            $Rank=$RankHigh;
                    }
            }
            foreach my $Difference (@Differences) {
                    next if ($Difference==0);
                    my $Rank=shift @RankOrder;
                    $PositiveRankSum+=$Rank if $Difference>0;
                    $NegativeRankSum+=$Rank if $Difference<0;
            }
            my $LowerRankSum=($PositiveRankSum>$NegativeRankSum)?($NegativeRankSum):($PositiveRankSum);
            if ($LowerRankSum<=$TValue{$NonzeroValueNum}) {
                    my $RankSumStat=($PositiveRankSum>$NegativeRankSum)?('gt'):('lt');
                    say OUT "$Gene1\t$Gene1Exp\t$Gene2\t$Gene2Exp\t$RankSumStat";
            } else {
                    say OUT "$Gene1\t$Gene1Exp\t$Gene2\t$Gene2Exp\tns";
            }
    }
    ```
    

### 4.2.2 画图

<aside>
💡 先把上一步输出的结果整理成`TPM.matrix.forheatmap.txt` 形式，然后有一个文件指定列的分组，如SampleGroup.txt，之后用TBtools画图即可，注意可以分一个maingroup和一个subgroup。

</aside>

- TPM.matrix.forheatmap.txt
    
    ![Untitled](%E7%94%A8DESEQ2%E5%81%9A%E5%B7%AE%E5%BC%82%2029da4/Untitled.png)
    
- SampleGroup.txt
    
    ![Untitled](%E7%94%A8DESEQ2%E5%81%9A%E5%B7%AE%E5%BC%82%2029da4/Untitled%201.png)
    

## 4.3 整理wilcoxon秩和检验和差异表达基因的结果

<aside>
💡 注意文件中基因名要一一对应，或者更改第一个循环中注释掉的那行中的正则。

</aside>

### 4.3.1 看秩和检验中有表达差异的基因对在盐胁迫下是否有差异表达

```bash
perl wilcoxon_deg.pl DEGoutput.txt TPM.matrix.forheatmap.txt
```

- wilcoxon_deg.pl
    
    ```perl
    open DEG,"<","$ARGV[0]";
    open WGD,"<","$ARGV[1]";
    my %ExpDifference;
    
    while (<DEG>) {
        next if /baseMean/;
    # s/evm\.model\.Bruguiera_gymnorrhiza_PacBio_hic_scaffold_/bgy\./;
        my ($genename,$log2fc)=(split/\s+/)[0,2];
    		$genename=~s/"//g;
        $ExpDifference{$genename}='Up' if $log2fc<-1;
        $ExpDifference{$genename}='Down' if $log2fc>1;
    }
    close DEG;
    
    open OUT,">DEG.wilcoxon.txt";
    
    while (<WGD>) {
        next if /GeneID/;
        my $genepair=(split/\s+/)[0];
        my ($gene1,$gene2)=split/_/,$genepair;
        print OUT  "$gene1\t$gene2\t";
        map {
            if (exists $ExpDifference{$_}) {print OUT "$ExpDifference{$_}\t"} else {print OUT "NS\t"};
        } ($gene1,$gene2);
            print OUT "\n";
    }
    close WGD;
    close OUT;
    ```
    

### 4.3.2 去除均不差异表达的基因对

```bash
perl remove.doubleNS.pl DEG.wilcoxon.txt >DEG.wilcoxon.filtered.txt
```

- remove.doubleNS.pl
    
    ```perl
    while (<>) {
    	my $line=$_;
    	s/\s+$//;
    	my ($sig1,$sig2)=(split/\s+/)[2,3];
    	next if ($sig1 eq 'NS' && $sig2 eq 'NS');
    	print $line
    }
    ```