# Demo and Examples

In the testing below, we assume that
1. You have installed DeepRepeat and its dependent packages. If not, please refer to `README.md` to install them
2. Your current folder have the DeepRepeat folder which have been well compiled according the last step.
3. The following folder organization is expected:
```
	./
		DeepRepeat/
		trainedmod_32_0.2/
		test_op/
		htt12/
		na12878_loci/
		strique_g4c2_data/
```

## Reproduce results for 9 loci and other loci

### For 9 loci
```
wget na12878_loci_nanopore.tar.gz
tar -xvf na12878_loci_nanopore.tar.gz

```

You should have a subfolder `na12878_loci/` which contains 9 sub-folders for each of 9 loci. If no, please check your downloading or check the folder organization to make it correct. After that, the following commands can be used to generate the prediction for 9 loci on NA12878.

For chr16:73546662-73546736 with TGC repeats
```
python DeepRepeat/bin/DeepRepeat.py Detect --gn hx1 --TempRem 0 --epchon 200 --repeat_relax_bp 20 --UniqueID TGC_chr16_73546662-73546736 --is_pcr 0 --repeatName TGC_chr16_73546662-73546736 --repeat chr16:73546662-73546736:TGC:3 --f5i na12878_loci/TGC_chr16_73546662-73546736/na.f5index --o test_op/TGC_chr16_73546662-73546736 --bam na12878_loci/TGC_chr16_73546662-73546736/TGC_chr16_73546662-73546736.bam --f5folder na12878_loci/TGC_chr16_73546662-73546736 --algLen 500 2>&1 | grep "GMM:"
```
For chr6:145272942-145273029 with AGAT repeats
```
python DeepRepeat/bin/DeepRepeat.py Detect --gn hx1 --TempRem 0 --epchon 200 --repeat_relax_bp 20 --UniqueID AGAT_chr6_145272942-145273029 --is_pcr 0 --repeatName AGAT_chr6_145272942-145273029 --repeat chr6:145272942-145273029:AGAT:4 --f5i na12878_loci/AGAT_chr6_145272942-145273029/na.f5index --o test_op/AGAT_chr6_145272942-145273029 --bam na12878_loci/AGAT_chr6_145272942-145273029/AGAT_chr6_145272942-145273029.bam --f5folder na12878_loci/AGAT_chr6_145272942-145273029 --algLen 500 2>&1 | grep "GMM:"
```
For chrX:55366379-55366473 with TATC repeats
```
python DeepRepeat/bin/DeepRepeat.py Detect --gn hx1 --TempRem 0 --epchon 200 --repeat_relax_bp 20 --UniqueID TATC_chrX_55366379-55366473 --is_pcr 0 --repeatName TATC_chrX_55366379-55366473 --repeat chrX:55366379-55366473:TATC:4 --f5i na12878_loci/TATC_chrX_55366379-55366473/na.f5index --o test_op/TATC_chrX_55366379-55366473 --bam na12878_loci/TATC_chrX_55366379-55366473/TATC_chrX_55366379-55366473.bam --f5folder na12878_loci/TATC_chrX_55366379-55366473 --algLen 500 2>&1 | grep "GMM:"
```
For chr10:107609688-107609784 with ATCT repeats
```
python DeepRepeat/bin/DeepRepeat.py Detect --gn hx1 --TempRem 0 --epchon 200 --repeat_relax_bp 20 --UniqueID ATCT_chr10_107609688-107609784 --is_pcr 0 --repeatName ATCT_chr10_107609688-107609784 --repeat chr10:107609688-107609784:ATCT:4 --f5i na12878_loci/ATCT_chr10_107609688-107609784/na.f5index --o test_op/ATCT_chr10_107609688-107609784 --bam na12878_loci/ATCT_chr10_107609688-107609784/ATCT_chr10_107609688-107609784.bam --f5folder na12878_loci/ATCT_chr10_107609688-107609784 --algLen 500 2>&1 | grep "GMM:"
```
For chr7:152384548-152384614 with CAC repeats
```
python DeepRepeat/bin/DeepRepeat.py Detect --gn hx1 --TempRem 0 --epchon 200 --repeat_relax_bp 20 --UniqueID CAC_chr7_152384548-152384614 --is_pcr 0 --repeatName CAC_chr7_152384548-152384614 --repeat chr7:152384548-152384614:CAC:3 --f5i na12878_loci/CAC_chr7_152384548-152384614/na.f5index --o test_op/CAC_chr7_152384548-152384614 --bam na12878_loci/CAC_chr7_152384548-152384614/CAC_chr7_152384548-152384614.bam --f5folder na12878_loci/CAC_chr7_152384548-152384614 --algLen 500 2>&1 | grep "GMM:"
```
For chr4:42555085-42555198 with TATC repeats
```
python DeepRepeat/bin/DeepRepeat.py Detect --gn hx1 --TempRem 0 --epchon 200 --repeat_relax_bp 20 --UniqueID TATC_chr4_42555085-42555198 --is_pcr 0 --repeatName TATC_chr4_42555085-42555198 --repeat chr4:42555085-42555198:TATC:4 --f5i na12878_loci/TATC_chr4_42555085-42555198/na.f5index --o test_op/TATC_chr4_42555085-42555198 --bam na12878_loci/TATC_chr4_42555085-42555198/TATC_chr4_42555085-42555198.bam --f5folder na12878_loci/TATC_chr4_42555085-42555198 --algLen 500 2>&1 | grep "GMM:"
```
For chr1:161051967-161052060 with CAC repeats
```
python DeepRepeat/bin/DeepRepeat.py Detect --gn hx1 --TempRem 0 --epchon 200 --repeat_relax_bp 20 --UniqueID CAC_chr1_161051967-161052060 --is_pcr 0 --repeatName CAC_chr1_161051967-161052060 --repeat chr1:161051967-161052060:CAC:3 --f5i na12878_loci/CAC_chr1_161051967-161052060/na.f5index --o test_op/CAC_chr1_161051967-161052060 --bam na12878_loci/CAC_chr1_161051967-161052060/CAC_chr1_161051967-161052060.bam --f5folder na12878_loci/CAC_chr1_161051967-161052060 --algLen 500 2>&1 | grep "GMM:"
```
For chr12:4702128-4702202 with CCA repeats
```
python DeepRepeat/bin/DeepRepeat.py Detect --gn hx1 --TempRem 0 --epchon 200 --repeat_relax_bp 20 --UniqueID CCA_chr12_4702128-4702202 --is_pcr 0 --repeatName CCA_chr12_4702128-4702202 --repeat chr12:4702128-4702202:CCA:3 --f5i na12878_loci/CCA_chr12_4702128-4702202/na.f5index --o test_op/CCA_chr12_4702128-4702202 --bam na12878_loci/CCA_chr12_4702128-4702202/CCA_chr12_4702128-4702202.bam --f5folder na12878_loci/CCA_chr12_4702128-4702202 --algLen 500 2>&1 | grep "GMM:"
```
For chr17:10995712-10995775 with CTG repeats
```
python DeepRepeat/bin/DeepRepeat.py Detect --gn hx1 --TempRem 0 --epchon 200 --repeat_relax_bp 20 --UniqueID CTG_chr17_10995712-10995775 --is_pcr 0 --repeatName CTG_chr17_10995712-10995775 --repeat chr17:10995712-10995775:CTG:3 --f5i na12878_loci/CTG_chr17_10995712-10995775/na.f5index --o test_op/CTG_chr17_10995712-10995775 --bam na12878_loci/CTG_chr17_10995712-10995775/CTG_chr17_10995712-10995775.bam --f5folder na12878_loci/CTG_chr17_10995712-10995775 --algLen 500 2>&1 | grep "GMM:"
```
Here, the options `--TempRem`, `--epchon` and `--is_pcr` can be ignored since they have default parameters.
Usually, you only need several minutes to finish each of the running. The dictionary after `allocr:` provides more detail, and you might need to check it if the number of supporting reads is smaller or the repeat count distribution is complicated. The format of the dictionary after `allocr:` is: XXX1:YYY2, XXX1:YYY2, .... XXXn:YYYn where XXX is the repeat count detected from a long read and YYY is the number of long reads with this repeat count. One or Two peaks are expected for human genome.


### For any other loci
If you want to run other loci on Na12878 or other genome, you need to (1) download the fast5 data and basecall them with Albacore v2.3*, and (2) generated alignment BAM files. Then, 
```
python DeepRepeat/bin/DeepRepeat.py Detect --gn "hx1-na12878" --TempRem 0 --epchon 200 --repeat_relax_bp 20 --UniqueID UserDefinedID --is_pcr 0 --repeatName UserDefinedID --repeat chr:start_pos-end_pos:Repeat_motif:Repeat_len --o test_op/ --bam UserAligned.bam --f5folder Basecalled_f5_folder --algLen 500 2>&1 
```


## Reproduce results for 12 HTT samples
The datasets of Nanopore long-read data for 12 HTT samples have been uploaded to SRA(PRJNA678742. The data will be available once the paper is accepted.). You need to download them or one of them firstly before you run the commands below. Similarly, I assume that you have a sub-folder `htt12` which has a sub-folder for one of 12 HTT samples. Then, the command below can be run to get the results. Please note that you should already have `htt.f5index` in each subfolder of `htt12/`. if no, please download them using `wget htt12.f5index.tar.gz` and `tar -xvf htt12.f5index.tar.gz`. **Without `htt.f5index`, you cannot run successfully DeepRepeat.**

Firstly, make a log folder.
```
mkdir log
```
For ND30047,
```
time python DeepRepeat/bin/DeepRepeat.py Detect --gn CAG_htt2_90 --TempRem 0 --epchon 50 --repeat_relax_bp 15 --UniqueID ND30047 --is_pcr 0 --repeatName ND30047_htt --repeat chr4:3074876-3074933:CAG:3 --f5i htt12/ND30047/htt.f5index --o test_op/ND30047 --bam htt12/ND30047/ND30047.bam --f5folder htt12/ND30047 --algLen 100 --rpg DeepRepeat/bin/data/wg_trf.cag.bed  > log/ND30047_dr012.log 2>&1
grep "GMM:" log/ND30047_dr012.log 
```
For ND40534, 
```
time python DeepRepeat/bin/DeepRepeat.py Detect --gn CAG_htt2_90 --TempRem 0 --epchon 50 --repeat_relax_bp 15 --UniqueID ND40534 --is_pcr 0 --repeatName ND40534_htt --repeat chr4:3074876-3074933:CAG:3 --f5i htt12/ND40534/htt.f5index --o test_op/ND40534 --bam htt12/ND40534/ND40534.bam --f5folder htt12/ND40534 --algLen 100 --rpg DeepRepeat/bin/data/wg_trf.cag.bed  > log/ND40534_dr012.log 2>&1
grep "GMM:" log/ND40534_dr012.log
```
For ND30015,
```
time python DeepRepeat/bin/DeepRepeat.py Detect --gn CAG_htt2_90 --TempRem 0 --epchon 50 --repeat_relax_bp 15 --UniqueID ND30015 --is_pcr 0 --repeatName ND30015_htt --repeat chr4:3074876-3074933:CAG:3 --f5i htt12/ND30015/htt.f5index --o test_op/ND30015 --bam htt12/ND30015/ND30015.bam --f5folder htt12/ND30015 --algLen 100 --rpg DeepRepeat/bin/data/wg_trf.cag.bed  > log/ND30015_dr012.log 2>&1
grep "GMM:" log/ND30015_dr012.log
```
For ND31551,
```
time python DeepRepeat/bin/DeepRepeat.py Detect --gn CAG_htt2_90 --TempRem 0 --epchon 50 --repeat_relax_bp 15 --UniqueID ND31551 --is_pcr 0 --repeatName ND31551_htt --repeat chr4:3074876-3074933:CAG:3 --f5i htt12/ND31551/htt.f5index --o test_op/ND31551 --bam htt12/ND31551/ND31551.bam --f5folder htt12/ND31551 --algLen 100 --rpg DeepRepeat/bin/data/wg_trf.cag.bed  > log/ND31551_dr012.log 2>&1
grep "GMM:" log/ND31551_dr012.log
```
For NA12878,
```
time python DeepRepeat/bin/DeepRepeat.py Detect --gn CAG_htt2_90 --TempRem 0 --epchon 50 --repeat_relax_bp 15 --UniqueID NA12878 --is_pcr 0 --repeatName NA12878_htt --repeat chr4:3074876-3074933:CAG:3 --f5i htt12/NA12878/htt.f5index --o test_op/NA12878 --bam htt12/NA12878/NA12878.bam --f5folder htt12/NA12878 --algLen 100 --rpg DeepRepeat/bin/data/wg_trf.cag.bed  > log/NA12878_dr012.log 2>&1
grep "GMM:" log/NA12878_dr012.log
```
For ND33947,
```
time python DeepRepeat/bin/DeepRepeat.py Detect --gn CAG_htt2_90 --TempRem 0 --epchon 50 --repeat_relax_bp 15 --UniqueID ND33947 --is_pcr 0 --repeatName ND33947_htt --repeat chr4:3074876-3074933:CAG:3 --f5i htt12/ND33947/htt.f5index --o test_op/ND33947 --bam htt12/ND33947/ND33947.bam --f5folder htt12/ND33947 --algLen 100 --rpg DeepRepeat/bin/data/wg_trf.cag.bed  > log/ND33947_dr012.log 2>&1
grep "GMM:" log/ND33947_dr012.log
```
For ND33392,
```
time python DeepRepeat/bin/DeepRepeat.py Detect --gn CAG_htt2_90 --TempRem 0 --epchon 50 --repeat_relax_bp 15 --UniqueID ND33392 --is_pcr 0 --repeatName ND33392_htt --repeat chr4:3074876-3074933:CAG:3 --f5i htt12/ND33392/htt.f5index --o test_op/ND33392 --bam htt12/ND33392/ND33392.bam --f5folder htt12/ND33392 --algLen 100 --rpg DeepRepeat/bin/data/wg_trf.cag.bed  > log/ND33392_dr012.log 2>&1
grep "GMM:" log/ND33392_dr012.log
```
For ND30422,
```
time python DeepRepeat/bin/DeepRepeat.py Detect --gn CAG_htt2_90 --TempRem 0 --epchon 50 --repeat_relax_bp 15 --UniqueID ND30422 --is_pcr 0 --repeatName ND30422_htt --repeat chr4:3074876-3074933:CAG:3 --f5i htt12/ND30422/htt.f5index --o test_op/ND30422 --bam htt12/ND30422/ND30422.bam --f5folder htt12/ND30422 --algLen 100 --rpg DeepRepeat/bin/data/wg_trf.cag.bed  > log/ND30422_dr012.log 2>&1
grep "GMM:" log/ND30422_dr012.log
```
For ND30626
```
time python DeepRepeat/bin/DeepRepeat.py Detect --gn CAG_htt2_90 --TempRem 0 --epchon 50 --repeat_relax_bp 15 --UniqueID ND30626 --is_pcr 0 --repeatName ND30626_htt --repeat chr4:3074876-3074933:CAG:3 --f5i htt12/ND30626/htt.f5index --o test_op/ND30626 --bam htt12/ND30626/ND30626.bam --f5folder htt12/ND30626 --algLen 100 --rpg DeepRepeat/bin/data/wg_trf.cag.bed  > log/ND30626_dr012.log 2>&1
grep "GMM:" log/ND30626_dr012.log
```
For ND30016
```
time python DeepRepeat/bin/DeepRepeat.py Detect --gn CAG_htt2_90 --TempRem 0 --epchon 50 --repeat_relax_bp 15 --UniqueID ND30016 --is_pcr 0 --repeatName ND30016_htt --repeat chr4:3074876-3074933:CAG:3 --f5i htt12/ND30016/htt.f5index --o test_op/ND30016 --bam htt12/ND30016/ND30016.bam --f5folder htt12/ND30016 --algLen 100 --rpg DeepRepeat/bin/data/wg_trf.cag.bed  > log/ND30016_dr012.log 2>&1
grep "GMM:" log/ND30016_dr012.log
```
For ND29970
```
time python DeepRepeat/bin/DeepRepeat.py Detect --gn CAG_htt2_90 --TempRem 0 --epchon 50 --repeat_relax_bp 15 --UniqueID ND29970 --is_pcr 0 --repeatName ND29970_htt --repeat chr4:3074876-3074933:CAG:3 --f5i htt12/ND29970/htt.f5index --o test_op/ND29970 --bam htt12/ND29970/ND29970.bam --f5folder htt12/ND29970 --algLen 100 --rpg DeepRepeat/bin/data/wg_trf.cag.bed  > log/ND29970_dr012.log 2>&1
grep "GMM:" log/ND29970_dr012.log
```
For GM04723
```
time python DeepRepeat/bin/DeepRepeat.py Detect --gn CAG_htt2_90 --TempRem 0 --epchon 50 --repeat_relax_bp 15 --UniqueID GM04723 --is_pcr 0 --repeatName GM04723_htt --repeat chr4:3074876-3074933:CAG:3 --f5i htt12/GM04723/htt.f5index --o test_op/GM04723 --bam htt12/GM04723/GM04723.bam --f5folder htt12/GM04723 --algLen 100 --rpg DeepRepeat/bin/data/wg_trf.cag.bed  > log/GM04723_dr012.log 2>&1
grep "GMM:" log/GM04723_dr012.log
```
Here, the options `--TempRem`, `--is_pcr` or  `--repeat_relax_bp` can be ignored since they have default parameters.
Usually, you only need several minutes to finish each of the running. The dictionary after `allocr:` provides more detail, and you might need to check it if the number of supporting reads is smaller or the repeat count distribution is complicated. The format of the dictionary after `allocr:` is: XXX1:YYY2, XXX1:YYY2, .... XXXn:YYYn where XXX is the repeat count detected from a long read and YYY is the number of long reads with this repeat count. One or Two peaks are expected for human genome.



## Reproduce results for G4C2 results

The steps below are used to generate results of G4C2 results. Please note that you need to install Albacore v2.3.* for basecalling of the data.
### 1. Preparation
```
mkdir strique_g4c2_data/
```
Then, download G4C2 data via `wget https://ndownloader.figshare.com/files/15298166`, and then use Albacore for basecalling. I recommended storing the basecalled data under `strique_g4c2_data/alb23/. Feel free to change the location of basecalled data but you need to change the location path below accordingly.

### 2. Basecalling and alignment
#### 2.1 Get fq files
```
mkdir -p strique_g4c2_data/m_fq
python DeepRepeat/tools/get_fq_barcode.py strique_g4c2_data/alb23/workspace/pass/ strique_g4c2_data/m_fq
```

#### 2.2 Prepare reference file
Download reference sequence for G4C2 samples and stored them under `strique_g4c2_data/pCRAmpBE.fa`

#### 2.3 Alingment
```
mkdir -p strique_g4c2_data/bam_results/
python DeepRepeat/tools/get_bam.py strique_g4c2_data/m_fq/ strique_g4c2_data/bam_results/ ".fq" strique_g4c2 strique_g4c2_data/pCRAmpBE.fa 
```

### 3. Run DeepRepeat
For barcode 10,
```
python DeepRepeat/bin/DeepRepeat.py Detect --gn CCCCGG_strique_g4c2_3 --TempRem 0 --epchon 500 --repeat_relax_bp 50 --UniqueID barcode10 --is_pcr 0 --repeatName barcode10 --repeat pCRAmpBE:1107-1135:GGGGCC:6 --o test_op/barcode10 --bam strique_g4c2_data/bam_results/barcode10.bam --f5folder strique_g4c2_data/alb23/ --algLen 100 --rpg strique_g4c2_data/plasmid_trf.bed  > log/barcode10.log 2>&1
grep "GMM:" log/barcode10.log
```
For barcode 11,
```
python DeepRepeat/bin/DeepRepeat.py Detect --gn CCCCGG_strique_g4c2_3 --TempRem 0 --epchon 500 --repeat_relax_bp 50 --UniqueID barcode11 --is_pcr 0 --repeatName barcode11 --repeat pCRAmpBE:1107-1135:GGGGCC:6 --o test_op/barcode11 --bam strique_g4c2_data/bam_results/barcode11.bam --f5folder strique_g4c2_data/alb23/ --algLen 100 --rpg strique_g4c2_data/plasmid_trf.bed  > log/barcode11.log 2>&1
grep "GMM:" log/barcode11.log
```
The format after `allocr:` is similar to the description above. This output have no clear peak, due to (1) the demultiplexing issues and thus many peaks are detected, and (2) the training data is smaller.




